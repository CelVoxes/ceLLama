# Load required library
library(httr)

#' Select Top Up- and Down-Regulated Genes
#'
#' @param degs A data frame containing gene markers with their log2 fold changes.
#' @param n Number of top genes to select for both up- and down-regulated lists.
#' @return A list with top up-regulated and down-regulated genes.
#'
#' @export
select_top_genes <- function(degs, n = 20) {
  up_genes <- degs[order(degs$avg_log2FC, decreasing = TRUE), "gene"][1:n]
  down_genes <- degs[order(degs$avg_log2FC, decreasing = FALSE), "gene"][1:n]
  list(up = up_genes, down = down_genes)
}

#' Format Annotation Data
#'
#' @param top_genes A list of top genes for each cluster.
#' @param base_prompt A base prompt to use for formatting the annotation data.
#' @return A list of formatted annotation prompts for each cluster.
#'
#' @export
format_annotation_data <- function(top_genes, base_prompt) {
  annotation_data <- lapply(names(top_genes), function(cluster) {
    up_genes <- paste(top_genes[[cluster]]$up, collapse = ", ")
    down_genes <- paste(top_genes[[cluster]]$down, collapse = ", ")
    prompt <- paste(
      "This cell cluster (", cluster, ") has up-regulated genes:", up_genes,
      "and down-regulated genes:", down_genes, "."
    )
    return(prompt)
  })
  return(annotation_data)
}

#' Get Annotation for a Cluster
#'
#' @param description The description for the cluster to be annotated.
#' @param model The model to use for annotation.
#' @param url The URL for the annotation service.
#' @return The annotation for the cluster.
#'
#' @export
get_annotation <- function(description, model, url, seed = 100, temperature = 0) {
  data <- list(
    model = model,
    stream = FALSE,
    prompt = description,

    # reproducible results; if seed given and temp 0.
    options = list(
      seed = seed,
      temperature = temperature
    )
  )
  response <- httr::POST(url, body = data, encode = "json")
  content <- httr::content(response, "parsed")
  return(content$response)
}

#' Get Reason for the Annotation
#'
#' @param annotation The annotation received for the cluster.
#' @param description The description used to generate the annotation.
#' @param model The model to use for the reason generation.
#' @param url The URL for the reason generation service.
#' @return The reason for the annotation.
#'
#' @export
get_reason <- function(annotation, description, model, url, seed = 100, temperature = 0, verbose = T) {
  reason_prompt <- paste(
    "The annotation for the cell cluster is:", annotation,
    ". Can you provide the reason for this annotation based on the following description:",
    description
  )
  data <- list(
    model = model,
    stream = FALSE,
    prompt = reason_prompt,

    # reproducible results; if seed given and temp 0.
    options = list(
      seed = seed,
      temperature = temperature
    )
  )
  response <- httr::POST(url, body = data, encode = "json")
  content <- httr::content(response, "parsed")
  if (verbose) message(">> Reason: ", content$response)

  return(content$response)
}

get_openai_annotation <- function(description, model = "gpt-4o-mini", api_key, temperature = 0) {
  url <- "https://api.openai.com/v1/chat/completions"

  headers <- c(
    "Content-Type" = "application/json",
    "Authorization" = paste("Bearer", api_key)
  )

  data <- list(
    model = model,
    messages = list(
      list(
        role = "system",
        content = "You are an expert immunologist. Please provide a concise cell type annotation based on the given gene expression information."
      ),
      list(
        role = "user",
        content = description
      )
    ),
    temperature = temperature
  )

  response <- httr::POST(url, httr::add_headers(.headers = headers), body = data, encode = "json")
  content <- httr::content(response, "parsed")

  if (httr::status_code(response) != 200) {
    stop("Error in OpenAI API call: ", content$error$message)
  }

  return(content$choices[[1]]$message$content)
}


#' Annotate Clusters
#'
#' This function annotates cell clusters based on their top up-regulated and down-regulated genes.
#'
#' @param marker.list A list of data frames containing gene markers for each cluster.
#' @param n_genes Number of top genes to select for both up- and down-regulated lists.
#' @param base_prompt A base prompt to use for formatting the annotation data.
#' @param model The model to use for annotation.
#' @param get_reason boolean. Re-prompts LLM for the annotation reason.
#' @param url The URL for the annotation service.
#' @param temperature a number between 0-1, creates diversity for the response
#' @param use_openai boolean. Use OpenAI API instead of local LLM.
#' @param openai_api_key string. OpenAI API key (required if use_openai is TRUE).
#' @return A list of annotations and reasons for each cell cluster.
#'
#' @export
ceLLama <- function(
    marker.list,
    n_genes = 20,
    seed = 101,
    base_prompt = "Act like an expert immunologist and give me the cell type annotation for this cluster. Please, reply only with the answer and nothing else! If you're not sure just label it as 'unsure'.",
    model = "llama3.1",
    get_reason = FALSE,
    url = "http://localhost:11434/api/generate",
    temperature = 0.1,
    use_openai = FALSE,
    openai_api_key = NULL) {
  # Select top genes
  top_genes <- lapply(marker.list, select_top_genes, n = n_genes)

  # Format annotation data
  annotation_data <- format_annotation_data(top_genes, base_prompt)

  # Annotate all clusters and get reasons
  annotations <- lapply(annotation_data, function(description) {
    annotation_prompt <- paste0(description, base_prompt)

    if (use_openai) {
      if (is.null(openai_api_key)) {
        stop("OpenAI API key is required when use_openai is set to TRUE")
      }
      message(model, " is being used!")
      annotation <- get_openai_annotation(annotation_prompt, api_key = openai_api_key, model, temperature = temperature)
    } else {
      annotation <- get_annotation(annotation_prompt, model, url, seed = seed, temperature = temperature)
    }

    message(">> Response: ", annotation)

    reason <- NULL
    if (get_reason) {
      if (use_openai) {
        reason <- get_openai_annotation(paste(
          "The annotation for the cell cluster is:", annotation,
          ". Can you provide the reason for this annotation based on the following description:",
          description
        ), api_key = openai_api_key, temperature = temperature)
      } else {
        reason <- get_reason(annotation, description, model, url, seed = seed, temperature = temperature)
      }
    }
    list(annotation = annotation, description = description, reason = reason)
  })

  return(annotations)
}
