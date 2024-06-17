#' Function to generate a simple report in Markdown format
#' 
#' @param annotations results of ceLLama, cell type annotations 
#' 
#' @export
generate_report_md <- function(annotations, output_file = "report.md") {
    report <- lapply(seq_along(annotations), function(i) {
        annotation <- annotations[[i]]$annotation
        description <- annotations[[i]]$description
        reason <- annotations[[i]]$reason
        
        paste("### Cluster", i-1, 
              "\n**Annotation:**", annotation, "\n",
              "\n\n**Description:**", description, "\n",
              "\n\n**Reason:**", reason, "\n")
    })
    
    report_content <- paste(paste(unlist(report), collapse = "\n"), "\n> created by [ceLLama](https://github.com/eonurk/ceLLama)")
    writeLines(report_content, con = output_file)
}

#' Function to convert Markdown to HTML
#' 
#' @export
create_html_report <- function(md_file = "report.md", html_file = "report.html") {
    rmarkdown::render(md_file, output_file = html_file, output_format = "html_document")
}
