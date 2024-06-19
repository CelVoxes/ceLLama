import pandas as pd
import requests

def select_top_genes(degs, n=20):
    """
    Select top up- and down-regulated genes.
    
    Parameters:
    degs (pd.DataFrame): A dataframe containing gene markers with their log2 fold changes.
    n (int): Number of top genes to select for both up- and down-regulated lists.
    
    Returns:
    dict: A dictionary with top up-regulated and down-regulated genes.
    """
    up_genes = degs.nlargest(n, 'avg_log2FC')['gene'].tolist()
    down_genes = degs.nsmallest(n, 'avg_log2FC')['gene'].tolist()
    return {'up': up_genes, 'down': down_genes}

def format_annotation_data(top_genes, base_prompt):
    """
    Format annotation data.
    
    Parameters:
    top_genes (dict): A dictionary of top genes for each cluster.
    base_prompt (str): A base prompt to use for formatting the annotation data.
    
    Returns:
    dict: A dictionary of formatted annotation prompts for each cluster.
    """
    annotation_data = {}
    for cluster, genes in top_genes.items():
        up_genes = ", ".join(genes['up'])
        down_genes = ", ".join(genes['down'])
        prompt = f"This cell cluster ({cluster}) has up-regulated genes: {up_genes} and down-regulated genes: {down_genes}."
        annotation_data[cluster] = prompt
    return annotation_data

def get_annotation(description, model, url, seed=100, temperature=0):
    """
    Get annotation for a cluster.
    
    Parameters:
    description (str): The description for the cluster to be annotated.
    model (str): The model to use for annotation.
    url (str): The URL for the annotation service.
    
    Returns:
    str: The annotation for the cluster.
    """
    data = {
        "model": model,
        "stream": False,
        "prompt": description,
        "options": {
            "seed": seed,
            "temperature": temperature
        }
    }
    response = requests.post(url, json=data)
    content = response.json()
    return content.get('response')

def get_reason_fn(annotation, description, model, url, seed=100, temperature=0, verbose=True):
    """
    Get reason for the annotation.
    
    Parameters:
    annotation (str): The annotation received for the cluster.
    description (str): The description used to generate the annotation.
    model (str): The model to use for the reason generation.
    url (str): The URL for the reason generation service.
    
    Returns:
    str: The reason for the annotation.
    """
    reason_prompt = f"The annotation for the cell cluster is: {annotation}. Can you provide the reason for this annotation based on the following description: {description}"
    data = {
        "model": model,
        "stream": False,
        "prompt": reason_prompt,
        "options": {
            "seed": seed,
            "temperature": temperature
        }
    }
    response = requests.post(url, json=data)
    content = response.json()
    if verbose:
        print(f">> Reason: {content.get('response')}")
    return content.get('response')

def ceLLama(marker_list, n_genes=20, seed=101, 
            base_prompt="Act like an expert immunologist and give me the cell type annotation for this cluster. Please, reply only with the answer and nothing else! If you're not sure just label it as 'unsure'.", 
            model="llama3", 
            get_reason=False, 
            url="http://localhost:11434/api/generate", 
            temperature=0.1):
    """
    Annotate cell clusters based on their top up-regulated and down-regulated genes.
    
    Parameters:
    marker_list (dict): A dictionary of dataframes containing gene markers for each cluster.
    n_genes (int): Number of top genes to select for both up- and down-regulated lists.
    seed (int): Seed for reproducibility.
    base_prompt (str): A base prompt to use for formatting the annotation data.
    model (str): The model to use for annotation.
    get_reason (bool): Whether to get the reason for the annotation.
    url (str): The URL for the annotation service.
    temperature (float): A number between 0-1 to create diversity in the response.
    
    Returns:
    dict: A dictionary of annotations and reasons for each cell cluster.
    """
    top_genes = {cluster: select_top_genes(degs, n=n_genes) for cluster, degs in marker_list.items()}
    annotation_data = format_annotation_data(top_genes, base_prompt)
    
    annotations = {}
    for cluster, description in annotation_data.items():
        annotation_prompt = f"{description} {base_prompt}"
        annotation = get_annotation(annotation_prompt, model, url, seed=seed, temperature=temperature)
        print(f">> Response: {annotation}")
        
        reason = None
        if get_reason:
            reason = get_reason_fn(annotation, description, model, url, seed=seed, temperature=temperature)
        
        annotations[cluster] = {'annotation': annotation, 'description': description, 'reason': reason}
    
    return annotations
