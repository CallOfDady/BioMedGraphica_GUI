# phenotype_process.py

import pandas as pd
import numpy as np

def process_phenotype(entity_type, id_type, file_path, selected_column, feature_label, matcher, selection_callback=None):
    """Process Phenotype data, allow user selection, and perform further processing."""
    print(f"Processing Phenotype: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}")

    # Determine the separator based on the file extension
    if file_path.endswith('.txt') or file_path.endswith('.tsv'):
        sep = '\t'
    elif file_path.endswith('.csv'):
        sep = ','
    else:
        raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")

    phenotype = pd.read_csv(file_path, sep=sep)

    all_topk_phenotypes = {}

    # Get topk phenotypes for each phenotype value in the selected column
    for phenotype_value in phenotype[selected_column].dropna():
        topk_phenotypes = matcher.get_topk_entities(query=phenotype_value, k=3, embeddings=phenotype_embeddings)
        all_topk_phenotypes[phenotype_value] = topk_phenotypes
    print(f"Phenotype: {all_topk_phenotypes}")

    # Allow user selection of the topk phenotypes
    if selection_callback:
        updated_topk_phenotypes = selection_callback(all_topk_phenotypes, entity_type, id_type, file_path, selected_column, feature_label, phenotype)
        if updated_topk_phenotypes is not None:
            all_topk_phenotypes = updated_topk_phenotypes
        else:
            print("Warning: selection_callback returned None")