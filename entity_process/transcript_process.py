# transcript_process.py

import pandas as pd
import numpy as np

def process_transcript(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Transcript data."""
    print(f"Processing Transcript: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}")
    
    transcript_entity_data = pd.read_csv('medgraphica_data/medgraphica_transcript.csv')

    # Determine the separator based on the file extension
    if file_path.endswith('.txt') or file_path.endswith('.tsv'):
        sep = '\t'
    elif file_path.endswith('.csv'):
        sep = ','
    else:
        raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.") 

    def process_gene_expression(transcript_entity_data):
        """Process Gene Expression data."""
        print(f"Processing Gene Expression for {entity_type} with Feature Label: {feature_label}")

    # Check if the feature_label contains "expression" or "cnv"
    if "expression" in feature_label.lower():
        process_gene_expression(transcript_entity_data)
    elif "xxx" in feature_label.lower():
        print(f"Processing xxx data for Transcript: {entity_type} with Feature Label: {feature_label}")
    else:
        print(f"Processing for Feature Label: {feature_label} is not implemented yet.")
        # Leave space for other feature_label handling in the future