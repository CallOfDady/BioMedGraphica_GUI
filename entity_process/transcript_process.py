# transcript_process.py

import pandas as pd
import numpy as np

def process_transcript(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Transcript data."""
    print(f"Processing Transcript: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}")
    
    transcript_entity_data = pd.read_csv('data/MedGraphica/Node/Transcript/medgraphica_transcript.csv')

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

        # Read the gene expression file
        gene_expression = pd.read_csv(file_path, sep=sep)

        # Retain only the necessary columns
        if id_type == 'Gene_Name':
            transcript_entity_data = transcript_entity_data[[id_type, 'MedGraphica_ID']].copy()
        else:
            transcript_entity_data = transcript_entity_data[[id_type, 'MedGraphica_ID', 'Gene_Name']].copy()

        # Merge the gene expression data with the gene entity data
        gene_expression_merged = pd.merge(
            gene_expression,
            transcript_entity_data,
            left_on=selected_column,
            right_on=id_type,
            how='inner'
        )

        # Drop rows where 'MedGraphica_ID' is NaN (if any)
        gene_expression_merged.dropna(subset=['MedGraphica_ID'], inplace=True)

        # Ensure 'MedGraphica_ID' and 'Gene_Name' are the first and second columns, respectively
        cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in gene_expression_merged.columns if col not in ['MedGraphica_ID', 'Gene_Name', selected_column, id_type]]
        gene_expression_data = gene_expression_merged[cols]
        
        # Select numeric columns
        numeric_cols = gene_expression_data.select_dtypes(include=[np.number]).columns.tolist()

        # Group by 'MedGraphica_ID' and calculate the mean for all numeric columns
        gene_expression_data = gene_expression_data.groupby('MedGraphica_ID')[numeric_cols].mean()

        # Reset the index to turn 'MedGraphica_ID' back into a column
        gene_expression_data.reset_index(inplace=True)

        # Merge back the gene names after grouping by 'MedGraphica_ID'
        gene_expression_data = pd.merge(gene_expression_data, transcript_entity_data[['MedGraphica_ID', 'Gene_Name']], on='MedGraphica_ID', how='left')

        # remove everything after the first hyphen
        gene_expression_data['Gene_Name'] = gene_expression_data['Gene_Name'].apply(lambda x: x.split('-')[0] if pd.notnull(x) else x)

        # Ensure 'MedGraphica_ID' is first, 'Gene_Name' is second, and other columns follow the original order
        final_cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in gene_expression_data.columns if col not in ['MedGraphica_ID', 'Gene_Name']]
        gene_expression_data = gene_expression_data[final_cols]

        # Export the result to a CSV file
        output_file_path = f'cache/{feature_label.lower()}.csv'
        gene_expression_data.to_csv(output_file_path, sep=",", index=False)

        print(f"Gene expression data processing completed. Output saved to {output_file_path}")


    # Check if the feature_label contains "expression" or "cnv"
    if "expression" in feature_label.lower():
        process_gene_expression(transcript_entity_data)
    elif "xxx" in feature_label.lower():
        print(f"Processing xxx data for Transcript: {entity_type} with Feature Label: {feature_label}")
    else:
        print(f"Processing for Feature Label: {feature_label} is not implemented yet.")
        # Leave space for other feature_label handling in the future