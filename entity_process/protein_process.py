# protein_process.py

import pandas as pd
import numpy as np

def process_protein(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Protein data."""
    print(f"Processing Protein: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}")
    
    # Load the protein entity data
    protein_entity_data = pd.read_csv('data/MedGraphica/Node/Protein/medgraphica_protein.csv')    
        
    # Determine the separator based on the file extension
    if file_path.endswith('.txt') or file_path.endswith('.tsv'):
        sep = '\t'
    elif file_path.endswith('.csv'):
        sep = ','
    else:
        raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")
    
    # Load the protein data file
    protein = pd.read_csv(file_path, sep=sep)

    if id_type == 'Gene_Name':
        protein_entity_data = protein_entity_data[['MedGraphica_ID', 'Gene_Name']].copy()
    else:
        protein_entity_data = protein_entity_data[[id_type, 'MedGraphica_ID', 'Gene_Name']].copy()

    # Merge with protein data on selected_column
    protein_merged = pd.merge(
        protein,
        protein_entity_data,
        left_on=selected_column,
        right_on=id_type,
        how='inner'
    )

    # Drop rows where 'MedGraphica_ID' is NaN (if any)
    protein_merged.dropna(subset=['MedGraphica_ID'], inplace=True)

    # Ensure 'MedGraphica_ID' and 'Gene_Name' are the first and second columns, respectively
    cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in protein_merged.columns if col not in ['MedGraphica_ID', 'Gene_Name', selected_column, id_type]]
    protein_data = protein_merged[cols]

    # Select numeric columns
    numeric_cols = protein_data.select_dtypes(include=[np.number]).columns.tolist()

    # Group by 'MedGraphica_ID' and calculate the mean for all numeric columns
    protein_data = protein_data.groupby('MedGraphica_ID')[numeric_cols].mean()

    # Reset the index to turn 'MedGraphica_ID' back into a column
    protein_data.reset_index(inplace=True)

    # Merge back the gene names after grouping by 'MedGraphica_ID'
    protein_data = pd.merge(protein_data, protein_entity_data[['MedGraphica_ID', 'Gene_Name']], on='MedGraphica_ID', how='left')

    # Ensure 'MedGraphica_ID' is the first column, "Gene_Name" is the second column
    final_cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in protein_data.columns if col not in ['MedGraphica_ID', 'Gene_Name']]
    protein_data = protein_data[final_cols]

    # Export the result to a CSV file
    output_file_path = f'cache/{feature_label.lower()}.csv'
    protein_data.to_csv(output_file_path, index=False)

    print(f"Processed protein data saved to {output_file_path}")