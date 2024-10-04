# protein_process.py

import pandas as pd
import numpy as np

def process_protein(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Protein data."""
    print(f"Processing Protein: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}")
    
    # Load the protein entity data
    protein_entity_data = pd.read_csv('medgraphica_data/medgraphica_protein.csv')    
        
    # Determine the separator based on the file extension
    if file_path.endswith('.txt') or file_path.endswith('.tsv'):
        sep = '\t'
    elif file_path.endswith('.csv'):
        sep = ','
    else:
        raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")
    
    # Load the protein data file
    protein = pd.read_csv(file_path, sep=sep)

    # Rename 'Unnamed: 0' to 'ID_Column' if it exists
    if 'Unnamed: 0' in protein.columns:
        protein.rename(columns={'Unnamed: 0': 'ID_Column'}, inplace=True)
    
    # Update selected_column if it was 'Unnamed: 0'
    if selected_column == 'Unnamed: 0':
        selected_column = 'ID_Column'

    # Split the selected_column into protein_names and protein_ids
    protein['Gene_Name'] = protein[selected_column].apply(lambda x: x.split('|')[0])  # Extract gene names
    protein['Protein_ID'] = protein[selected_column].apply(lambda x: x.split('|')[1].split('-')[0])  # Extract protein ids and handle cases with '-'
    
    # Select numeric columns for aggregation
    numeric_cols = protein.select_dtypes(include=[float, int]).columns.tolist()

    # Group by 'Gene_Name' and calculate the mean of each numeric column
    protein_grouped = protein.groupby('Gene_Name')[numeric_cols].mean().reset_index()
    
    # Merge back 'protein_id' to keep them as markers
    # Use drop_duplicates to ensure we don't multiply the data unnecessarily
    protein_grouped = pd.merge(protein_grouped, protein[['Gene_Name', 'Protein_ID']].drop_duplicates(), on='Gene_Name', how='left')

    # Merge with protein_entity_data to add 'MedGraphica_ID'
    protein_grouped = pd.merge(protein_grouped, 
                               protein_entity_data[[id_type, 'MedGraphica_ID']], 
                               left_on='Protein_ID', 
                               right_on=id_type, 
                               how='left')
    # Drop rows where 'MedGraphica_ID' is NaN
    protein_grouped.dropna(subset=['MedGraphica_ID'], inplace=True)
    # Drop the extra merge key (id_type from protein_entity_data)
    protein_grouped.drop(columns=[id_type], inplace=True)

    # Ensure 'MedGraphica_ID' is the first column, 'protein_id' second, 'Gene_Name' third
    final_cols = ['MedGraphica_ID', 'Protein_ID', 'Gene_Name'] + [col for col in protein_grouped.columns if col not in ['MedGraphica_ID', 'Protein_ID', 'Gene_Name']]
    protein_grouped = protein_grouped[final_cols]
    print("Grouped protein data:", protein_grouped.shape)

    # Export the result to a CSV file
    output_file_path = f'cache/{entity_type.lower()}_{feature_label.lower()}.csv'
    protein_grouped.to_csv(output_file_path, sep=",", index=False)

    print(f"Protein data processing completed. Output saved to {output_file_path}")