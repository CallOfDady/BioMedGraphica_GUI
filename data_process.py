# data_process.py

import os
import pandas as pd
import numpy as np
from collections import defaultdict

def process_and_merge_data(cache_folder, clinical_data_path):
    """
    Main function to process files, sort by Patient_ID, find common Patient_IDs, filter and merge data.
    
    Args:
        cache_folder (str): The path to the cache directory.
    """
    # Step 1: Initialize dictionary and list for storing dataframes and sample_id sets
    all_sample_dataframes = {}
    sample_id_sets = []
    
    print('Reading clinical data...')
    clinical_data = pd.read_csv(clinical_data_path)
    
    edge_data = pd.read_csv('data/BioMedGraphica/Interaction/biomedgraphica_edge.csv')
    filtered_edge_data = edge_data[['From_ID', 'To_ID', 'Type']].copy()
    
    # Step 2: Traverse all files in the cache folder
    for file_name in os.listdir(cache_folder):
        if file_name.endswith('.csv'):
            file_path = os.path.join(cache_folder, file_name)
            
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path)
            
            # Create a variable name based on the file name (removing the file extension and adding '_df')
            variable_name = os.path.splitext(file_name)[0] + '_df'
            print(f"Processing file: {file_name} as variable: {variable_name}")
            
            # Sort the DataFrame by 'Patient_ID'
            df = df.sort_values(by='Patient_ID')
            
            # Add the DataFrame to the dictionary with the variable name as the key
            all_sample_dataframes[variable_name] = df
            
            # Extract the first column (assumed to be 'Patient_ID') and store it in the sample_id set
            sample_ids = set(df.iloc[:, 0])  # Assuming first column contains the 'Patient_ID'
            sample_id_sets.append(sample_ids)
    
    # Step 3: Calculate the intersection of all sample_id sets
    if sample_id_sets:
        common_sample_ids = set.intersection(*sample_id_sets)
        print(f"Found {len(common_sample_ids)} common Patient_IDs across all files.")
        
        # Export the common_sample_ids to a CSV file
        common_id_df = pd.DataFrame({'Patient_ID': list(common_sample_ids)})
        common_id_df = common_id_df.sort_values(by='Patient_ID')
        common_id_df.to_csv('common_sample_ids.csv', index=False)
        print("Common Patient_IDs have been saved to 'common_sample_ids.csv'.")
    else:
        print("No files found or no Patient_IDs extracted.")
        return
    
    clinical_data = clinical_data[clinical_data['Patient_ID'].isin(common_sample_ids)]
    clinical_data = clinical_data.sort_values(by='Patient_ID')
    clinical_data.to_csv('filtered_clinical_data.csv', index=False)
    clinical_data_array = clinical_data.to_numpy()
    clinical_data_array = clinical_data_array[:, 1:].astype(np.float64)
    np.save('yAll.npy', clinical_data_array)

    # Step 4: Filter the dataframes in the dictionary based on common_sample_ids
    for key in all_sample_dataframes.keys():
        df = all_sample_dataframes[key]
        all_sample_dataframes[key] = df[df['Patient_ID'].isin(common_sample_ids)]
        print(f"Filtered {key} based on common Patient_IDs.")
    
    # Step 5: Merge files in the specified order and save as a NumPy file
    file_order = ['protein', 'cnv', 'gene_expression', 'methylation']
    merged_data = None

    gene_id_mapping = [] 
    
    for file_name in file_order:
        variable_name = file_name + '_df'
        if variable_name in all_sample_dataframes:
            df = all_sample_dataframes[variable_name]

            gene_id_mapping.extend(df.columns[1:])

            data_array = df.values[:, 1:].astype(np.float64)
            
            # Convert the DataFrame to NumPy array, excluding the 'Patient_ID' column
            if merged_data is None:
                merged_data = data_array # Start with the first file's data (excluding 'Patient_ID' column)
            else:
                # Concatenate the remaining data (excluding 'Patient_ID' column) along axis=1 (columns)
               merged_data = np.concatenate((merged_data, data_array), axis=1)

    np.save('xAll.npy', merged_data)
    print("Merged data has been saved.")

    gene_id_mapping_df = pd.DataFrame({
        'Index': range(len(gene_id_mapping)),
        'BioMedGraphica_ID': gene_id_mapping
    })
    gene_id_mapping_df.to_csv('gene_id_mapping.csv', index=False)

    gene_id_mapping_index = pd.Series(gene_id_mapping_df['Index'].values, index=gene_id_mapping_df['BioMedGraphica_ID'])

    # Filter the edge data
    filtered_edge_data = filtered_edge_data[
    filtered_edge_data['From_ID'].isin(gene_id_mapping) & filtered_edge_data['To_ID'].isin(gene_id_mapping)
    ]
    filtered_edge_data.to_csv('filtered_edge_id_data.csv', index=False)

    # Map the gene IDs to their respective indices
    filtered_edge_data['From_ID'] = filtered_edge_data['From_ID'].map(gene_id_mapping_index)
    filtered_edge_data['To_ID'] = filtered_edge_data['To_ID'].map(gene_id_mapping_index)
    filtered_edge_data.to_csv('filtered_edge_index_data.csv', index=False)

    filtered_edge_data = filtered_edge_data.drop(columns=['Type'])

    edge_index = filtered_edge_data.T

    edge_index_array = edge_index.to_numpy().astype(np.int64)
    np.save('edge_index.npy', edge_index_array)



# process_and_merge_data('./cache', './test_data/clinical_data.csv')  # Test the data processing pipeline
