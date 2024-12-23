# data_process.py

import os
import pandas as pd
import numpy as np
import glob
from collections import defaultdict

def merge_data_and_generate_entity_mapping(cache_folder, file_order):
    # Initialize dictionary and list for storing dataframes and sample_id sets
    all_sample_dataframes = {}
    sample_id_sets = []

    # Define processed data path based on the given cache_folder
    processed_data_path = os.path.join(cache_folder, 'processed_data/')
    os.makedirs(processed_data_path, exist_ok=True)

    print('File order:', file_order)

    # Find the clinical data path within the cache folder
    def get_clinical_data_path(cache_folder):
        """
        Locate the clinical data file in the cache folder by identifying the
        unique folder that starts with '_y' and contains a single .csv file.
        """
        # Locate '_y' prefixed folders within cache_folder
        y_folders = [f for f in os.listdir(cache_folder) if f.startswith('_y') and os.path.isdir(os.path.join(cache_folder, f))]
        
        # Ensure there is exactly one '_y' folder
        if len(y_folders) != 1:
            raise ValueError("There should be exactly one '_y' folder in the cache folder.")

        # Path to the identified '_y' folder
        y_folder_path = os.path.join(cache_folder, y_folders[0])

        # Find all .csv files within the '_y' folder
        csv_files = glob.glob(os.path.join(y_folder_path, '*.csv'))

        # Ensure there is exactly one .csv file
        if len(csv_files) != 1:
            raise ValueError("The '_y' folder should contain exactly one .csv file.")

        # Return the path to the unique .csv file
        return csv_files[0]

    # Integrate the clinical data path retrieval and read the data
    try:
        clinical_data_path = get_clinical_data_path(cache_folder)
        print('Reading clinical data...')
        clinical_data = pd.read_csv(clinical_data_path)
        print(f"Clinical data successfully loaded from: {clinical_data_path}")
    except ValueError as e:
        print(f"Error: {e}")

    # Traverse all files in the cache folder
    for file_name in os.listdir(cache_folder):
        if file_name.endswith('.csv'):
            file_path = os.path.join(cache_folder, file_name)
            df = pd.read_csv(file_path).sort_values(by='Sample_ID')
            variable_name = os.path.splitext(file_name)[0] + '_df'
            all_sample_dataframes[variable_name] = df
            sample_ids = set(df.iloc[:, 0])
            sample_id_sets.append(sample_ids)

    # Calculate the intersection of all sample_id sets
    if sample_id_sets:
        common_sample_ids = set.intersection(*sample_id_sets)
        print(f"Found {len(common_sample_ids)} common Sample_IDs across all files.")
        
        # Export the common_sample_ids to a CSV file
        common_id_df = pd.DataFrame({'Sample_ID': sorted(common_sample_ids)})
        common_id_df.to_csv(os.path.join(processed_data_path, 'common_sample_ids.csv'), index=False)
    else:
        print("No files found or no Sample_IDs extracted.")
        return None, None, processed_data_path
    
    # Filter and save clinical data
    clinical_data = clinical_data[clinical_data['Sample_ID'].isin(common_sample_ids)].sort_values(by='Sample_ID')
    clinical_data.to_csv(os.path.join(processed_data_path, 'filtered_clinical_data.csv'), index=False)
    np.save(os.path.join(processed_data_path, 'yAll.npy'), clinical_data.iloc[:, 1:].astype(np.float64).to_numpy())

    # Filter the dataframes in the dictionary based on common_sample_ids and merge them
    merged_data, entity_index_id_mapping = None, []
    for file_name in file_order:
        variable_name = file_name + '_df'
        if variable_name in all_sample_dataframes:
            df = all_sample_dataframes[variable_name]
            df = df[df['Sample_ID'].isin(common_sample_ids)]
            entity_index_id_mapping.extend(df.columns[1:])
            data_array = df.iloc[:, 1:].values.astype(np.float64)
            merged_data = data_array if merged_data is None else np.concatenate((merged_data, data_array), axis=1)

    # Save the merged data to a NumPy file
    np.save(os.path.join(processed_data_path, 'xAll.npy'), merged_data)
    
    # Create and save the entity index mapping
    entity_index_id_mapping_df = pd.DataFrame({
        'Index': range(len(entity_index_id_mapping)),
        'BioMedGraphica_ID': entity_index_id_mapping
    })
    entity_index_id_mapping_df.to_csv(os.path.join(processed_data_path, 'entity_index_id_mapping.csv'), index=False)

    return entity_index_id_mapping_df, merged_data, processed_data_path


def filter_and_save_edge_data(database_path, entity_index_id_mapping):
    """Filter edge data based on entity index and return unique types."""
    edge_csv_path = os.path.join(database_path, 'Relation', 'biomedgraphica_relation.csv')
    edge_data_raw = pd.read_csv(edge_csv_path)

    edge_data = edge_data_raw[['From_ID', 'To_ID', 'Type']].copy()

    # Filter edge data based on entity_index_id_mapping
    filtered_edge_data = edge_data[
        edge_data['From_ID'].isin(entity_index_id_mapping['BioMedGraphica_ID']) &
        edge_data['To_ID'].isin(entity_index_id_mapping['BioMedGraphica_ID'])
    ]

    # Get unique types from the filtered edge data
    unique_types = filtered_edge_data['Type'].unique().tolist()
    
    return filtered_edge_data, unique_types

def process_edge_data_with_selected_types(filtered_edge_data, selected_types, entity_index_id_mapping, processed_data_path):
    """Process edge data based on selected types and save filtered data."""
    
    # Filter edge data based on selected types
    filtered_edge_data = filtered_edge_data[filtered_edge_data['Type'].isin(selected_types)]
    
    # Create From_Index and To_Index columns
    entity_index = pd.Series(entity_index_id_mapping['Index'].values, index=entity_index_id_mapping['BioMedGraphica_ID'])
    filtered_edge_data.loc[:, 'From_Index'] = filtered_edge_data['From_ID'].map(entity_index)
    filtered_edge_data.loc[:, 'To_Index'] = filtered_edge_data['To_ID'].map(entity_index)
    
    # Sort the edge data by From_Index
    filtered_edge_data = filtered_edge_data.sort_values(by=['From_Index'])

    # Save edge_index as a NumPy array
    edge_index = filtered_edge_data[['From_Index', 'To_Index']].T
    np.save(os.path.join(processed_data_path, 'edge_index.npy'), edge_index.to_numpy().astype(np.int64))

    # Save the filtered edge data with BioMedGraphica_ID
    filtered_edge_data.to_csv(os.path.join(processed_data_path, 'filtered_edge_id_index_data.csv'), index=False)


