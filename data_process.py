# data_process.py

import os
import pandas as pd
from collections import defaultdict

def process_all_data(cache_path):
    """
    Process all data from the specified cache directory.
    
    Args:
        cache_path (str): The path to the cache directory.
    """
    if not os.path.exists(cache_path):
        print(f"Cache path does not exist: {cache_path}")
        return

    all_gene_dataframes = {}  # Dictionary to hold the dataframes
    grouped_gene_sets = defaultdict(list)  # Dictionary to group gene sets by file prefix

    # Iterate through all files in the cache directory
    for filename in os.listdir(cache_path):
        file_path = os.path.join(cache_path, filename)

        # Ensure it's a file and not a directory
        if os.path.isfile(file_path):
            # Get the filename without extension and append '_df' to create a variable name
            base_name = os.path.splitext(filename)[0]
            variable_name = f"{base_name}_df"
            print(f"Processing file: {file_path}")
            
            # Read the CSV file into a DataFrame
            dataframe = pd.read_csv(file_path)
            all_gene_dataframes[variable_name] = dataframe
            print(f"Loaded {variable_name} with {dataframe.shape[0]} records.")

            # Extract the file prefix for grouping (take the first two parts of the name separated by underscore)
            prefix_parts = base_name.split("_")
            prefix = "_".join(prefix_parts[:2])  # Take only the first two parts
            if 'Gene_Name' in dataframe.columns:
                gene_set = set(dataframe['Gene_Name'].str.strip().str.upper())  # Convert to uppercase for consistency
                grouped_gene_sets[prefix].append(gene_set)
    
    # Now for each group, take the union of the gene sets
    combined_gene_sets = {}
    for prefix, gene_set_list in grouped_gene_sets.items():
        combined_gene_sets[prefix] = set.union(*gene_set_list)  # Take union within each group
        print(f"Combined gene set for {prefix}: {len(combined_gene_sets[prefix])} genes")
    
    # After combining by prefix, now take the intersection of all combined gene sets
    if combined_gene_sets:
        common_genes = set.intersection(*combined_gene_sets.values())
        print(f"Number of common genes after intersection: {len(common_genes)}")
        print(common_genes)

        # Convert the intersection back to a list, if needed
        common_genes_list = list(common_genes)
        print(f"Common genes list has {len(common_genes_list)} items.")

    common_genes_df = pd.DataFrame(common_genes_list, columns=['Gene_Name'])
    print(common_genes_df)
    common_genes_df = common_genes_df.sort_values('Gene_Name').reset_index(drop=True)
    common_genes_df.to_csv('common_genes.csv', index=False)

    # print(all_gene_dataframes)

    for variable_name, df in all_gene_dataframes.items():
        # Perform an inner merge to keep only common genes
        if 'Gene_Name' in df.columns:
            # Merge with `common_genes_df` to filter only the common genes
            filtered_df = pd.merge(common_genes_df, df, on='Gene_Name', how='inner')
            all_gene_dataframes[variable_name] = filtered_df  # Replace the original dataframe with the filtered one
            print(f"Filtered {variable_name}, new shape: {filtered_df.shape}")

    # print('Filtered all gene dataframes:')
    # print(all_gene_dataframes)
    
    print("Finished processing all files in the cache.")

    ROSMAP_biospecimen = pd.read_csv('E:\LabWork\mosGraphFlow\ROSMAP-raw\Meta-Data\ROSMAP_biospecimen_metadata.csv', sep=',')
    
    # Split the 'specimenID' and construct the matching ID
    ROSMAP_biospecimen['matching_id'] = ROSMAP_biospecimen['specimenID'].apply(lambda x: '.'.join(x.split('.')[2:4]))

    # TODO: Load the clinical data
    survival = pd.read_csv('E:\LabWork\mosGraphFlow\ROSMAP-raw\ROSMAP_clinical\ROSMAP_clinical.csv', sep=',')
    
    original_column_mappings = {}

    # Process Protein DataFrames and store the original column mappings
    matching_id_to_individual_protein = dict(zip(ROSMAP_biospecimen['matching_id'], ROSMAP_biospecimen['individualID']))
    for variable_name, df in all_gene_dataframes.items():
        if "protein" in variable_name.lower():  # Check if the variable name contains "protein"
            print(f"Processing DataFrame: {variable_name}")

            # Save original column mappings (before updating)
            original_column_mappings[variable_name] = {
                'original_columns': df.columns.tolist(),
                'individual_ids': [matching_id_to_individual_protein.get(col, col) for col in df.columns]
            }

            # Update the column names using the dictionary
            df.columns = [matching_id_to_individual_protein.get(col, col) for col in df.columns]

            # Update the DataFrame in the dictionary
            all_gene_dataframes[variable_name] = df
            print(f"Updated column names for {variable_name}")

    # Process Methylation DataFrames and store the original column mappings
    matching_id_to_individual_methy = dict(zip(ROSMAP_biospecimen['specimenID'], ROSMAP_biospecimen['individualID']))
    for variable_name, df in all_gene_dataframes.items():
        if "methy" in variable_name.lower():  # Check if the variable name contains "methy"
            print(f"Processing DataFrame: {variable_name}")

            # Save original column mappings (before updating)
            original_column_mappings[variable_name] = {
                'original_columns': df.columns.tolist(),
                'individual_ids': [matching_id_to_individual_methy.get(col, col) for col in df.columns]
            }

            # Update the column names using the dictionary
            df.columns = [matching_id_to_individual_methy.get(col, col) for col in df.columns]

            # Update the DataFrame in the dictionary
            all_gene_dataframes[variable_name] = df
            print(f"Updated column names for {variable_name}")

    
    # Process Gene Expression DataFrames and store the original column mappings
    def process_gene_column_name(col_name):
        parts = col_name.split('_')
        return '_'.join(parts[:2]) if len(parts) > 1 else col_name

    matching_id_to_individual_gene = dict(zip(ROSMAP_biospecimen['specimenID'], ROSMAP_biospecimen['individualID']))

    for variable_name, df in all_gene_dataframes.items():
        if "expression" in variable_name.lower():  # Check if the variable name contains "expression"
            print(f"Processing DataFrame: {variable_name}")

            # First process the gene column name
            processed_columns = [process_gene_column_name(col) for col in df.columns]
            
            # Save original column mappings (before updating)
            original_column_mappings[variable_name] = {
                'original_columns': processed_columns,  # Use processed column names for mapping
                'individual_ids': [matching_id_to_individual_gene.get(col, col) for col in processed_columns]
            }

            # Then update the column names using the dictionary
            df.columns = [matching_id_to_individual_gene.get(col, col) for col in processed_columns]

            # Update the DataFrame in the dictionary
            all_gene_dataframes[variable_name] = df
            print(f"Updated column names for {variable_name}")

    # Process CNV DataFrames and store the original column mappings
    matching_id_to_individual_cnv = dict(zip(ROSMAP_biospecimen['specimenID'], ROSMAP_biospecimen['individualID']))
    for variable_name, df in all_gene_dataframes.items():
        if "cnv" in variable_name.lower():  # Check if the variable name contains "cnv"
            print(f"Processing DataFrame: {variable_name}")

            # Save original column mappings (before updating)
            original_column_mappings[variable_name] = {
                'original_columns': df.columns.tolist(),
                'individual_ids': [matching_id_to_individual_cnv.get(col, col) for col in df.columns]
            }

            # Update the column names using the dictionary
            df.columns = [matching_id_to_individual_cnv.get(col, col) for col in df.columns]

            # Update the DataFrame in the dictionary
            all_gene_dataframes[variable_name] = df
            print(f"Updated column names for {variable_name}")

    # At this point, you can create the mapping DataFrames using the original_column_mappings dictionary
    map_dataframes = {}

    # Create mapping DataFrames based on original column mappings
    for variable_name, mapping in original_column_mappings.items():
        original_cols = mapping['original_columns']
        individual_ids = mapping['individual_ids']
        
        # Create the mapping DataFrame and add to map_dataframes dictionary
        mapping_df = pd.DataFrame({
            'original_column': original_cols,
            'individualID': individual_ids
        })
        
        # Save to the map_dataframes dictionary
        map_dataframes[f"{variable_name}_map"] = mapping_df

    # View mapping dataframes in map_dataframes dictionary
    for map_name, map_df in map_dataframes.items():
        print(f"{map_name}:\n{map_df.head()}")

    # Collect all 'individualID' columns from map_dataframes and survival
    combined_individualID = pd.concat([
        map_df['individualID'] for map_df in map_dataframes.values()
    ] + [survival['individualID']]).unique()

    # Create a DataFrame to hold the union of individualIDs
    union_map = pd.DataFrame(combined_individualID, columns=['individualID'])

    # Merge with each original map DataFrame in map_dataframes, avoiding column name conflicts
    for map_name, map_df in map_dataframes.items():
        # Merge each map DataFrame with union_map based on 'individualID'
        # Use custom suffixes to avoid conflicts with duplicate column names
        union_map = union_map.merge(map_df, on='individualID', how='outer', suffixes=('', f'_{map_name}'))
        print(f"Merged {map_name} with union_map, current shape: {union_map.shape}")

    # Merge with the survival data (including 'projid' and 'Study' columns)
    union_map = union_map.merge(survival[['individualID', 'projid', 'Study']], on='individualID', how='outer')
    union_map.to_csv('union_map.csv', index=False)
    print(union_map)
    print(all_gene_dataframes)

    union_map_cleaned = union_map.dropna()
    union_map_cleaned.reset_index(drop=True, inplace=True)
    print(union_map_cleaned)
    union_map_cleaned.to_csv('union_map_cleaned.csv', index=False)