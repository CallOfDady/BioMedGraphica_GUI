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

    # Update the dictionary mapping with the new matching IDs
    matching_id_to_individual_protein = dict(zip(ROSMAP_biospecimen['matching_id'], ROSMAP_biospecimen['individualID']))

    # Iterate through all dataframes and update the column names
    for variable_name, df in all_gene_dataframes.items():
        if "protein" in variable_name.lower():  # Check if the variable name contains "protein"
            print(f"Processing DataFrame: {variable_name}")

            # Update the column names using the dictionary
            df.columns = [matching_id_to_individual_protein.get(col, col) for col in df.columns]

            # Update the DataFrame in the dictionary
            all_gene_dataframes[variable_name] = df
            print(f"Updated column names for {variable_name}")

    matching_id_to_individual_methy = dict(zip(ROSMAP_biospecimen['specimenID'], ROSMAP_biospecimen['individualID']))
    for variable_name, df in all_gene_dataframes.items():
        if "methy" in variable_name.lower():  # Check if the variable name contains "methy"
            print(f"Processing DataFrame: {variable_name}")

            # Update the column names using the dictionary
            df.columns = [matching_id_to_individual_methy.get(col, col) for col in df.columns]

            # Update the DataFrame in the dictionary
            all_gene_dataframes[variable_name] = df
            print(f"Updated column names for {variable_name}")

    def process_gene_column_name(col_name):
        parts = col_name.split('_')
        return '_'.join(parts[:2]) if len(parts) > 1 else col_name
    
    matching_id_to_individual_gene = dict(zip(ROSMAP_biospecimen['specimenID'], ROSMAP_biospecimen['individualID']))
    for variable_name, df in all_gene_dataframes.items():
        if "expression" in variable_name.lower():  # Check if the variable name contains "expression"
            print(f"Processing DataFrame: {variable_name}")

            df.columns = [process_gene_column_name(col) for col in df.columns]

            # Update the column names using the dictionary
            df.columns = [matching_id_to_individual_gene.get(col, col) for col in df.columns]

            # Update the DataFrame in the dictionary
            all_gene_dataframes[variable_name] = df
            print(f"Updated column names for {variable_name}")
    
    matching_id_to_individual_cnv = dict(zip(ROSMAP_biospecimen['specimenID'], ROSMAP_biospecimen['individualID']))
    for variable_name, df in all_gene_dataframes.items():
        if "cnv" in variable_name.lower():  # Check if the variable name contains "cnv"
            print(f"Processing DataFrame: {variable_name}")

            # Update the column names using the dictionary
            df.columns = [matching_id_to_individual_cnv.get(col, col) for col in df.columns]

            # Update the DataFrame in the dictionary
            all_gene_dataframes[variable_name] = df
            print(f"Updated column names for {variable_name}")

    # print(all_gene_dataframes)

    # TODO: Load the clinical data
    survival = pd.read_csv('E:\LabWork\mosGraphFlow\ROSMAP-raw\ROSMAP_clinical\ROSMAP_clinical.csv', sep=',')

    R_columns_all = []

    for variable_name, df in all_gene_dataframes.items():
        # TODO: Adjust the keywords as needed
        if any(keyword in variable_name.lower() for keyword in ['methy', 'cnv', 'expression', 'protein']): 
            print(f"Extracting 'R' columns from DataFrame: {variable_name}")
        
            # Convert all column names to strings
            R_columns = [str(col) for col in df.columns if str(col).startswith('R')]
     
            R_columns_all.append(set(R_columns))

            print(f"Finished processing DataFrame: {variable_name}")

    # Extract 'R' columns from survival['individualID']
    R_columns_survival = [col for col in survival['individualID'] if col.startswith('R')]

    # Add survival's 'R' columns as a set
    R_columns_all.append(set(R_columns_survival))

    # Find the common 'R' columns across all datasets
    common_R_columns = set.intersection(*R_columns_all)
    print(f"Number of common R columns: {len(common_R_columns)}")
    common_R_columns_list = list(common_R_columns)
    # print(common_R_columns_list)

    patient_id_df = pd.DataFrame(common_R_columns_list, columns=['Patient_ID'])
    patient_id_list = patient_id_df['Patient_ID'].tolist()

    sorted_gene = common_genes_df
    sorted_gene_list = sorted_gene['Gene_Name'].tolist()
    print(sorted_gene)



    # Promoter
    sorted_gene_methy = [gene + '-Promoter' for gene in sorted_gene_list]
    sorted_gene_methy_df = pd.DataFrame(sorted_gene_methy, columns=['Gene'])
    print(sorted_gene_methy_df)

    # Gene
    sorted_gene_gene = [gene + '-Gene' for gene in sorted_gene_list]
    sorted_gene_gene_df = pd.DataFrame(sorted_gene_gene, columns=['Gene'])
    print(sorted_gene_gene_df)

    # Tran
    sorted_gene_tran = [gene + '-Transcript' for gene in sorted_gene_list]
    sorted_gene_tran_df = pd.DataFrame(sorted_gene_tran, columns=['Gene'])
    print(sorted_gene_tran_df)

    # Protein
    sorted_gene_protein = [gene + '-Protein' for gene in sorted_gene_list]
    sorted_gene_protein_df = pd.DataFrame(sorted_gene_protein, columns=['Gene'])
    print(sorted_gene_protein_df)

    # all-gene
    sorted_gene_all = sorted_gene_methy + sorted_gene_gene + sorted_gene_tran  + sorted_gene_protein
    sorted_all_gene_df = pd.DataFrame(sorted_gene_all, columns=['Gene'])
    print(sorted_all_gene_df)

    sorted_patient_id_list = patient_id_df.sort_values(by='Patient_ID')['Patient_ID'].tolist()
    # print(sorted_patient_id_list)
    sorted_patient_id_df = patient_id_df.sort_values(by='Patient_ID').reset_index(drop=True)
    # print(sorted_patient_id_df)

    for variable_name, df in all_gene_dataframes.items():
        if 'MedGraphica_ID' in df.columns and 'Gene_Name' in df.columns:
            available_patient_samples = [col for col in sorted_patient_id_list if col in df.columns]
            
            new_column_order = ['MedGraphica_ID', 'Gene_Name'] + available_patient_samples
            
            df = df[new_column_order]
            
            all_gene_dataframes[variable_name] = df
            
            print(f"Updated column order for {variable_name}:")
            print(df.columns.tolist())

    print(all_gene_dataframes)

    # TODO: user input for the output folder
    # Create the output folder if it doesn't exist
    output_folder = './output'
    processed_data_folder = output_folder + '/ROSMAP_processed_data'
    graph_output_folder = output_folder + '/ROSMAP_graph_data'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(processed_data_folder):
        os.makedirs(processed_data_folder)
    if not os.path.exists(graph_output_folder):
        os.makedirs(graph_output_folder)

    # map_all_gene_df = sorted_all_gene_df.copy()

    # map_all_gene_df['Gene_Num'] = range(1, len(map_all_gene_df) + 1)
    # map_all_gene_df['Node_Type'] = map_all_gene_df['Gene'].apply(lambda x: x.split('-')[0] + '-' + x.split('-')[-1])
    # map_all_gene_df = map_all_gene_df[['Gene_Num', 'Gene', 'Node_Type']]
    # map_all_gene_df.to_csv(os.path.join(graph_output_folder, 'map-all-gene.csv'), index=False)

    sorted_all_gene_dict = sorted_all_gene_df['Gene'].to_dict()
    sorted_all_gene_name_dict = {value: key for key, value in sorted_all_gene_dict.items()}
    num_gene = sorted_gene_tran_df.shape[0]
    num_gene_protein = sorted_gene_protein_df.shape[0]
    nodetype_list =  ['Promoter'] * num_gene + ['Gene'] * num_gene + ['Transcript'] * num_gene + ['Protein'] * num_gene_protein
    map_all_gene_df = pd.DataFrame({'Gene_Num': sorted_all_gene_dict.keys(), 'Gene': sorted_all_gene_dict.values(), 'Node_Type': nodetype_list})
    print(map_all_gene_df)
    map_all_gene_df.to_csv(os.path.join(graph_output_folder, 'map_all_gene.csv'), index=False)


    sorted_methy_list = sorted_gene_methy_df['Gene'].tolist()
    sorted_gene_list = sorted_gene_gene_df['Gene'].tolist()
    sorted_tran_list = sorted_gene_tran_df['Gene'].tolist()
    sorted_protein_list = sorted_gene_protein_df['Gene'].tolist()


    # [Promoter - Gene]
    meth_gene_edge_df = pd.DataFrame({'From': sorted_methy_list, 'To': sorted_gene_list})
    print(meth_gene_edge_df)


    # [Gene - Transcript]
    gene_tran_edge_df = pd.DataFrame({'From': sorted_gene_list, 'To': sorted_tran_list})
    print(gene_tran_edge_df)

    # [Transcript - ProteinEIN]
    tran_protein_edge_df = pd.DataFrame({'From': sorted_tran_list, 'To': sorted_protein_list})
    print(tran_protein_edge_df)

    # replace gene name with gene number
    meth_gene_num_edge_df = meth_gene_edge_df.copy()
    meth_gene_num_edge_df['From'] = meth_gene_edge_df['From'].map(sorted_all_gene_name_dict)
    meth_gene_num_edge_df['To'] = meth_gene_edge_df['To'].map(sorted_all_gene_name_dict)
    print(meth_gene_num_edge_df)

    gene_tran_num_edge_df = gene_tran_edge_df.copy()
    gene_tran_num_edge_df['From'] = gene_tran_edge_df['From'].map(sorted_all_gene_name_dict)
    gene_tran_num_edge_df['To'] = gene_tran_edge_df['To'].map(sorted_all_gene_name_dict)
    print(gene_tran_num_edge_df)

    tran_protein_num_edge_df = tran_protein_edge_df.copy()
    tran_protein_num_edge_df['From'] = tran_protein_edge_df['From'].map(sorted_all_gene_name_dict)
    tran_protein_num_edge_df['To'] = tran_protein_edge_df['To'].map(sorted_all_gene_name_dict)
    print(tran_protein_num_edge_df)

    protein_df = all_gene_dataframes['protein_protein_df']

    # Read the medgraphica_protein_protein.csv file
    medgraphica_protein_protein_df = pd.read_csv('medgraphica_data\Interaction\medgraphica_protein_protein.csv')

    # Extract all MedGraphica_IDs from protein_protein_df
    protein_medgraphica_ids = protein_df['MedGraphica_ID'].unique()

    # Filter interactions in medgraphica_protein_protein_df where Protein_1 matches MedGraphica_IDs
    protein_protein_filtered_interactions = medgraphica_protein_protein_df[medgraphica_protein_protein_df['Protein_1'].isin(protein_medgraphica_ids)]

    # Merge interaction data
    # Map Protein_1 and Protein_2 to their corresponding Gene_Name, then add '-Protein' suffix
    protein_protein_interaction_with_gene_names = protein_protein_filtered_interactions.copy()
    protein_protein_interaction_with_gene_names = protein_protein_interaction_with_gene_names.merge(protein_df[['MedGraphica_ID', 'Gene_Name']],
                                                                    left_on='Protein_1', right_on='MedGraphica_ID',
                                                                    how='left')
    protein_protein_interaction_with_gene_names.rename(columns={'Gene_Name': 'Gene_1'}, inplace=True)

    protein_protein_interaction_with_gene_names = protein_protein_interaction_with_gene_names.merge(protein_df[['MedGraphica_ID', 'Gene_Name']],
                                                                    left_on='Protein_2', right_on='MedGraphica_ID',
                                                                    how='left')
    protein_protein_interaction_with_gene_names.rename(columns={'Gene_Name': 'Gene_2'}, inplace=True)

    # Add '-Protein' suffix to Gene_1 and Gene_2
    protein_protein_interaction_with_gene_names['Gene_1'] = protein_protein_interaction_with_gene_names['Gene_1'] + '-Protein'
    protein_protein_interaction_with_gene_names['Gene_2'] = protein_protein_interaction_with_gene_names['Gene_2'] + '-Protein'

    # Map Gene_1 and Gene_2 to corresponding numeric IDs using sorted_all_gene_name_dict
    protein_protein_interaction_with_gene_names['Gene_1_Num'] = protein_protein_interaction_with_gene_names['Gene_1'].map(sorted_all_gene_name_dict)
    protein_protein_interaction_with_gene_names['Gene_2_Num'] = protein_protein_interaction_with_gene_names['Gene_2'].map(sorted_all_gene_name_dict)

    # Handle NaN values before converting to integers
    protein_protein_interaction_with_gene_names['Gene_1_Num'] = protein_protein_interaction_with_gene_names['Gene_1_Num'].fillna(0).astype(int)
    protein_protein_interaction_with_gene_names['Gene_2_Num'] = protein_protein_interaction_with_gene_names['Gene_2_Num'].fillna(0).astype(int)

    # Select relevant columns and add Gene_1_Name, Gene_2_Name, Gene_1_ID, Gene_2_ID
    protein_protein_interaction_with_gene_names['Gene_1_ID'] = protein_protein_interaction_with_gene_names['Protein_1']
    protein_protein_interaction_with_gene_names['Gene_2_ID'] = protein_protein_interaction_with_gene_names['Protein_2']

    protein_protein_interaction_with_gene_names['Edge_Type'] = 'Protein-Protein'

    protein_protein_edge_num_df = protein_protein_interaction_with_gene_names[
        ['Gene_1_Num', 'Gene_2_Num', 'Gene_1', 'Gene_2', 'Gene_1_ID', 'Gene_2_ID', 'Edge_Type']
    ].dropna()

    protein_protein_edge_num_df.rename(columns={'Gene_1_Num': 'From'}, inplace=True)
    protein_protein_edge_num_df.rename(columns={'Gene_2_Num': 'To'}, inplace=True)

    protein_protein_edge_num_df = protein_protein_edge_num_df[['From', 'To', 'Edge_Type']]

    # Save the data to protein_edge_num.csv file
    protein_protein_edge_num_df.to_csv(os.path.join(graph_output_folder, 'protein_protein_edge_num.csv'), index=False)

    print("Processing complete, the relationships have been saved to the 'protein_edge_num.csv' file.")

    
    # Load the transcript_gene_expression_df and medgraphica_transcript_protein.csv
    transcript_df = all_gene_dataframes['transcript_gene_expression_df']

    # Read the medgraphica_transcript_protein.csv file
    medgraphica_transcript_protein_df = pd.read_csv('medgraphica_data\Interaction\medgraphica_transcript_protein.csv')

    # Extract MedGraphica_IDs from both transcript_df and protein_df
    transcript_medgraphica_ids = transcript_df['MedGraphica_ID'].unique()
    protein_medgraphica_ids = protein_df['MedGraphica_ID'].unique()

    # Filter interactions in medgraphica_transcript_protein_df where Transcript or Protein matches MedGraphica_IDs
    transcript_protein_filtered_interactions = medgraphica_transcript_protein_df[
        (medgraphica_transcript_protein_df['Transcript'].isin(transcript_medgraphica_ids)) & 
        (medgraphica_transcript_protein_df['Protein'].isin(protein_medgraphica_ids))
    ]

    # Merge interaction data for Transcript and Protein
    # Map Transcript and Protein to their corresponding Gene_Name, then add '-Transcript' for transcripts and '-Protein' for proteins
    transcript_protein_interaction_with_gene_names = transcript_protein_filtered_interactions.copy()

    # Map Transcript to Gene_1 (from transcript_df)
    transcript_protein_interaction_with_gene_names = transcript_protein_interaction_with_gene_names.merge(transcript_df[['MedGraphica_ID', 'Gene_Name']], 
                                                                                                          left_on='Transcript', 
                                                                                                          right_on='MedGraphica_ID', 
                                                                                                          how='left')
    transcript_protein_interaction_with_gene_names.rename(columns={'Gene_Name': 'Gene_1'}, inplace=True)

    # Map Protein to Gene_2 (from protein_df)
    transcript_protein_interaction_with_gene_names = transcript_protein_interaction_with_gene_names.merge(protein_df[['MedGraphica_ID', 'Gene_Name']], 
                                                                                                          left_on='Protein', 
                                                                                                          right_on='MedGraphica_ID', 
                                                                                                          how='left')
    transcript_protein_interaction_with_gene_names.rename(columns={'Gene_Name': 'Gene_2'}, inplace=True)

    # Add '-Transcript' suffix to Gene_1 (for Transcript) and '-Protein' to Gene_2 (for Protein)
    transcript_protein_interaction_with_gene_names['Gene_1'] = transcript_protein_interaction_with_gene_names['Gene_1'] + '-Transcript'
    transcript_protein_interaction_with_gene_names['Gene_2'] = transcript_protein_interaction_with_gene_names['Gene_2'] + '-Protein'

    # Map Gene_1 (Transcript) and Gene_2 (Protein) to corresponding numeric IDs using sorted_all_gene_name_dict
    transcript_protein_interaction_with_gene_names['Gene_1_Num'] = transcript_protein_interaction_with_gene_names['Gene_1'].map(sorted_all_gene_name_dict)
    transcript_protein_interaction_with_gene_names['Gene_2_Num'] = transcript_protein_interaction_with_gene_names['Gene_2'].map(sorted_all_gene_name_dict)

    # Handle NaN values before converting to integers
    transcript_protein_interaction_with_gene_names['Gene_1_Num'] = transcript_protein_interaction_with_gene_names['Gene_1_Num'].fillna(0).astype(int)
    transcript_protein_interaction_with_gene_names['Gene_2_Num'] = transcript_protein_interaction_with_gene_names['Gene_2_Num'].fillna(0).astype(int)

    # Select relevant columns and add Gene_1_Name, Gene_2_Name, Gene_1_ID, Gene_2_ID
    transcript_protein_interaction_with_gene_names['Gene_1_ID'] = transcript_protein_interaction_with_gene_names['Transcript']
    transcript_protein_interaction_with_gene_names['Gene_2_ID'] = transcript_protein_interaction_with_gene_names['Protein']

    transcript_protein_interaction_with_gene_names['Edge_Type'] = 'Transcript-Protein'

    transcript_protein_interaction_with_gene_names.rename(columns={'Gene_1_Num': 'From'}, inplace=True)
    transcript_protein_interaction_with_gene_names.rename(columns={'Gene_2_Num': 'To'}, inplace=True)

    transcript_protein_edge_num_df = transcript_protein_interaction_with_gene_names[
        ['From', 'To', 'Edge_Type']
    ].dropna()

    # Save the data to transcript_protein_edge_num.csv file
    transcript_protein_edge_num_df.to_csv(os.path.join(graph_output_folder, 'transcript_protein_edge_num.csv'), index=False)

    print("Processing complete, the transcript-protein relationships have been saved to the 'transcript_protein_edge_num.csv' file.")

    # Load the gene_cnv_df and medgraphica_gene_transcript.csv
    gene_df = pd.concat([df for key, df in all_gene_dataframes.items() if 'cnv' in key.lower()])
    medgraphica_gene_transcript_df = pd.read_csv('medgraphica_data\Interaction\medgraphica_gene_transcript.csv')

    # Extract MedGraphica_IDs from both gene_df and transcript_df
    gene_medgraphica_ids = gene_df['MedGraphica_ID'].unique()
    transcript_medgraphica_ids = transcript_df['MedGraphica_ID'].unique()

    # Filter interactions in medgraphica_gene_transcript_df where Gene or Transcript matches MedGraphica_IDs
    gene_transcript_filtered_interactions = medgraphica_gene_transcript_df[
        (medgraphica_gene_transcript_df['Gene'].isin(gene_medgraphica_ids)) &
        (medgraphica_gene_transcript_df['Transcript'].isin(transcript_medgraphica_ids))
    ]

    # Merge interaction data for Gene and Transcript
    gene_transcript_filtered_interactions_with_gene_names = gene_transcript_filtered_interactions.copy()

    gene_transcript_filtered_interactions_with_gene_names = gene_transcript_filtered_interactions_with_gene_names.merge(gene_df[['MedGraphica_ID', 'Gene_Name']],
                                                                                                                        left_on='Gene',
                                                                                                                        right_on='MedGraphica_ID',
                                                                                                                        how='left')
    gene_transcript_filtered_interactions_with_gene_names.rename(columns={'Gene_Name': 'Gene_1'}, inplace=True)
    
    gene_transcript_filtered_interactions_with_gene_names = gene_transcript_filtered_interactions_with_gene_names.merge(transcript_df[['MedGraphica_ID', 'Gene_Name']],
                                                                                                                        left_on='Transcript',
                                                                                                                        right_on='MedGraphica_ID',
                                                                                                                        how='left')
    gene_transcript_filtered_interactions_with_gene_names.rename(columns={'Gene_Name': 'Gene_2'}, inplace=True)
    
    # Add '-Transcript' suffix to Gene_1 (for Transcript) and '-Protein' to Gene_2 (for Protein)
    gene_transcript_filtered_interactions_with_gene_names['Gene_1'] = gene_transcript_filtered_interactions_with_gene_names['Gene_1'] + '-Gene'
    gene_transcript_filtered_interactions_with_gene_names['Gene_2'] = gene_transcript_filtered_interactions_with_gene_names['Gene_2'] + '-Transcript'

    # Map Gene_1 (Transcript) and Gene_2 (Protein) to corresponding numeric IDs using sorted_all_gene_name_dict
    gene_transcript_filtered_interactions_with_gene_names['Gene_1_Num'] = gene_transcript_filtered_interactions_with_gene_names['Gene_1'].map(sorted_all_gene_name_dict)
    gene_transcript_filtered_interactions_with_gene_names['Gene_2_Num'] = gene_transcript_filtered_interactions_with_gene_names['Gene_2'].map(sorted_all_gene_name_dict)

    # Handle NaN values before converting to integers
    gene_transcript_filtered_interactions_with_gene_names['Gene_1_Num'] = gene_transcript_filtered_interactions_with_gene_names['Gene_1_Num'].fillna(0).astype(int)
    gene_transcript_filtered_interactions_with_gene_names['Gene_2_Num'] = gene_transcript_filtered_interactions_with_gene_names['Gene_2_Num'].fillna(0).astype(int)

    # Select relevant columns and add Gene_1_Name, Gene_2_Name, Gene_1_ID, Gene_2_ID
    gene_transcript_filtered_interactions_with_gene_names['Gene_1_ID'] = gene_transcript_filtered_interactions_with_gene_names['Gene']
    gene_transcript_filtered_interactions_with_gene_names['Gene_2_ID'] = gene_transcript_filtered_interactions_with_gene_names['Transcript']

    gene_transcript_filtered_interactions_with_gene_names['Edge_Type'] = 'Gene-Transcript'

    gene_transcript_filtered_interactions_with_gene_names.rename(columns={'Gene_1_Num': 'From'}, inplace=True)
    gene_transcript_filtered_interactions_with_gene_names.rename(columns={'Gene_2_Num': 'To'}, inplace=True)

    gene_transcript_edge_num_df = gene_transcript_filtered_interactions_with_gene_names[
        ['From', 'To', 'Edge_Type']
    ].dropna()

    # Save the data to gene_transcript_edge_num.csv file
    gene_transcript_edge_num_df.to_csv(os.path.join(graph_output_folder, 'gene_transcript_edge_num.csv'), index=False)

    print("Processing complete, the gene-transcript relationships have been saved to the 'gene_transcript_edge_num.csv' file.")

    # Load the promoter_df and medgraphica_promoter_gene.csv
    promoter_df = pd.concat([df for key, df in all_gene_dataframes.items() if 'methy' in key.lower()])
    medgraphica_promoter_gene_df = pd.read_csv('medgraphica_data\Interaction\medgraphica_promoter_gene.csv')

    # Extract MedGraphica_IDs from both promoter_df and gene_df
    promoter_medgraphica_ids = promoter_df['MedGraphica_ID'].unique()
    gene_medgraphica_ids = gene_df['MedGraphica_ID'].unique()

    # Filter interactions in medgraphica_promoter_gene_df where Promoter or Gene matches MedGraphica_IDs
    promoter_gene_filtered_interactions = medgraphica_promoter_gene_df[
        (medgraphica_promoter_gene_df['From_ID'].isin(promoter_medgraphica_ids)) &
        (medgraphica_promoter_gene_df['To_ID'].isin(gene_medgraphica_ids))
    ]

    # Merge interaction data for Promoter and Gene
    promoter_gene_filtered_interactions_with_gene_names = promoter_gene_filtered_interactions.copy()

    promoter_gene_filtered_interactions_with_gene_names = promoter_gene_filtered_interactions_with_gene_names.merge(promoter_df[['MedGraphica_ID', 'Gene_Name']],
                                                                                                                        left_on='From_ID',
                                                                                                                        right_on='MedGraphica_ID',
                                                                                                                        how='left')
    promoter_gene_filtered_interactions_with_gene_names.rename(columns={'Gene_Name': 'From_Gene'}, inplace=True)
    
    promoter_gene_filtered_interactions_with_gene_names = promoter_gene_filtered_interactions_with_gene_names.merge(gene_df[['MedGraphica_ID', 'Gene_Name']],
                                                                                                                        left_on='To_ID',
                                                                                                                        right_on='MedGraphica_ID',
                                                                                                                        how='left')
    promoter_gene_filtered_interactions_with_gene_names.rename(columns={'Gene_Name': 'To_Gene'}, inplace=True)
    
    # Add suffix
    promoter_gene_filtered_interactions_with_gene_names['From_Gene'] = promoter_gene_filtered_interactions_with_gene_names['From_Gene'] + '-Promoter'
    promoter_gene_filtered_interactions_with_gene_names['To_Gene'] = promoter_gene_filtered_interactions_with_gene_names['To_Gene'] + '-Gene'

    # Map Gene_1 (Transcript) and Gene_2 (Protein) to corresponding numeric IDs using sorted_all_gene_name_dict
    promoter_gene_filtered_interactions_with_gene_names['From_Gene_Num'] = promoter_gene_filtered_interactions_with_gene_names['From_Gene'].map(sorted_all_gene_name_dict)
    promoter_gene_filtered_interactions_with_gene_names['To_Gene_Num'] = promoter_gene_filtered_interactions_with_gene_names['To_Gene'].map(sorted_all_gene_name_dict)

    # Handle NaN values before converting to integers
    promoter_gene_filtered_interactions_with_gene_names['From_Gene_Num'] = promoter_gene_filtered_interactions_with_gene_names['From_Gene_Num'].fillna(0).astype(int)
    promoter_gene_filtered_interactions_with_gene_names['From_Gene_Num'] = promoter_gene_filtered_interactions_with_gene_names['To_Gene_Num'].fillna(0).astype(int)

    # Select relevant columns and add Gene_1_Name, Gene_2_Name, Gene_1_ID, Gene_2_ID
    promoter_gene_filtered_interactions_with_gene_names['From_Gene_ID'] = promoter_gene_filtered_interactions_with_gene_names['From_ID']
    promoter_gene_filtered_interactions_with_gene_names['From_Gene_ID'] = promoter_gene_filtered_interactions_with_gene_names['To_ID']

    promoter_gene_filtered_interactions_with_gene_names['Edge_Type'] = 'Promoter-Gene'

    promoter_gene_filtered_interactions_with_gene_names.rename(columns={'From_Gene_Num': 'From'}, inplace=True)
    promoter_gene_filtered_interactions_with_gene_names.rename(columns={'To_Gene_Num': 'To'}, inplace=True)

    promoter_gene_edge_num_df = promoter_gene_filtered_interactions_with_gene_names[
        ['From', 'To', 'Edge_Type']
    ].dropna()

    # Save the data to promoter_gene_edge_num.csv file
    promoter_gene_edge_num_df.to_csv(os.path.join(graph_output_folder, 'promoter_gene_edge_num.csv'), index=False)

    print("Processing complete, the promoter-gene relationships have been saved to the 'promoter_gene_edge_num.csv' file.")

    all_gene_edge_num_df = pd.concat([protein_protein_edge_num_df, transcript_protein_edge_num_df, gene_transcript_edge_num_df, promoter_gene_edge_num_df])
    # num_protein_protein_edge = protein_protein_edge_num_df.shape[0]
    # num_transcript_protein_edge = transcript_protein_edge_num_df.shape[0]
    # num_gene_transcript_edge = gene_transcript_edge_num_df.shape[0]
    # num_promoter_gene_edge = promoter_gene_edge_num_df.shape[0]

    # edge_type_list = ['Protein-Protein'] * num_protein_protein_edge + ['Transcript-Protein'] * num_transcript_protein_edge + ['Gene-Transcript'] * num_gene_transcript_edge + ['Promoter-Gene'] * num_promoter_gene_edge
    # all_gene_edge_num_df['Edge_Type'] = edge_type_list
    # all_gene_edge_num_df = all_gene_edge_num_df.sort_values(by=['From', 'To']).reset_index(drop=True)

    all_gene_edge_num_df.to_csv(os.path.join(graph_output_folder, 'all_gene_edge_num.csv'), index=False)


    # --------------------------------------------------------------------------------------------------------------------------------------
    # Survival data processing
    survival_filtered = survival[survival['individualID'].isin(sorted_patient_id_list)]
    
    # Calculate the proportion of NaN values in each column
    survival_nan_column_proportions = survival_filtered.isna().mean()

    # Identify columns to be dropped (where proportion of NaN values is greater than 1/3)
    columns_to_drop = survival_nan_column_proportions[survival_nan_column_proportions > 1/3].index.tolist()

    # Drop these columns from the DataFrame
    survival_filtered = survival_filtered.drop(columns=columns_to_drop)

    # List of columns that were dropped
    print("Columns dropped:", columns_to_drop)

    # Change it to binary
    def modify_ceradsc(value):
        if value in [1, 2]:
            return 0
        elif value in [3, 4]:
            return 1
        else:
            return value  # Keeps other values as they are, if there are any

    survival_filtered['ceradsc'] = survival_filtered['ceradsc'].apply(modify_ceradsc)
    
    count_female = (survival_filtered['msex'] == 0).sum()
    count_male = (survival_filtered['msex'] == 1).sum()
    print(count_female, count_male)

    count_AD = ((survival_filtered['ceradsc'] == 0)).sum()
    count_NOAD = ((survival_filtered['ceradsc'] == 1)).sum()
    count_AD, count_NOAD
    print(count_AD, count_NOAD)

    count_female_AD = ((survival_filtered['msex'] == 0) & ((survival_filtered['ceradsc'] == 0))).sum()
    count_male_AD = ((survival_filtered['msex'] == 1) & ((survival_filtered['ceradsc'] == 0))).sum()
    print('AD Female:', count_female_AD)
    print('AD Male:', count_male_AD)

    from sklearn.utils import resample

    # Calculate the number of AD and non-AD samples in the original dataset
    count_AD = (((survival_filtered['ceradsc'] == 0)) ).sum()
    count_NOAD = ((survival_filtered['ceradsc'] == 1) ).sum()

    # Get the AD and non-AD datasets
    df_AD = survival_filtered[((survival_filtered['ceradsc'] == 0))]
    df_NOAD = survival_filtered[(survival_filtered['ceradsc'] == 1)]

    # Determine the number of samples after downsampling, taking the smaller value between female and male AD sample counts
    n_samples = min(count_AD, count_NOAD)

    # Downsample the AD dataset
    df_AD_downsampled = resample(df_AD, 
                                        replace=False,
                                        n_samples=n_samples,
                                        random_state=123)

    # Downsample the non-AD dataset
    df_NOAD_downsampled = resample(df_NOAD, 
                                    replace=False,
                                    n_samples=n_samples,
                                    random_state=123)

    # Combine the downsampled AD and non-AD datasets
    df_balanced_ds = pd.concat([df_AD_downsampled, df_NOAD_downsampled]).reset_index(drop=True)

    # Print the counts of the downsampled female and male AD samples
    count_AD_downsampled = len(df_AD_downsampled)
    count_NOAD_downsampled = len(df_NOAD_downsampled)
    print('AD:', count_AD_downsampled)
    print('non-AD:', count_NOAD_downsampled)

    # Overwrite the original survival_filtered with the balanced dataset
    survival_filtered = df_balanced_ds
    print(survival_filtered)

    # Filter the dataset to keep only AD samples
    df_AD = survival_filtered[((survival_filtered['ceradsc'] == 0))]
    df_NOAD = survival_filtered[(survival_filtered['ceradsc'] == 1)]

    # Calculate the number of female and male AD samples
    count_female_AD = (df_AD['msex'] == 0).sum()
    count_male_AD = (df_AD['msex'] == 1).sum()

    # Print the counts of female and male AD samples
    print('AD Female:', count_female_AD)
    print('AD Male:', count_male_AD)

    # Calculate the number of female and male NOAD samples
    count_female_NOAD = (df_NOAD['msex'] == 0).sum()
    count_male_NOAD = (df_NOAD['msex'] == 1).sum()

    # Print the counts of female and male NOAD samples
    print('non-AD Female:',  count_female_NOAD)
    print('non-AD Male:', count_male_NOAD)

    # Overwrite the original survival_filtered with the balanced dataset
    survival_filtered = df_balanced_ds.reset_index(drop=True)
    print(survival_filtered)


    survival_filtered_feature_df = survival_filtered.copy()

    # TODO
    survival_filtered_feature_df = survival_filtered_feature_df[['individualID', 'ceradsc']]
    # survival_filtered_feature_df = survival_filtered_feature_df[['individualID', 'msex']]
    print(survival_filtered_feature_df)

    nan_counts = survival_filtered_feature_df.isna().sum()  # or df.isnull()

    # TODO
    # print(survival_filtered_feature_df['ceradsc'].unique())
    # print(survival_filtered_feature_df['msex'].unique())

    survival_filtered_feature_df.to_csv(os.path.join(graph_output_folder, 'survival_label.csv'), index=False)

    # Randomize the survival label
    def input_random(randomized, graph_output_folder):
        if randomized == True:
            random_survival_filtered_feature_df = survival_filtered_feature_df.sample(frac = 1).reset_index(drop=True)
            random_survival_filtered_feature_df.to_csv(os.path.join(graph_output_folder, 'random_survival_label.csv'), index=False)
        else:
            random_survival_filtered_feature_df = pd.read_csv(os.path.join(graph_output_folder, 'random_survival_label.csv'))
        print(random_survival_filtered_feature_df)

    # TODO
    input_random(randomized=True, graph_output_folder=graph_output_folder)
    # input_random(randomized=False, graph_output_folder=graph_output_folder)
    
    # Split deep learning input into training and test
    def split_k_fold(k, graph_output_folder):
        random_survival_filtered_feature_df = pd.read_csv(os.path.join(graph_output_folder, 'random_survival_label.csv'))
        num_points = random_survival_filtered_feature_df.shape[0]
        num_div = int(num_points / k)
        num_div_list = [i * num_div for i in range(0, k)]
        num_div_list.append(num_points)
        # Split [random_survival_filtered_feature_df] into [k] folds
        for place_num in range(k):
            low_idx = num_div_list[place_num]
            high_idx = num_div_list[place_num + 1]
            print('\n--------TRAIN-TEST SPLIT WITH TEST FROM ' + str(low_idx) + ' TO ' + str(high_idx) + '--------')
            split_input_df = random_survival_filtered_feature_df[low_idx : high_idx]
            split_input_df.to_csv(os.path.join(graph_output_folder, 'split_random_survival_label_' + str(place_num + 1) + '.csv'), index=False)
            print(split_input_df.shape)

    split_k_fold(k=5, graph_output_folder=graph_output_folder)

    # Process the edge_index file

    gene_edge_num_df = pd.read_csv(os.path.join(graph_output_folder, 'all_gene_edge_num.csv'))
    src_gene_list = list(gene_edge_num_df['From'])
    dest_gene_list = list(gene_edge_num_df['To'])
    edge_type_list = list(gene_edge_num_df['Edge_Type'])
    gene_edge_num_reverse_df = pd.DataFrame({'From': dest_gene_list, 'To': src_gene_list, 'Edge_Type': edge_type_list})
    # gene_edge_num_reverse_wointernal_df = gene_edge_num_reverse_df[~gene_edge_num_reverse_df['Path'].str.contains('internal link')]
    print(gene_edge_num_df)
    gene_edge_num_all_df = gene_edge_num_df.drop_duplicates().sort_values(by=['From', 'To']).reset_index(drop=True)
    print(gene_edge_num_all_df)
    gene_edge_num_all_df.to_csv(os.path.join(graph_output_folder, 'gene_edge_num_all_df.csv'), index=False)

    gene_edge_name_all_df = gene_edge_num_all_df.replace(sorted_all_gene_dict)
    print(gene_edge_name_all_df)
    gene_edge_name_all_df.to_csv(os.path.join(graph_output_folder, 'gene_edge_name_all_df.csv'), index=False)

    # remove internal link
    kegg_path_gene_edge_num_all_df = gene_edge_num_all_df
    kegg_path_gene_edge_num_all_df.to_csv(os.path.join(graph_output_folder, 'path_gene_edge_num_all.csv'), index=False)
    # keep only internal link
    internal_gene_edge_num_all_df = gene_edge_num_all_df

    internal_gene_edge_num_all_df.to_csv(os.path.join(graph_output_folder, 'internal_gene_edge_num_all.csv'), index=False)

    kegg_path_gene_edge_name_all_df = kegg_path_gene_edge_num_all_df.replace(sorted_all_gene_dict)
    kegg_path_gene_edge_name_all_df.to_csv(os.path.join(graph_output_folder, 'path_gene_edge_name_all.csv'), index=False)
    internal_gene_edge_name_all_df = internal_gene_edge_num_all_df.replace(sorted_all_gene_dict)
    internal_gene_edge_name_all_df.to_csv(os.path.join(graph_output_folder, 'internal_gene_edge_name_all.csv'), index=False)