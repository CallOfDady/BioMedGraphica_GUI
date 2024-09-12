# data_process.py

import pandas as pd
import numpy as np

def process_gene(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Gene data."""
    print(f"Processing Gene: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}, Feature Label: {feature_label}")
    
    # Load the gene entity data
    gene_entity_data = pd.read_csv('medgraphica_data/medgraphica_gene.csv')    
        
    # Determine the separator based on the file extension
    if file_path.endswith('.txt') or file_path.endswith('.tsv'):
        sep = '\t'
    elif file_path.endswith('.csv'):
        sep = ','
    else:
        raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")

    def process_gene_expression(gene_entity_data):
        """Process Gene Expression data."""
        print(f"Processing Gene Expression for {entity_type} with Feature Label: {feature_label}")

        # Read the gene expression file
        gene_expression = pd.read_csv(file_path, sep=sep)

        # Process gene_id by stripping the version number
        gene_ids = gene_expression[selected_column].apply(lambda x: x.split('.')[0]).tolist()
        gene_expression[selected_column] = gene_ids

        # Retain only the necessary columns and rename them
        gene_entity_data = gene_entity_data[[id_type, 'medgraphica_id', 'hgnc_symbol']].copy()
        gene_entity_data.rename(columns={id_type: 'gene_id', 'hgnc_symbol': 'gene_name'}, inplace=True)

        # # Drop NaN and duplicates from gene_entity_data
        # gene_entity_data = gene_entity_data.dropna(subset=['medgraphica_id']).drop_duplicates()
        # print("gene_entity_data shape after dropna and drop_duplicates:", gene_entity_data.shape)

        # Merge the gene expression data with the transcript entity data on 'gene_id' and keep 'medgraphica_id' and 'gene_name'
        merged_data = pd.merge(gene_expression, gene_entity_data, on='gene_id', how='inner')

        # Replace 'Nan' strings with actual NaN values
        merged_data.replace('Nan', pd.NA, inplace=True)

        # Drop rows where 'medgraphica_id' is NaN (if any)
        merged_data.dropna(subset=['medgraphica_id'], inplace=True)

        # Ensure that 'medgraphica_id' and 'gene_name' are the first columns, followed by the rest in their original order
        cols = ['medgraphica_id', 'gene_name'] + [col for col in merged_data.columns if col not in ['medgraphica_id', 'gene_name']]
        gene_expression_new = merged_data[cols]
        
        # Select numeric columns
        numeric_cols = gene_expression_new.select_dtypes(include=[np.number]).columns.tolist()

        # Group by 'medgraphica_id' and calculate the mean for all numeric columns
        gene_expression_new = gene_expression_new.groupby('medgraphica_id')[numeric_cols].mean()

        # Reset the index to turn 'medgraphica_id' back into a column
        gene_expression_new.reset_index(inplace=True)

        # Merge back the gene names after grouping by 'medgraphica_id'
        gene_expression_new = pd.merge(gene_expression_new, gene_entity_data[['medgraphica_id', 'gene_name']], on='medgraphica_id', how='left')

        # Ensure 'medgraphica_id' is first, 'gene_name' is second, and other columns follow the original order
        final_cols = ['medgraphica_id', 'gene_name'] + [col for col in gene_expression_new.columns if col not in ['medgraphica_id', 'gene_name']]
        gene_expression_new = gene_expression_new[final_cols]

        # Export the result to a CSV file
        output_file_path = f'cache/{feature_label.lower()}_{entity_type.lower()}.csv'
        gene_expression_new.to_csv(output_file_path, sep=",", index=False)

        print(f"Gene expression data processing completed. Output saved to {output_file_path}")

    def process_gene_cnv(gene_entity_data):
        """Process Gene CNV data."""
        print(f"Processing Gene CNV for {entity_type} with Feature Label: {feature_label}")
        
        def get_gene_names(chrom, start, end):
            """Retrieve gene names and medgraphica_ids if the given chrom, start, and end positions are within the gene's chromosome, start, and end positions."""
            
            # Convert 'Gene start (bp)', 'Gene end (bp)', and 'chromosome' columns to numeric
            gene_entity_data['gene_start'] = pd.to_numeric(gene_entity_data['gene_start'], errors='coerce')
            gene_entity_data['gene_end'] = pd.to_numeric(gene_entity_data['gene_end'], errors='coerce')
            gene_entity_data['chromosome'] = pd.to_numeric(gene_entity_data['chromosome'], errors='coerce')

            # Filter rows by chromosome first
            matching_genes = gene_entity_data[gene_entity_data['chromosome'] == chrom]
            
            # Further filter rows where the start and end positions fall within the gene's start and end positions
            matching_genes = matching_genes[
                (matching_genes['gene_start'] <= start) & 
                (matching_genes['gene_end'] >= end)
            ]

            # Extract the HGNC symbol (gene name) and medgraphica_id from the matching genes
            gene_names = matching_genes['hgnc_symbol'].tolist()
            medgraphica_ids = matching_genes['medgraphica_id'].tolist()

            # Ensure all gene names and medgraphica_ids are strings
            gene_names = [str(gene) for gene in gene_names]
            medgraphica_ids = [str(id) for id in medgraphica_ids]

            # Return both gene names and medgraphica_ids
            return ', '.join(medgraphica_ids) if medgraphica_ids else None, ', '.join(gene_names) if gene_names else None
        
        # Read the CNV file
        cnv = pd.read_csv(file_path, sep=sep)
        print(f"CNV data loaded. Columns available: {cnv.columns.tolist()}")

        # Split the selected_column (e.g., "START, END") into two column names
        start_col, end_col = selected_column.split(', ')
        print(f"Selected columns for START: {start_col}, END: {end_col}")
        
        # Find the chromosome column in the CNV DataFrame
        chrom_col = [col for col in cnv.columns if 'CHROM' in col.upper()]
        if not chrom_col:
            raise ValueError("No chromosome column found with 'CHROM' in the name.")
        chrom_col = chrom_col[0]
        print(f"Chromosome column identified: {chrom_col}")
        
        # Initialize empty lists to store gene names and medgraphica_ids
        gene_names = []
        medgraphica_ids = []

        # Iterate over each row in the CNV DataFrame
        for index, row in cnv.iterrows():
            chrom = row[chrom_col]   # Get the chromosome value from the specified chromosome column
            start = row[start_col]   # Get the start value from the specified start column
            end = row[end_col]       # Get the end value from the specified end column
            
            # print(f"Processing row {index}: CHROM = {chrom}, START = {start}, END = {end}")
            
            # Call the get_gene_names function to retrieve the medgraphica_id and gene_name for the current chrom, start, and end range
            try:
                medgraphica_id, gene_name = get_gene_names(chrom, start, end)
                # print(f"Medgraphica ID(s) found: {medgraphica_id}, Gene name(s) found: {gene_name}")
            except Exception as e:
                print(f"Error during gene name lookup for row {index}: {e}")
                # Set only missing values to None (i.e., handle partial matches)
                if 'medgraphica_id' not in locals():
                    medgraphica_id = None
                if 'gene_name' not in locals():
                    gene_name = None
            
            # Ensure medgraphica_id and gene_name are strings; if None, convert them to empty strings
            medgraphica_ids.append(medgraphica_id if medgraphica_id is not None else "")
            gene_names.append(gene_name if gene_name is not None else "")
        
        print("Finished processing all rows. Adding medgraphica_id and gene_name columns to CNV data.")
        
        # Create a copy of CNV data and add the medgraphica_id and gene_name columns
        cnv_data = cnv.copy()
        cnv_data['medgraphica_id'] = medgraphica_ids
        cnv_data['gene_name'] = gene_names

        # Reorder columns to place 'medgraphica_id' and 'gene_name' at the front
        cols = ['medgraphica_id', 'gene_name'] + [col for col in cnv_data.columns if col not in ['medgraphica_id', 'gene_name']]
        cnv_data = cnv_data[cols]

        # Split 'medgraphica_id' and 'gene_name' by comma and ensure they're synchronized
        cnv_data['medgraphica_id'] = cnv_data['medgraphica_id'].str.split(', ')
        cnv_data['gene_name'] = cnv_data['gene_name'].str.split(', ')

        # Check if both columns have the same number of elements before exploding
        assert cnv_data['medgraphica_id'].apply(len).equals(cnv_data['gene_name'].apply(len)), "Mismatch in number of elements between 'medgraphica_id' and 'gene_name'"

        # Explode both columns so that each pair of medgraphica_id and gene_name is on its own row
        cnv_data_exploded = cnv_data.explode(['medgraphica_id', 'gene_name'])

        # Strip whitespace to ensure clean gene names and medgraphica_ids
        cnv_data_exploded['gene_name'] = cnv_data_exploded['gene_name'].str.strip()
        cnv_data_exploded['medgraphica_id'] = cnv_data_exploded['medgraphica_id'].str.strip()


        # Filter out rows where 'medgraphica_id' or 'gene_name' is empty
        cnv_data_filtered = cnv_data_exploded[(cnv_data_exploded['medgraphica_id'] != '') & (cnv_data_exploded['gene_name'] != '')]
        sample_columns = cnv_data_filtered.columns[8:]
        # print(f"Sample columns: {sample_columns}")
        cnv_aggregated = cnv_data_filtered.groupby(['medgraphica_id', 'gene_name', 'SVTYPE'])[sample_columns].sum().reset_index()
        
        # Split the CNV data into DEL, DUP, and mCNV
        cnv_aggregated_del = cnv_aggregated[cnv_aggregated['SVTYPE'] == 'DEL']
        cnv_aggregated_dup = cnv_aggregated[cnv_aggregated['SVTYPE'] == 'DUP']
        cnv_aggregated_mcnv = cnv_aggregated[cnv_aggregated['SVTYPE'] == 'mCNV']
        print(f"Shape of cnv_aggregated_del: {cnv_aggregated_del.shape}")
        print(f"Shape of cnv_aggregated_dup: {cnv_aggregated_dup.shape}")
        print(f"Shape of cnv_aggregated_mcnv: {cnv_aggregated_mcnv.shape}")


        # Create a DataFrame from the union set of all gene names
        cnv_all_genes = pd.DataFrame({'medgraphica_id': list(set(cnv_aggregated_del['medgraphica_id']).union(set(cnv_aggregated_dup['medgraphica_id']), set(cnv_aggregated_mcnv['medgraphica_id'])))})

        # Merge each aggregated DataFrame with the all genes DataFrame
        # This will align all genes across all DataFrames, and fill missing entries with 0
        cnv_aggregated_del = pd.merge(cnv_all_genes, cnv_aggregated_del, on='medgraphica_id', how='outer').fillna(-1)
        cnv_aggregated_dup = pd.merge(cnv_all_genes, cnv_aggregated_dup, on='medgraphica_id', how='outer').fillna(-1)
        cnv_aggregated_mcnv = pd.merge(cnv_all_genes, cnv_aggregated_mcnv, on='medgraphica_id', how='outer').fillna(-1)

        # Drop the 'SVTYPE' column from each merged DataFrame
        cnv_aggregated_del.drop(columns='SVTYPE', inplace=True)
        cnv_aggregated_dup.drop(columns='SVTYPE', inplace=True)
        cnv_aggregated_mcnv.drop(columns='SVTYPE', inplace=True)
        print(f"Shape of cnv_aggregated_del: {cnv_aggregated_del.shape}")
        print(f"Shape of cnv_aggregated_dup: {cnv_aggregated_dup.shape}")
        print(f"Shape of cnv_aggregated_mcnv: {cnv_aggregated_mcnv.shape}")
    
        output_file_path = f'cache/{feature_label.lower()}_{entity_type.lower()}.csv'
        cnv_aggregated.to_csv(output_file_path, sep=",", index=False)
        cnv_aggregated_del.to_csv(f'cache/{feature_label.lower()}_{entity_type.lower()}_del.csv', sep=",", index=False)
        cnv_aggregated_dup.to_csv(f'cache/{feature_label.lower()}_{entity_type.lower()}_dup.csv', sep=",", index=False)
        cnv_aggregated_mcnv.to_csv(f'cache/{feature_label.lower()}_{entity_type.lower()}_mcnv.csv', sep=",", index=False)
        
        print(f"CNV data processing completed. Output saved to {output_file_path}")
        
        return cnv_data

    # Check if the feature_label contains "expression" or "cnv"
    if "expression" in feature_label.lower():
        process_gene_expression(gene_entity_data)
    elif "cnv" in feature_label.lower():
        process_gene_cnv(gene_entity_data)
    else:
        print(f"Processing for Feature Label: {feature_label} is not implemented yet.")
        # Leave space for other feature_label handling in the future

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

        # Read the gene expression file
        gene_expression = pd.read_csv(file_path, sep=sep)
        print("Original gene_expression shape:", gene_expression.shape)

        # Process gene_id by stripping the version number
        gene_ids = gene_expression[selected_column].apply(lambda x: x.split('.')[0]).tolist()
        gene_expression[selected_column] = gene_ids

        # Retain only the necessary columns and rename them
        transcript_entity_data = transcript_entity_data[[id_type, 'medgraphica_id', 'gene_name']].copy()
        transcript_entity_data.rename(columns={id_type: 'gene_id'}, inplace=True)
        print("transcript_entity_data shape after selection and rename:", transcript_entity_data.shape)

        # Drop NaN and duplicates from transcript_entity_data
        transcript_entity_data = transcript_entity_data.dropna(subset=['medgraphica_id']).drop_duplicates()
        print("transcript_entity_data shape after dropna and drop_duplicates:", transcript_entity_data.shape)

        # Merge the gene expression data with the transcript entity data on 'gene_id' and keep 'medgraphica_id' and 'gene_name'
        merged_data = pd.merge(gene_expression, transcript_entity_data, on='gene_id', how='inner')
        print("Merged data shape:", merged_data.shape)

        # Replace 'Nan' strings with actual NaN values
        merged_data.replace('Nan', pd.NA, inplace=True)

        # Drop rows where 'medgraphica_id' is NaN (if any)
        merged_data.dropna(subset=['medgraphica_id'], inplace=True)

        # Ensure that 'medgraphica_id' and 'gene_name' are the first columns, followed by the rest in their original order
        cols = ['medgraphica_id', 'gene_name'] + [col for col in merged_data.columns if col not in ['medgraphica_id', 'gene_name']]
        gene_expression_new = merged_data[cols]
        print("gene_expression_new shape after column reordering:", gene_expression_new.shape)

        # Select numeric columns
        numeric_cols = gene_expression_new.select_dtypes(include=[np.number]).columns.tolist()

        # Group by 'medgraphica_id' and calculate the mean for all numeric columns
        gene_expression_new = gene_expression_new.groupby('medgraphica_id')[numeric_cols].mean()
        print("gene_expression_new shape after grouping by 'medgraphica_id':", gene_expression_new.shape)

        # Reset the index to turn 'medgraphica_id' back into a column
        gene_expression_new.reset_index(inplace=True)

        # Merge back the gene names after grouping by 'medgraphica_id'
        gene_expression_new = pd.merge(gene_expression_new, transcript_entity_data[['medgraphica_id', 'gene_name']], on='medgraphica_id', how='left')
        print("Final gene_expression_new shape after merging back gene_name:", gene_expression_new.shape)

        # Ensure 'medgraphica_id' is first, 'gene_name' is second, and other columns follow the original order
        final_cols = ['medgraphica_id', 'gene_name'] + [col for col in gene_expression_new.columns if col not in ['medgraphica_id', 'gene_name']]
        gene_expression_new = gene_expression_new[final_cols]

        # Export the result to a CSV file
        output_file_path = f'cache/{feature_label.lower()}_{entity_type.lower()}.csv'
        gene_expression_new.to_csv(output_file_path, sep=",", index=False)


        print(f"Gene expression data processing completed. Output saved to {output_file_path}")

    # Check if the feature_label contains "expression" or "cnv"
    if "expression" in feature_label.lower():
        process_gene_expression(transcript_entity_data)
    elif "xxx" in feature_label.lower():
        print(f"Processing xxx data for Transcript: {entity_type} with Feature Label: {feature_label}")
    else:
        print(f"Processing for Feature Label: {feature_label} is not implemented yet.")
        # Leave space for other feature_label handling in the future

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
    protein.rename(columns={'Unnamed: 0': 'id_column'}, inplace=True)
    selected_column
    
    # Split the selected_column into protein_names and protein_ids
    protein['protein_name'] = protein[selected_column].apply(lambda x: x.split('|')[0])  # Extract protein names
    protein['protein_id'] = protein[selected_column].apply(lambda x: x.split('|')[1])    # Extract protein ids
    
    # Merge with protein_entity_data on protein_id and selected_column
    # Keep only the 'individualID' column from protein_entity_data
    protein_merged = pd.merge(protein, 
                              protein_entity_data[[selected_column, 'individualID']], 
                              left_on='protein_id', 
                              right_on=selected_column, 
                              how='left')
    
    # Drop the extra merge key (selected_column from protein_entity_data)
    protein_merged.drop(columns=[selected_column], inplace=True)

    # Fill any missing individualID values with NaN
    protein_merged['individualID'].fillna('NaN', inplace=True)

    # Calculate the NaN proportion of each row
    nan_proportions = protein.isna().mean(axis=1)
    print("NaN proportions by row:")
    print(nan_proportions)
    
    # Group by 'gene_name' and calculate the mean of each numeric column
    protein_grouped = protein_merged.groupby('protein_name', as_index=False).mean()
    
    # Fill NaN values with 0
    protein_grouped = protein_grouped.fillna(0)

    # Export the result to a CSV file
    output_file_path = f'cache/{feature_label.lower()}_{entity_type.lower()}.csv'
    protein_grouped.to_csv(output_file_path, sep=",", index=False)

    print(f"Protein data processing completed. Output saved to {output_file_path}")

    return protein_grouped

def process_drug(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Drug data."""
    print(f"Processing Drug: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}")
    # TODO

def process_disease(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Disease data."""
    print(f"Processing Disease: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}")
    # TODO

def process_phenotype(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Phenotype data."""
    print(f"Processing Phenotype: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}")
    # TODO
