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
        gene_expression[selected_column] = gene_expression[selected_column].apply(lambda x: x.split('.')[0])

        # Retain only the necessary columns
        gene_entity_data = gene_entity_data[[id_type, 'MedGraphica_ID', 'HGNC_Symbol']].copy()

        # Merge the gene expression data with the gene entity data
        merged_data = pd.merge(
            gene_expression,
            gene_entity_data,
            left_on=selected_column,
            right_on=id_type,
            how='inner'
        )

        # Replace 'Nan' strings with actual NaN values
        merged_data.replace('Nan', pd.NA, inplace=True)

        # Drop rows where 'MedGraphica_ID' is NaN (if any)
        merged_data.dropna(subset=['MedGraphica_ID'], inplace=True)

        # Ensure 'MedGraphica_ID' and 'HGNC_Symbol' are the first and second columns, respectively
        cols = ['MedGraphica_ID', 'HGNC_Symbol'] + [col for col in merged_data.columns if col not in ['MedGraphica_ID', 'HGNC_Symbol', selected_column, id_type]]
        gene_expression_new = merged_data[cols]
        
        # Select numeric columns
        numeric_cols = gene_expression_new.select_dtypes(include=[np.number]).columns.tolist()

        # Group by 'MedGraphica_ID' and calculate the mean for all numeric columns
        gene_expression_new = gene_expression_new.groupby('MedGraphica_ID')[numeric_cols].mean()

        # Reset the index to turn 'MedGraphica_ID' back into a column
        gene_expression_new.reset_index(inplace=True)

        # Merge back the gene names after grouping by 'MedGraphica_ID'
        gene_expression_new = pd.merge(gene_expression_new, gene_entity_data[['MedGraphica_ID', 'HGNC_Symbol']], on='MedGraphica_ID', how='left')

        # Rename 'HGNC_Symbol' to 'Gene_Name'
        gene_expression_new.rename(columns={'HGNC_Symbol': 'Gene_Name'}, inplace=True)

        # Ensure 'MedGraphica_ID' is first, 'Gene_Name' is second, and other columns follow the original order
        final_cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in gene_expression_new.columns if col not in ['MedGraphica_ID', 'Gene_Name']]
        gene_expression_new = gene_expression_new[final_cols]

        # Export the result to a CSV file
        output_file_path = f'cache/{entity_type.lower()}_{feature_label.lower()}.csv'
        gene_expression_new.to_csv(output_file_path, sep=",", index=False)

        print(f"Gene expression data processing completed. Output saved to {output_file_path}")

    def process_gene_cnv(gene_entity_data):
        """Process Gene CNV data."""
        print(f"Processing Gene CNV for {entity_type} with Feature Label: {feature_label}")
        
        def get_gene_names(chrom, start, end):
            """Retrieve gene names and medgraphica_ids if the given chrom, start, and end positions are within the gene's chromosome, start, and end positions."""
            
            # Convert 'Gene start (bp)', 'Gene end (bp)', and 'chromosome' columns to numeric
            gene_entity_data['Gene_Start'] = pd.to_numeric(gene_entity_data['Gene_Start'], errors='coerce')
            gene_entity_data['Gene_End'] = pd.to_numeric(gene_entity_data['Gene_End'], errors='coerce')
            gene_entity_data['Chromosome'] = pd.to_numeric(gene_entity_data['Chromosome'], errors='coerce')

            # Filter rows by chromosome first
            matching_genes = gene_entity_data[gene_entity_data['Chromosome'] == chrom]
            
            # Further filter rows where the start and end positions fall within the gene's start and end positions
            matching_genes = matching_genes[
                (matching_genes['Gene_Start'] <= start) & 
                (matching_genes['Gene_End'] >= end)
            ]

            # Extract the HGNC symbol (gene name) and MedGraphica_ID from the matching genes
            gene_names = matching_genes['HGNC_Symbol'].tolist()
            medgraphica_ids = matching_genes['MedGraphica_ID'].tolist()

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
            
            # Call the get_gene_names function to retrieve the MedGraphica_ID and Gene_Name for the current chrom, start, and end range
            try:
                medGraphica_id, gene_name = get_gene_names(chrom, start, end)
            except Exception as e:
                print(f"Error during gene name lookup for row {index}: {e}")
                # Set only missing values to None (i.e., handle partial matches)
                if 'MedGraphica_ID' not in locals():
                    medGraphica_id = None
                if 'Gene_Name' not in locals():
                    gene_name = None
            
            # Ensure MedGraphica_ID and gene_name are strings; if None, convert them to empty strings
            medgraphica_ids.append(medGraphica_id if medGraphica_id is not None else "")
            gene_names.append(gene_name if gene_name is not None else "")
        
        print("Finished processing all rows. Adding MedGraphica_ID and gene_name columns to CNV data.")
        
        # Create a copy of CNV data and add the MedGraphica_ID and gene_name columns
        cnv_data = cnv.copy()
        cnv_data['MedGraphica_ID'] = medgraphica_ids
        cnv_data['Gene_Name'] = gene_names

        # Reorder columns to place 'MedGraphica_ID' and 'Gene_Name' at the front
        cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in cnv_data.columns if col not in ['MedGraphica_ID', 'Gene_Name']]
        cnv_data = cnv_data[cols]
        # cnv_data.to_csv('cache/cnv_data.csv', index=False)

        # Split 'MedGraphica_ID' and 'Gene_Name' by comma and ensure they're synchronized
        cnv_data['MedGraphica_ID'] = cnv_data['MedGraphica_ID'].str.split(', ')
        cnv_data['Gene_Name'] = cnv_data['Gene_Name'].str.split(', ')

        # Check if both columns have the same number of elements before exploding
        assert cnv_data['MedGraphica_ID'].apply(len).equals(cnv_data['Gene_Name'].apply(len)), "Mismatch in number of elements between 'MedGraphica_ID' and 'Gene_Name'"

        # Explode both columns so that each pair of MedGraphica_ID and gene_name is on its own row
        cnv_data_exploded = cnv_data.explode(['MedGraphica_ID', 'Gene_Name'])

        # Strip whitespace to ensure clean gene names and medgraphica_ids
        cnv_data_exploded['Gene_Name'] = cnv_data_exploded['Gene_Name'].str.strip()
        cnv_data_exploded['MedGraphica_ID'] = cnv_data_exploded['MedGraphica_ID'].str.strip()

        # Filter out rows where 'MedGraphica_ID' or 'Gene_Name' is empty
        cnv_data_filtered = cnv_data_exploded[(cnv_data_exploded['MedGraphica_ID'] != '') & (cnv_data_exploded['Gene_Name'] != '')]
        sample_columns = cnv_data_filtered.columns[8:]
        
        # Group by 'MedGraphica_ID', 'Gene_Name', and 'SVTYPE' and sum the sample columns
        cnv_aggregated = cnv_data_filtered.groupby(['MedGraphica_ID', 'Gene_Name', 'SVTYPE'])[sample_columns].sum().reset_index()
        
        # Get unique values in 'SVTYPE'
        unique_svtypes = cnv_aggregated['SVTYPE'].unique()
        print(f"Unique SVTYPEs: {unique_svtypes}")

        # Iterate over unique SVTYPEs and create separate files for each
        for svtype in unique_svtypes:
            svtype_df = cnv_aggregated[cnv_aggregated['SVTYPE'] == svtype]

            # Remove the 'SVTYPE' column as it's no longer needed
            svtype_df = svtype_df.drop(columns='SVTYPE')

            # Define the output file path
            output_file_path = f'cache/{entity_type.lower()}_{feature_label.lower()}_{svtype.lower()}.csv'
            svtype_df.to_csv(output_file_path, sep=",", index=False)
            print(f"Exported {svtype} data to {output_file_path}")

        print(f"CNV data processing completed.")
    
    def process_gene_methyl(gene_entity_data):
        """Process Gene Methylation data."""
        # Read the methylation data
        methylation = pd.read_csv(file_path, sep=sep)
        
        # Print the shape of the original methylation data
        print("Original methylation data shape:", methylation.shape)

        # Expand gene_entity_data so that each gene name in id_type gets its own row
        gene_entity_expanded = gene_entity_data.copy()
        gene_entity_expanded = gene_entity_expanded.assign(
            **{id_type: gene_entity_expanded[id_type].str.split(';')}
        ).explode(id_type)

        # Merge with methylation data on selected_column
        methylation_merged = pd.merge(
            methylation,
            gene_entity_expanded[[id_type, 'MedGraphica_ID']],
            left_on=selected_column,
            right_on=id_type,
            how='left'
        )

        # Drop the extra merge key (id_type)
        methylation_merged.drop(columns=[id_type], inplace=True)

        # Rename selected_column to Gene_Name
        methylation_merged.rename(columns={selected_column: 'Gene_Name'}, inplace=True)

        # Ensure MedGraphica_ID is the first column
        cols = ['MedGraphica_ID'] + [col for col in methylation_merged.columns if col != 'MedGraphica_ID']
        methylation_merged = methylation_merged[cols]

        # Drop rows where MedGraphica_ID is NaN or empty
        methylation_merged.dropna(subset=['MedGraphica_ID'], inplace=True)
        methylation_merged = methylation_merged[methylation_merged['MedGraphica_ID'] != '']

        # Print the shape of the merged methylation data
        print("Merged methylation data shape:", methylation_merged.shape)

        # Find unique values in the 'Region' column
        unique_regions = methylation_merged['Region'].unique()
        print("Unique Regions:", unique_regions)

        # Iterate over unique regions and create a new DataFrame for each region, then export
        for region in unique_regions:
            # Create a DataFrame specific to the current region
            region_df = methylation_merged[methylation_merged['Region'] == region]
            
            # Convert the region name to lowercase and replace spaces with underscores for the file name
            region_name = region.lower().replace(' ', '_')
            
            # Define the output file path
            output_file_path = f'cache/{entity_type.lower()}_{feature_label.lower()}_{region_name}.csv'
            
            # Export the region-specific DataFrame to a CSV file
            region_df.to_csv(output_file_path, sep=",", index=False)
            print(f"Exported {region} data to {output_file_path}")


    # Check if the feature_label contains "expression" or "cnv"
    if "expression" in feature_label.lower():
        process_gene_expression(gene_entity_data)
    elif "cnv" in feature_label.lower():
        process_gene_cnv(gene_entity_data)
    elif "meth" in feature_label.lower():
        process_gene_methyl(gene_entity_data)
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

        # # Read the gene expression file
        # gene_expression = pd.read_csv(file_path, sep=sep)
        # print("Original gene_expression shape:", gene_expression.shape)

        # # Process gene_id by stripping the version number
        # gene_ids = gene_expression[selected_column].apply(lambda x: x.split('.')[0]).tolist()
        # gene_expression[selected_column] = gene_ids

        # # Retain only the necessary columns and rename them
        # transcript_entity_data = transcript_entity_data[[id_type, 'MedGraphica_ID', 'Gene_Name']].copy()
        # transcript_entity_data.rename(columns={id_type: 'Gene_ID'}, inplace=True)
        # print("transcript_entity_data shape after selection and rename:", transcript_entity_data.shape)

        # # Drop NaN and duplicates from transcript_entity_data
        # transcript_entity_data = transcript_entity_data.dropna(subset=['MedGraphica_ID']).drop_duplicates()
        # print("transcript_entity_data shape after dropna and drop_duplicates:", transcript_entity_data.shape)

        # # Merge the gene expression data with the transcript entity data on 'gene_id' and keep 'MedGraphica_ID' and 'Gene_Name'
        # merged_data = pd.merge(gene_expression, transcript_entity_data, on='Gene_ID', how='inner')
        # print("Merged data shape:", merged_data.shape)

        # # Replace 'Nan' strings with actual NaN values
        # merged_data.replace('Nan', pd.NA, inplace=True)

        # # Drop rows where 'MedGraphica_ID' is NaN (if any)
        # merged_data.dropna(subset=['MedGraphica_ID'], inplace=True)

        # # Ensure that 'MedGraphica_ID' and 'Gene_Name' are the first columns, followed by the rest in their original order
        # cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in merged_data.columns if col not in ['MedGraphica_ID', 'Gene_Name']]
        # gene_expression_new = merged_data[cols]
        # print("gene_expression_new shape after column reordering:", gene_expression_new.shape)

        # # Select numeric columns
        # numeric_cols = gene_expression_new.select_dtypes(include=[np.number]).columns.tolist()

        # # Group by 'MedGraphica_ID' and calculate the mean for all numeric columns
        # gene_expression_new = gene_expression_new.groupby('MedGraphica_ID')[numeric_cols].mean()
        # print("gene_expression_new shape after grouping by 'MedGraphica_ID':", gene_expression_new.shape)

        # # Reset the index to turn 'MedGraphica_ID' back into a column
        # gene_expression_new.reset_index(inplace=True)

        # # Merge back the gene names after grouping by 'MedGraphica_ID'
        # gene_expression_new = pd.merge(gene_expression_new, transcript_entity_data[['MedGraphica_ID', 'Gene_Name']], on='MedGraphica_ID', how='left')
        # print("Final gene_expression_new shape after merging back gene_name:", gene_expression_new.shape)

        # # Ensure 'MedGraphica_ID' is first, 'Gene_Name' is second, and other columns follow the original order
        # final_cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in gene_expression_new.columns if col not in ['MedGraphica_ID', 'Gene_Name']]
        # gene_expression_new = gene_expression_new[final_cols]

        # # Export the result to a CSV file
        # output_file_path = f'cache/{feature_label.lower()}_{entity_type.lower()}.csv'
        # gene_expression_new.to_csv(output_file_path, sep=",", index=False)


        # print(f"Gene expression data processing completed. Output saved to {output_file_path}")

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
    final_cols = ['MedGraphica_ID', 'Protein_ID', 'Gene_Name'] + [col for col in protein_grouped.columns if col not in ['MedGraphica_ID', 'protein_id', 'Gene_Name']]
    protein_grouped = protein_grouped[final_cols]
    print("Grouped protein data:", protein_grouped.shape)

    # Export the result to a CSV file
    output_file_path = f'cache/{entity_type.lower()}_{feature_label.lower()}.csv'
    protein_grouped.to_csv(output_file_path, sep=",", index=False)

    print(f"Protein data processing completed. Output saved to {output_file_path}")

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
