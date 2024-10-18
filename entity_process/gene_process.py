# gene_process.py

import pandas as pd
import numpy as np

def process_gene(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Gene data."""
    print(f"Processing Gene: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}, Feature Label: {feature_label}")
    
    # Load the gene entity data
    gene_entity_data = pd.read_csv('data/MedGraphica/Node/Gene/medgraphica_gene.csv')    
        
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

    def process_gene_cnv(gene_entity_data):
        """Process Gene CNV data."""
        print(f"Processing Gene CNV for {entity_type} with Feature Label: {feature_label}")
        
        # Read the CNV file
        cnv = pd.read_csv(file_path, sep=sep)
        print(f"CNV data loaded. Columns available: {cnv.columns.tolist()}")

        if id_type == 'Locus-based ID':
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
            
            # Split the selected_column (e.g., "START, END") into two column names
            start_col, end_col = selected_column.split(', ')
            print(f"Selected columns for START: {start_col}, END: {end_col}")
            
            # Find the chromosome column in the CNV DataFrame
            chrom_col = [col for col in cnv.columns if 'CHROM' in col.upper()]
            if not chrom_col:
                raise ValueError("No chromosome column found with 'CHROM' in the name.")
            chrom_col = chrom_col[0]
            print(f"Chromosome column identified: {chrom_col}")
            
            # Initialize columns in the original CNV DataFrame for MedGraphica_ID and Gene_Name
            cnv['MedGraphica_ID'] = None
            cnv['Gene_Name'] = None

            # Iterate over each row in the CNV DataFrame and perform gene name lookup
            for index, row in cnv.iterrows():
                chrom = row[chrom_col]  # Get the chromosome value from the specified chromosome column
                start = row[start_col]  # Get the start value from the specified start column
                end = row[end_col]      # Get the end value from the specified end column
                
                # Call the get_gene_names function to retrieve the MedGraphica_ID and Gene_Name for the current chrom, start, and end range
                try:
                    medGraphica_id, gene_name = get_gene_names(chrom, start, end)
                except Exception as e:
                    print(f"Error during gene name lookup for row {index}: {e}")
                    medGraphica_id, gene_name = None, None
                
                # Assign MedGraphica_ID and Gene_Name directly to the CNV DataFrame
                cnv.at[index, 'MedGraphica_ID'] = medGraphica_id
                cnv.at[index, 'Gene_Name'] = gene_name

            print("Finished processing all rows. MedGraphica_ID and Gene_Name columns have been added to CNV data.")

        
            cnv.dropna(subset=['MedGraphica_ID'], inplace=True)

            # Drop the chromosome, start, and end columns
            cnv.drop([chrom_col, start_col, end_col], axis=1, inplace=True)

            # Split 'MedGraphica_ID' and 'Gene_Name' by comma and ensure they're synchronized
            cnv['MedGraphica_ID'] = cnv['MedGraphica_ID'].str.split(', ')
            cnv['Gene_Name'] = cnv['Gene_Name'].str.split(', ')

            # Check if both columns have the same number of elements before exploding
            assert cnv['MedGraphica_ID'].apply(len).equals(cnv['Gene_Name'].apply(len)), "Mismatch in number of elements between 'MedGraphica_ID' and 'Gene_Name'"

            # Explode both columns so that each pair of MedGraphica_ID and Gene_Name is on its own row
            cnv_exploded = cnv.explode(['MedGraphica_ID', 'Gene_Name'])

            numeric_cols = cnv_exploded.select_dtypes(include=[np.number]).columns.tolist()

            # Group by 'MedGraphica_ID' and calculate the mean for all numeric columns
            cnv_exploded = cnv_exploded.groupby('MedGraphica_ID')[numeric_cols].mean()

            # Reset the index to turn 'MedGraphica_ID' back into a column
            cnv_exploded.reset_index(inplace=True)

            # Merge with gene_entity_data to retrieve the Gene_Name (if needed)
            cnv_exploded = pd.merge(cnv_exploded, gene_entity_data[['MedGraphica_ID', 'HGNC_Symbol']], on='MedGraphica_ID', how='left')

            cnv_exploded.rename(columns={'HGNC_Symbol': 'Gene_Name'}, inplace=True)

            # Remove everything after the first hyphen in Gene_Name (optional, if necessary)
            cnv_exploded['Gene_Name'] = cnv_exploded['Gene_Name'].apply(lambda x: x.split('-')[0] if pd.notnull(x) else x)

            # Print the shape of the exploded CNV data
            print(f"Exploded CNV data shape: {cnv_exploded.shape}")

            # Reorder columns
            final_cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in cnv_exploded.columns if col not in ['MedGraphica_ID', 'Gene_Name']]
            cnv_exploded = cnv_exploded[final_cols]

            # Save the processed CNV data to a CSV file
            output_file_path = f'cache/{feature_label.lower()}.csv'
            cnv_exploded.to_csv(output_file_path, index=False)
            print(f"CNV data processing completed.")
        else:
            print(f"ID Type: {id_type} is not supported for CNV data processing.")
    
    # Check if the feature_label contains "expression" or "cnv"
    if "expression" in feature_label.lower():
        process_gene_expression(gene_entity_data)
    elif "cnv" in feature_label.lower():
        process_gene_cnv(gene_entity_data)
    else:
        print(f"Processing for Feature Label: {feature_label} is not implemented yet.")
        # Leave space for other feature_label handling in the future