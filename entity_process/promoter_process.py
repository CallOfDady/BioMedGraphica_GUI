# promoter_process.py

import pandas as pd
import numpy as np

def process_promoter(entity_type, id_type, file_path, selected_column, feature_label):
    """Process Gene data."""
    print(f"Processing Gene: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}, Feature Label: {feature_label}")
    
    # Load the gene entity data
    promoter_entity_data = pd.read_csv('data/MedGraphica/Node/Gene/medgraphica_promoter.csv')    
        
    # Determine the separator based on the file extension
    if file_path.endswith('.txt') or file_path.endswith('.tsv'):
        sep = '\t'
    elif file_path.endswith('.csv'):
        sep = ','
    else:
        raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")

    def process_promoter_methyl(promoter_entity_data):
        """Process Gene Methylation data."""
        # Read the methylation data
        methylation = pd.read_csv(file_path, sep=sep)
        
        # Print the shape of the original methylation data
        print("Original methylation data shape:", methylation.shape)

        if id_type == 'HGNC_Symbol':
            promoter_entity_data = promoter_entity_data[['MedGraphica_ID', 'HGNC_Symbol']].copy()
        else:
            promoter_entity_data = promoter_entity_data[[id_type, 'MedGraphica_ID', 'HGNC_Symbol']].copy()

        # Expand promoter_entity_data so that each gene name in id_type gets its own row
        promoter_entity_data_expanded = promoter_entity_data.copy()
        promoter_entity_data_expanded = promoter_entity_data_expanded.assign(
            **{id_type: promoter_entity_data_expanded[id_type].str.split(';')}
        ).explode(id_type)

        # Merge with methylation data on selected_column
        methylation_merged = pd.merge(
            methylation,
            promoter_entity_data_expanded,
            left_on=selected_column,
            right_on=id_type,
            how='inner'
        )

        methylation_merged.dropna(subset=['MedGraphica_ID'], inplace=True)

        promoter_entity_data.rename(columns={'HGNC_Symbol': 'Gene_Name'}, inplace=True)
        methylation_merged.rename(columns={'HGNC_Symbol': 'Gene_Name'}, inplace=True)

        # Ensure MedGraphica_ID is the first column
        cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in methylation_merged.columns if col not in ['MedGraphica_ID', 'Gene_Name', selected_column, id_type]]
        methylation_data = methylation_merged[cols]

        numeric_cols = methylation_data.select_dtypes(include=[np.number]).columns.tolist()

        # Group by 'MedGraphica_ID' and calculate the mean for all numeric columns
        methylation_data = methylation_data.groupby('MedGraphica_ID')[numeric_cols].mean()

        # Reset the index to turn 'MedGraphica_ID' back into a column
        methylation_data.reset_index(inplace=True)

        # Merge back the gene names after grouping by 'MedGraphica_ID'
        methylation_data = pd.merge(methylation_data, promoter_entity_data[['MedGraphica_ID', 'Gene_Name']], on='MedGraphica_ID', how='left')

        # Remove (-) after the first hyphen
        methylation_data['Gene_Name'] = methylation_data['Gene_Name'].apply(lambda x: x.split('-')[0] if pd.notnull(x) else x)

        # Print the shape of the merged methylation data
        print("Merged methylation data shape:", methylation_merged.shape)
        
        # Ensure 'MedGraphica_ID' is first, 'Gene_Name' is second, and other columns follow the original order
        final_cols = ['MedGraphica_ID', 'Gene_Name'] + [col for col in methylation_data.columns if col not in ['MedGraphica_ID', 'Gene_Name']]
        methylation_data = methylation_data[final_cols]

        output_file_path = f'cache/{feature_label.lower()}.csv'
        methylation_data.to_csv(output_file_path, index=False)


    
    # Check if the feature_label contains "expression" or "cnv"
    if "meth" in feature_label.lower():
        process_promoter_methyl(promoter_entity_data)
    elif "xxx" in feature_label.lower():
        print(f"Processing for Feature Label: {feature_label} is not implemented yet.")
    else:
        print(f"Processing for Feature Label: {feature_label} is not implemented yet.")
        # Leave space for other feature_label handling in the future