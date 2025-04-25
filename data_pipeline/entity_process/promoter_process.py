# promoter_process.py

import pandas as pd
import numpy as np
import os
import pymysql
from sqlalchemy import create_engine, types
import sys
import logging

def process_promoter(entity_type, id_type, file_path, selected_column, feature_label, database_path, fill0=False, sample_ids=None):
    print(f"Processing Promoter - fill0={fill0}, Feature Label: {feature_label}")

    # Load BioMedGraphica mapping table
    promoter_csv_path = os.path.join(database_path, "Entity", "Promoter", "BioMedGraphica_Conn_Promoter.csv")
    entity_data = pd.read_csv(promoter_csv_path)
    bmg_ids = entity_data["BioMedGraphica_Conn_ID"].drop_duplicates().tolist()

    # ========== FILL0 ==========
    if fill0:
        if sample_ids is None:
            raise ValueError("sample_ids must be provided when fill0=True")

        data_matrix = pd.DataFrame(0, index=sample_ids, columns=bmg_ids)
        data_matrix.index.name = "Sample_ID"
        data_matrix.reset_index(inplace=True)

        # output_file_path = f'cache/{feature_label.lower()}.csv'
        # data_matrix.to_csv(output_file_path, index=False)

        npy_output_path = f'cache/{feature_label.lower()}.npy'
        np.save(npy_output_path, data_matrix.drop(columns=["Sample_ID"]).values)
        print(f"[fill0] Expression matrix saved to: {npy_output_path}")

        # Save empty mapping
        map_output_file = f'cache/raw_id_mapping/{feature_label.lower()}_id_map.csv'
        mapping_df = pd.DataFrame({
            "BioMedGraphica_Conn_ID": bmg_ids,
            "Original_ID": ["" for _ in bmg_ids]
        })
        mapping_df.to_csv(map_output_file, index=False)
        print(f"[fill0] Mapping saved to: {map_output_file}")
        return

    # ========== NORMAL ==========
    sep = '\t' if file_path.endswith(('.tsv', '.txt')) else ','
    df = pd.read_csv(file_path, sep=sep)

    # Rename sample ID column
    df.rename(columns={selected_column: "Sample_ID"}, inplace=True)

    # Melt to long format
    melted = df.melt(id_vars="Sample_ID", var_name="Original_ID", value_name="value")

    # Identify columns (genes) that were used in input file
    used_ids = set(df.columns) - {"Sample_ID"}

    # Expand mapping: split id_type field on ';' and explode
    mapping_raw = entity_data[[id_type, "BioMedGraphica_Conn_ID"]].dropna()
    mapping_raw[id_type] = mapping_raw[id_type].astype(str).str.strip()
    mapping_expanded = mapping_raw.assign(
        Original_ID=mapping_raw[id_type].str.split(";")
    ).explode("Original_ID")
    mapping_expanded["Original_ID"] = mapping_expanded["Original_ID"].str.strip()

    # Filter mapping to only those IDs actually used in the input file
    mapping_df = mapping_expanded[mapping_expanded["Original_ID"].isin(used_ids)]
    mapping_df = mapping_df[["Original_ID", "BioMedGraphica_Conn_ID"]].drop_duplicates()

    # Merge data with mapping
    merged = pd.merge(melted, mapping_df, on="Original_ID", how="inner")

    # Pivot to wide matrix
    expr = merged.pivot_table(index="Sample_ID", columns="BioMedGraphica_Conn_ID", values="value", fill_value=0)
    expr = expr.reindex(columns=bmg_ids, fill_value=0)  # Ensure full coverage
    expr = expr.copy()
    expr.reset_index(inplace=True)

    # Save expression matrix
    # output_file_path = f'cache/{feature_label.lower()}.csv'
    # expr.to_csv(output_file_path, index=False)
    # print(f"[normal] Promoter data saved to: {output_file_path}")

    npy_output_path = f'cache/{feature_label.lower()}.npy'
    np.save(npy_output_path, expr.drop(columns="Sample_ID").values)
    print(f"[normal] promoter data saved to: {npy_output_path} (as .npy)")

    # Final mapping: group used Original_IDs by BMG ID
    grouped_mapping_df = (
        mapping_df.groupby("BioMedGraphica_Conn_ID")["Original_ID"]
        .apply(lambda x: ";".join(sorted(set(str(i) for i in x if pd.notna(i) and str(i).strip()))))
        .reset_index()
    )

    # Ensure all BMG IDs included (fill empty if not used)
    final_mapping_df = pd.DataFrame({"BioMedGraphica_Conn_ID": bmg_ids}).merge(
        grouped_mapping_df, on="BioMedGraphica_Conn_ID", how="left"
    ).fillna({"Original_ID": ""})

    # Save mapping table
    map_output_file = f'cache/raw_id_mapping/{feature_label.lower()}_id_map.csv'
    final_mapping_df.to_csv(map_output_file, index=False)
    print(f"[normal] Mapping saved to: {map_output_file}")


# def process_promoter_no_T(entity_type, id_type, file_path, selected_column, feature_label, database_path):
#     """Process Gene data."""
#     print(f"Processing Gene: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}, Feature Label: {feature_label}")

#     promoter_csv_path = os.path.join(database_path, "Entity", "Promoter", "biomedgraphica_promoter.csv")
#     promoter_entity_data = pd.read_csv(promoter_csv_path)
        
#     # Determine the separator based on the file extension
#     if file_path.endswith('.txt') or file_path.endswith('.tsv'):
#         sep = '\t'
#     elif file_path.endswith('.csv'):
#         sep = ','
#     else:
#         raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")

#     # Step 1: Read the promoter file (in the format Sample_ID | Gene1 | Gene2 | ...)
#     promoter = pd.read_csv(file_path, sep=sep)

#     # Rename the first column to 'Sample_ID' if it is not already named as such
#     if promoter.columns[0] != id_type:
#         promoter.rename(columns={promoter.columns[0]: id_type}, inplace=True)

#     # Step 3: Merge with gene_entity_data to map to MedGraphica_ID
#     promoter_entity_data = promoter_entity_data[[id_type, 'BioMedGraphica_Conn_ID']].copy()

#     # Merge the promoter data with the gene entity data on id_type
#     promoter_merged = pd.merge(
#         promoter,
#         promoter_entity_data,
#         on=id_type,
#         how='inner'
#     )

#     # Step 4: Extract the id_type and MedGraphica_ID columns as the mapping table
#     # This will create a mapping between the original id_type (e.g., Gene_Name) and MedGraphica_ID
#     mapping_table = promoter_merged[[id_type, 'BioMedGraphica_Conn_ID']].drop_duplicates()
#     mapping_table = mapping_table.rename(columns={id_type: 'Original_ID'})

#     # Save the mapping table to a separate CSV file
#     map_output_file = f'cache/raw_id_mapping/{feature_label.lower()}_id_map.csv'
#     mapping_table.to_csv(map_output_file, sep=",", index=False)
#     print(f"Mapping saved to {map_output_file}")

#     # Step 5: Drop the 'id_type' column after merging
#     promoter_merged.drop(columns=[id_type], inplace=True)

#     # Step 6: Ensure 'BioMedGraphica_Conn_ID' is the first column
#     cols = ['BioMedGraphica_Conn_ID'] + [col for col in promoter_merged.columns if col not in ['BioMedGraphica_Conn_ID']]
#     promoter_data = promoter_merged[cols]

#     # Select numeric columns (which are now samples, i.e., S1, S2, etc.)
#     numeric_cols = promoter_data.select_dtypes(include=[np.number]).columns.tolist()

#     # Group by 'BioMedGraphica_Conn_ID' and calculate the mean for all numeric columns
#     promoter_data.loc[:, numeric_cols] = promoter_data[numeric_cols].fillna(0)
#     promoter_data = promoter_data.groupby('BioMedGraphica_Conn_ID')[numeric_cols].mean()

#     # Reset the index to turn 'BioMedGraphica_Conn_ID' back into a column
#     promoter_data.reset_index(inplace=True)

#     # Step 7: Transpose the data back to the original format (Sample_ID | Gene1 | Gene2 | ...)
#     promoter_data_final = promoter_data.set_index('BioMedGraphica_Conn_ID').T.reset_index()
#     promoter_data_final.rename(columns={'index': 'Sample_ID'}, inplace=True)

#     # Export the final processed promoter data to a CSV file
#     output_file_path = f'cache/{feature_label.lower()}.csv'
#     promoter_data_final.to_csv(output_file_path, sep=",", index=False)

#     print(f"Promoter data processing completed. Output saved to {output_file_path}")


# def process_promoter_with_mysql(entity_type, id_type, file_path, selected_column, feature_label, database_config):
#     """
#     Process promoter data using SQL database for id_type and BioMedGraphica_Conn_ID retrieval.
#     """
#     print(f"Processing Gene: {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}, Feature Label: {feature_label}")

#     # Create a SQLAlchemy engine for MySQL
#     engine = create_engine(
#         f"mysql+pymysql://{database_config['user']}:{database_config['password']}@{database_config['host']}:{database_config['port']}/{database_config['database']}"
#     )

#     # Step 1: Retrieve id_type and BioMedGraphica_Conn_ID from the database
#     query = f"SELECT {id_type}, BioMedGraphica_Conn_ID FROM promoter"
#     promoter_entity_data = pd.read_sql(query, engine)
#     print("Data retrieved from SQL database.")

#     # Determine the separator based on the file extension
#     if file_path.endswith('.txt') or file_path.endswith('.tsv'):
#         sep = '\t'
#     elif file_path.endswith('.csv'):
#         sep = ','
#     else:
#         raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")

#     # Step 2: Read the promoter file (in the format Sample_ID | Gene1 | Gene2 | ...)
#     promoter = pd.read_csv(file_path, sep=sep)

#     # Step 3: Transpose the promoter data to Gene_Name | S1 | S2 | ... format
#     promoter_transposed = promoter.set_index(selected_column).T.reset_index()
#     promoter_transposed.rename(columns={'index': id_type}, inplace=True)

#     # Step 4: Merge with gene_entity_data to map to BioMedGraphica_Conn_ID
#     promoter_entity_data = promoter_entity_data[[id_type, 'BioMedGraphica_Conn_ID']].copy()

#     # Merge the promoter data with the gene entity data on id_type
#     promoter_merged = pd.merge(
#         promoter_transposed,
#         promoter_entity_data,
#         on=id_type,
#         how='inner'
#     )

#     # Step 5: Extract the id_type and BioMedGraphica_Conn_ID columns as the mapping table
#     # This will create a mapping between the original id_type (e.g., Gene_Name) and BioMedGraphica_Conn_ID
#     mapping_table = promoter_merged[[id_type, 'BioMedGraphica_Conn_ID']].drop_duplicates()
#     mapping_table = mapping_table.rename(columns={id_type: 'Original_ID'})

#     # Save the mapping table to a separate CSV file
#     map_output_file = f'cache/raw_id_mapping/{feature_label.lower()}_id_map.csv'
#     mapping_table.to_csv(map_output_file, sep=",", index=False)
#     print(f"Mapping saved to {map_output_file}")

#     # Step 6: Drop the 'id_type' column after merging
#     promoter_merged.drop(columns=[id_type], inplace=True)

#     # Step 7: Ensure 'BioMedGraphica_Conn_ID' is the first column
#     cols = ['BioMedGraphica_Conn_ID'] + [col for col in promoter_merged.columns if col not in ['BioMedGraphica_Conn_ID']]
#     promoter_data = promoter_merged[cols]

#     # Step 8: Select numeric columns (which are now samples, i.e., S1, S2, etc.)
#     numeric_cols = promoter_data.select_dtypes(include=[np.number]).columns.tolist()

#     # Group by 'BioMedGraphica_Conn_ID' and calculate the mean for all numeric columns
#     promoter_data.loc[:, numeric_cols] = promoter_data[numeric_cols].fillna(0)
#     promoter_data = promoter_data.groupby('BioMedGraphica_Conn_ID')[numeric_cols].mean()

#     # Reset the index to turn 'BioMedGraphica_Conn_ID' back into a column
#     promoter_data.reset_index(inplace=True)

#     # Step 9: Transpose the data back to the original format (Sample_ID | Gene1 | Gene2 | ...)
#     promoter_data_final = promoter_data.set_index('BioMedGraphica_Conn_ID').T.reset_index()
#     promoter_data_final.rename(columns={'index': 'Sample_ID'}, inplace=True)

#     # Step 10: Export the final processed promoter data to a CSV file
#     output_file_path = f'cache/{feature_label.lower()}.csv'
#     promoter_data_final.to_csv(output_file_path, sep=",", index=False)

#     print(f"Promoter data processing completed. Output saved to {output_file_path}")

# if __name__ == "__main__":
#     # Define MySQL database connection configuration
#     database_config = {
#         "host": "localhost",           
#         "user": "root",             
#         "password": "Xtq!001028",    
#         "database": "biomedgraphica",  
#         "port": 3306                    
#     }

#     # Define test parameters
#     entity_type = "Promoter"  
#     id_type = "HGNC_Symbol"                
#     selected_column = "Patient_ID"       
#     # file_path = rf"cache\_x_141\methylation_141.csv"
#     file_path = rf"input_data\UCSC_test_data\methylation.csv"
#     feature_label = "Promoter"       

#     # Call the function
#     process_promoter_with_mysql(entity_type, id_type, file_path, selected_column, feature_label, database_config)