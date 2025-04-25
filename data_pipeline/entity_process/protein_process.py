# protein_process.py

import pandas as pd
import numpy as np
import os

# def process_protein(entity_type, id_type, file_path, selected_column, feature_label, database_path, fill0=False):
#     """Process Protein data, map Sample_ID to MedGraphica_ID, and export mapping."""
#     print(f"Processing Protein for {entity_type}, ID Type: {id_type}, File: {file_path}, Column: {selected_column}, Feature Label: {feature_label}")

#     protein_csv_path = os.path.join(database_path, "Entity", "Protein", "BioMedGraphica_Conn_Protein.csv")
#     protein_entity_data = pd.read_csv(protein_csv_path)
        
#     # Determine the separator based on the file extension
#     if file_path.endswith('.txt') or file_path.endswith('.tsv'):
#         sep = '\t'
#     elif file_path.endswith('.csv'):
#         sep = ','
#     else:
#         raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")
    
#     protein = pd.read_csv(file_path, sep=sep)

#     # Step 2: Transpose the protein data to Protein_Name | S1 | S2 | ... format
#     protein_transposed = protein.set_index(selected_column).T.reset_index()
#     protein_transposed.rename(columns={'index': id_type}, inplace=True)

#     # Step 3: Merge with gene_entity_data to map to MedGraphica_ID
#     protein_entity_data = protein_entity_data[[id_type, 'BioMedGraphica_Conn_ID']].copy()

#     # Merge the protein data with the gene entity data on id_type
#     protein_merged = pd.merge(
#         protein_transposed,
#         protein_entity_data,
#         on=id_type,
#         how='inner'
#     )

#     # Step 4: Extract the id_type and MedGraphica_ID columns as the mapping table
#     # This will create a mapping between the original id_type (e.g., Protein_Name) and MedGraphica_ID
#     mapping_table = protein_merged[[id_type, 'BioMedGraphica_Conn_ID']].drop_duplicates()
#     mapping_table = mapping_table.rename(columns={id_type: 'Original_ID'})

#     # Save the mapping table to a separate CSV file
#     map_output_file = f'cache/raw_id_mapping/{feature_label.lower()}_id_map.csv'
#     mapping_table.to_csv(map_output_file, sep=",", index=False)
#     print(f"Mapping saved to {map_output_file}")

#     # Step 5: Drop the 'id_type' column after merging
#     protein_merged.drop(columns=[id_type], inplace=True)

#     # Step 6: Ensure 'BioMedGraphica_Conn_ID' is the first column
#     cols = ['BioMedGraphica_Conn_ID'] + [col for col in protein_merged.columns if col not in ['BioMedGraphica_Conn_ID']]
#     protein_data = protein_merged[cols]

#     # Select numeric columns (which are now samples, i.e., S1, S2, etc.)
#     numeric_cols = protein_data.select_dtypes(include=[np.number]).columns.tolist()

#     # Group by 'BioMedGraphica_Conn_ID' and calculate the mean for all numeric columns
#     protein_data.loc[:, numeric_cols] = protein_data[numeric_cols].fillna(0)
#     protein_data = protein_data.groupby('BioMedGraphica_Conn_ID')[numeric_cols].mean()

#     # Reset the index to turn 'BioMedGraphica_Conn_ID' back into a column
#     protein_data.reset_index(inplace=True)

#     # Step 7: Transpose the data back to the original format (Sample_ID | Protein1 | Protein2 | ...)
#     protein_data_final = protein_data.set_index('BioMedGraphica_Conn_ID').T.reset_index()
#     protein_data_final.rename(columns={'index': 'Sample_ID'}, inplace=True)

#     # Export the final processed protein data to a CSV file (without additional suffix)
#     output_file_path = f'cache/{feature_label.lower()}.csv'
#     protein_data_final.to_csv(output_file_path, sep=",", index=False)

#     print(f"Protein data processing completed. Output saved to {output_file_path}")

def process_protein(entity_type, id_type, file_path, selected_column, feature_label, database_path, fill0=False, sample_ids=None):
    print(f"Processing protein - fill0={fill0}, Feature Label: {feature_label}")

    # Load BioMedGraphica mapping table
    protein_csv_path = os.path.join(database_path, "Entity", "Protein", "BioMedGraphica_Conn_Protein.csv")
    entity_data = pd.read_csv(protein_csv_path)
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
    # print(f"[normal] protein data saved to: {output_file_path}")

    npy_output_path = f'cache/{feature_label.lower()}.npy'
    np.save(npy_output_path, expr.drop(columns="Sample_ID").values)
    print(f"[normal] protein data saved to: {npy_output_path} (as .npy)")

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