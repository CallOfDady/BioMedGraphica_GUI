# transcript_process.py

import pandas as pd
import numpy as np
import os

def process_transcript(entity_type, id_type, file_path, selected_column, feature_label, database_path, fill0=False, sample_ids=None):
    print(f"Processing transcript - fill0={fill0}, Feature Label: {feature_label}")

    # Load BioMedGraphica mapping table
    transcript_csv_path = os.path.join(database_path, "Entity", "Transcript", "BioMedGraphica_Conn_Transcript.csv")
    entity_data = pd.read_csv(transcript_csv_path)
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

    # # Save expression matrix
    # output_file_path = f'cache/{feature_label.lower()}.csv'
    # expr.to_csv(output_file_path, index=False)
    # print(f"[normal] transcript data saved to: {output_file_path}")

    npy_output_path = f'cache/{feature_label.lower()}.npy'
    np.save(npy_output_path, expr.drop(columns="Sample_ID").values)
    print(f"[normal] transcript data saved to: {npy_output_path} (as .npy)")

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
