import pandas as pd

def process_entities(entity_type, id_type, file_path, selected_column, feature_label, matcher, embeddings, selection_callback=None):
    """
    Process entities (Phenotype, Disease, Drug, etc.), allow user selection, and perform further processing.
    """
    print(f"Processing {entity_type}: ID Type: {id_type}, File: {file_path}, Column: {selected_column}")

    # Determine the separator based on the file extension
    if file_path.endswith('.txt') or file_path.endswith('.tsv'):
        sep = '\t'
    elif file_path.endswith('.csv'):
        sep = ','
    else:
        raise ValueError("Unsupported file type. Please provide a .txt, .tsv, or .csv file.")

    data_frame = pd.read_csv(file_path, sep=sep)
    all_topk_entities = {}

    # Get topk matches for each entity value in the selected column
    for entity_value in data_frame[selected_column].dropna():
        topk_entities = matcher.get_topk_entities(query=entity_value, k=3, embeddings=embeddings)
        all_topk_entities[entity_value] = topk_entities
    print(f"{entity_type}: {all_topk_entities}")

    # Allow user selection of the topk entities
    if selection_callback:
        updated_topk_entities = selection_callback(all_topk_entities, entity_type, id_type, file_path, selected_column, feature_label, data_frame)
        if updated_topk_entities is not None:
            all_topk_entities = updated_topk_entities
        else:
            print("Warning: selection_callback returned None")

    return all_topk_entities  # Return selected entities for further processing


def process_entities_file(entity_type, id_type, feature_label, data_frame, selected_column, selected_entities):
    """
    Process the selected entities and update the DataFrame.
    Replace the selected_column's entity_value with the selected match term (e.g., hpo_term or drug term) and
    add a MedGraphica_ID column with corresponding med_id for each selected match.
    """
    data_frame_copy = data_frame.copy()

    for entity_value, selected_matches in selected_entities.items():
        med_id, term = selected_matches[0]

        # Replace entity_value with term in the selected column
        data_frame_copy[selected_column] = data_frame_copy[selected_column].replace(entity_value, term)
        print(f"Replaced {entity_value} with {term} in column {selected_column}")

        # Add MedGraphica_ID column if not exists
        if 'MedGraphica_ID' not in data_frame_copy.columns:
            data_frame_copy['MedGraphica_ID'] = None

        # Insert med_id where the term was inserted
        data_frame_copy.loc[data_frame_copy[selected_column] == term, 'MedGraphica_ID'] = med_id
        print(f"Inserted MedGraphica_ID for rows where {selected_column} is {term}")

    # Reorder the columns to make MedGraphica_ID the first column
    cols = ['MedGraphica_ID'] + [col for col in data_frame_copy.columns if col != 'MedGraphica_ID']
    data_frame_copy = data_frame_copy[cols]

    # Save the updated DataFrame to a CSV file
    output_file_path = f'cache/{feature_label.lower()}.csv'
    data_frame_copy.to_csv(output_file_path, index=False)

    print(f"Processed entity data saved to {output_file_path}")
