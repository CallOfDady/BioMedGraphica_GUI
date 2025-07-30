# BioMedGraphica_GUI

## ⚠️ This project is no longer maintained

Please check the latest version here: [BioMedGraphica_WebApp](https://github.com/CallOfDady/BioMedGraphica_WebApp)

## Usage Notes

### Database Preparation and Access
Due to access control restrictions on the collected database, we provide users with data download links for convenient access. Additionally, we offer a Jupyter notebook to help with the merging of entities and the harmonization of relations into designated folders.

Once the data has been curated and placed in the appropriate folder, users can run the software locally on their machine in a client-based environment to begin processing.

### Required File Preparation
Before running the software, prepare the necessary files, including:
- **Entity files**
- **Clinical data file**

Each file should follow these guidelines:
- Use a standardized naming convention for sample IDs.
- Follow a consistent database identifier format for features, such as the Ensembl stable Gene ID.

#### Example Format
The first row of each data file should be structured as:
| Sample_ID | Ensembl_Gene_ID_1 | Ensembl_Gene_ID_2 | ... |
|-----------|-------------------|-------------------|-----|
| R0001     | 0.1234            | 0.1345            | ... |
| R0002     | 0.1455            | 0.3123            | ... |
| ...       | ...               | ...               | ... |

Files will be intersected in the final step to obtain a common sample set.

### Starting the GUI
Once files are prepared, launch the GUI. In the **Welcome** tab:
1. Locate the **BioMedGraphica database path** and verify its integrity for smooth processing.
2. Move to the **Import** tab, where each entity file should be provided with a unique label for identification during subsequent processing.
3. Select the appropriate entity type and ID type for each file from the dropdown menu.
4. Use the **file path** button to select each file’s path.
5. When all inputs are complete, click the **Export** button at the top of the page to save inputs as `config.csv`, making future processing easier.
6. Click **Next** to proceed to the **Read** tab.

### File Reading and Validation
In the **Read** tab, the software will read column names to identify the sample ID column in each file. For simplicity and readability:
- Place the sample ID column as the first column in each file and label it as `id`.
- The software will also perform an intersection of sample IDs across all entity files to obtain a common set of samples, reducing storage requirements.

After confirming that the data is read correctly, click **Next** to proceed to the **Process** tab.

### Processing and Finalizing Files
In the **Process** tab:
1. Each entity file will be processed individually and saved to the `./cache` directory in the root directory.
2. Sort entity files by clicking and arranging them from top to bottom in the dialog box.
3. Import clinical data, and aggregate all individually processed entity files into the required format for GNN training with the **Finalize** function.

**Note:** If files contain many columns, the Finalize operation may require significant RAM (20GB+). In case of issues or unexpected exits during individual file processing, reprocess the specific entity file, then use **Finalize** to aggregate it with the others. If starting a new process, click the **Clear Cache Folder** button in the top-right corner to clear the cache.

### Exporting Processed Files
Once processing is complete, click **Next** to enter the **Export** tab:
1. Preview all processed files in the cache folder.
2. After confirming accuracy, set the save path and click **Export**.

When all tasks are finished, click **Exit** to close the software.
