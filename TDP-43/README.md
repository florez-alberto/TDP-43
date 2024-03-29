# GISAID Data Processing Workflow
DOI: 10.5281/zenodo.10652456

This Python script treats and scans the data downloaded from GISAID website. It uses functions from the `main_utils.py` and `gisaid_database_merger.py` libraries to perform such steps.

## Workflow

1. **Step 1**: Download the files from GISAID. This will create several CSV files. This is done by running `gisaid_downloaderv2stable.py` and `gisaid_meta_downloaderv2stable.py`. These are found in the directory `/GISAID-crawler`

2. **Step 2**: Merge every set into a single file using `gisaid_database_merger.merge_downloads()` and `gisaid_database_merger.merge_metadata_to_a_single_csv()`.

3. **Step 3**: Create and populate the database `db2` using `main_utils.create_gisaid_DB2()`, `main_utils.write_metadata_to_db2()`, and `main_utils.write_merged_sequences_to_db2()`.

4. **Step 4**: Create and populate a new database. This step focuses on preparing the sequences to map YG TDP43 binding motifs on the DNA sequences. To have a clear database, a new db file is created and only the sequences that correspond to clean segments based on the canonical max and min lengths are passed. Sequences that contain all the 8 segments from the same isolateId are also added into an object for further analysis.

5. **Step 5**: Search the sequences by the motifs "GTGTGT" and "TGTGTG". This is done on complete genomes and on all clean sequences in the database. The function `main_utils.search_isolateId_with_at_least_one_match()` is used to find how many sequences do not contain at least one of the motifs in any of the gene symbols.

6. **Step 6**: Search the sequences by all YG motifs. This is similar to Step 5.

7. **Step 7**: Search for the exact occurrence of motifs in WSN using `main_utils.search_cannonical_model_by_string(fastas_directory,search_strings)`.

The console output derived from each step will print the information that was used to make the images.