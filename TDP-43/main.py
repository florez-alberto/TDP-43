import main_utils
import gisaid_database_merger

# ? STEP 1 download the files from gisaid (will create several csv files)
#! run gisaid_downloaderv2stable.py
#! run gisaid_meta_downloaderv2stable.py


# ? STEP2 merge the downloaded files in the download folder to single files
gisaid_database_merger.merge_downloads()
gisaid_database_merger.merge_metadata_to_a_single_csv()


# ? STEP 3 create and populate the database 2
main_utils.create_gisaid_DB2()
main_utils.write_metadata_to_db2()
main_utils.write_merged_sequences_to_db2()

# ? STEP 4 create and populate the database 3
# * we will focus on matching the motifs for TDP43 on the DNA sequences
# * to have a clear database we will create a new db3 and pass only the sequences that correspond to clean segments based on the cannonical max lengths

cannonical_max_lengths= [ { "segment": "HA", "minSize": 1695  , "maxSize":   1777 }, { "segment": "NA", "minSize":  1359 , "maxSize": 1467}, { "segment": "MP", "minSize":979 , "maxSize":  1027 }, { "segment": "NP", "minSize": 1494 , "maxSize": 1565 }, { "segment": "NS", "minSize": 835, "maxSize": 890 }, { "segment": "PA", "minSize":  2148, "maxSize": 2233}, { "segment": "PB1", "minSize":  2271, "maxSize": 2341 }, { "segment": "PB2", "minSize":  2271, "maxSize": 2341 } ]

# # *populate the database3 from db2 data with the cannonical_max_lengths and without using illegal characters
main_utils.create_gisaid_DB3()
main_utils.populate_db3_with_sequences(cannonical_max_lengths)

# # # * match sequences that contain all the proteins form same isolateId
main_utils.add_matching_isolateId_to_db3()
print("db3 created and populated")

# ? STEP 5 search the sequences by the motifs "GTGTGT" and "TGTGTG"
# only the first two YG motifs
search_strings= ["GTGTGT","TGTGTG"]

# # * only on complete genomes
main_utils.search_sequence_by_multiple_strings_db3_onlyMatched_genesymbols(search_strings,"all" )

# # * on all clean sequences in the database
print("on all genesymbols")
main_utils.search_sequence_by_multiple_strings_db3(search_strings,"all" )

#find how many sequences do not contain at least one of the motifs any of the genesymbols
main_utils.search_isolateId_with_at_least_one_match()


# ? STEP 6 search the sequences by the motifs all YG motifs
# all the YG motifs
search_strings= ["GTGTGT","TGTGTG", "TGTGCG", "TGCGTG", "CGTGTG", "GTGTGC" ]

# # * only on complete genomes
main_utils.search_sequence_by_multiple_strings_db3_onlyMatched_genesymbols(search_strings,"all" )

# # * on all clean sequences in the database
print("on all genesymbols")
main_utils.search_sequence_by_multiple_strings_db3(search_strings,"all" )

#find how many sequences do not contain at least one of the motifs any of the genesymbols
main_utils.search_isolateId_with_at_least_one_match()

# ? STEP 7: search for exact occurrence of motifs in WSN
fastas_directory="cannonical_models/"
main_utils.search_cannonical_model_by_string(fastas_directory,search_strings)

