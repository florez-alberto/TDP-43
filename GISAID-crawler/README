# GISAID Data Downloader Scripts

This repository contains two Python scripts that automate the process of downloading data from the GISAID website.

Please replace "your_username" and "your_password" with your appropriate credentials granted by GISAID.

## Scripts Description

1. **gisaid_downloaderv2stable.py**: This script logs into the GISAID website, navigates to a specific page, and downloads data based on specific dates. The dates are read from a CSV file named 'datesDF.csv'. The script uses a while loop to process each date in the CSV file. For each date, it checks if there are any files to download. If there are, it downloads the files, renames them, and updates the 'processed' status in the CSV file. Note: the last column should be set to 'False' in order for the script to process it.

2. **gisaid_meta_downloaderv2stable.py**: This script performs the same tasks as the previous script, but it also downloads all associated metadata, including the acknowledgements and host by isolateId.

Both scripts include error handling to deal with situations where certain elements on the webpage might not be available. They also include several `time.sleep()` calls to pause execution of the script, allowing for delays in page loading or file downloading. The files will be programmatically downloaded into the specified download folder.


