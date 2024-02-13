# TDP-43 Motifs in the GISAID database
<a href="https://zenodo.org/doi/10.5281/zenodo.10652455"><img src="https://zenodo.org/badge/756591816.svg" alt="DOI"></a>
  
This repository contains Python scripts that automate the process of downloading, merging, and processing data from the GISAID website. The scripts are organized into two directories: `GISAID-crawler` and `TDP-43`, each with its own README file detailing the specific operations performed by the scripts within.

To run each script, please set your working directory to the following:

## GISAID-crawler

To automate the process of downloading data and associated metadata from the GISAID website.

1. **gisaid_downloaderv2stable.py**: This script logs into the GISAID website, navigates to a specific page, and downloads data based on specific dates. The dates are read from a CSV file named 'datesDF.csv'. 

2. **gisaid_meta_downloaderv2stable.py**: This script performs the same tasks as the previous script, but it also downloads all associated metadata, including the acknowledgements and host by isolateId.

For more details, refer to the README file in the `GISAID-crawler` directory.

**Disclaimer**: in order to access the information in the GISAID database you must have your own access by creating a username and being given a password. In these scripts, I did not include any data contained in the database, in accordance with the GISAID terms and conditions. This data has not been shared to anyone nor cross-examined with any other influenza database. A separate table acknowledging all sources of the original data will be added.

## TDP-43

This directory contains scripts that handle the merging and processing of the downloaded data.

1. **main.py**: This is the main script that contains the workflow for processing the data. It uses functions from the `main_utils.py` and `gisaid_database_merger.py` libraries to perform such steps.

It handles creating and populating databases, preparing sequences to map YG TDP43 binding motifs on the DNA sequences, and searching the sequences by specific motifs.

For more details, refer to the README file in the `TDP-43` directory.



License: MIT

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
