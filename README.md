# Fall 2021 Bioinformatics Algorithms Group Project

## Team Members
Kian Banihashemi 
Augie Doebling
Kaanan Kharwa
Julian Menendez

## The Plastics
By analyzing the genomic data of several enzymes capable of breaking down PHB plastic, we hope to find similarities suggesting which areas are the most important for this task. This knowledge would allow future researchers to insert this active site into other bacteria in hopes of creating new organisms capable of breaking down PHB.


## Repository Organization
### BioInformatics_Group_Project_The_Plastics.ipynb
Main development notebook. Includes scripts for parsing and processing data as well as intrpeting results from `heatmap_generation.py`.

### heatmap_generation.py
Python program to generate heatmap data from a fasta file. 
Takes in a path to a fasta file as an argument and saves both the alignment cache and the heatmap data to pickle files with timestamps.
usage : `python3 heatmap_generation.py input_data/uniprot_keyword_search.fasta`

### data_caches **dir**
Saved data caches from intermediate runs of heatmap generation code.

### input_data **dir**
Input data files for use by heatmap_generation

### pickle files
Generated heatmap output pickle files for final dataset runs.

