Merging 2 Datasets based on protein changes 

This script was designed to merge two datasets from two sequencing assays: a Foundation Medicine assay and TruSight Oncology 500. The Python script synchronizes the protein change format between the two datasets and parses the protein change information to merge the datasets on. Other information such as amplifications, tumor mutation burden, and microsatellite status are also parsed from the datasets.

Pandas is used for the data clean-up and merging and Openpyxl is used to write the output to an Excel file and to create a more visually appealing appearance.

Due to the irregularity of the datasets, 2 functions are used in order to manually clean up some of the data that could not be regularly accounted for and identify rows where specific data stopped and started.


 