# Statistical_Analysis_of_MSdata_obtained_through_AntagonismAssay

Statistical Analysis of MS data obtained by single bacterial cultures and Interaction Assay

This script contains the statistical analysis of LC-MS/MS assessed data. Some parts of the R-script were copied from the “A hitchhiker's guide to statistical analysis of Feature-based Molecular Networks” (https://github.com/Functional-Metabolomics-Lab/FBMN-STATS) and further developed. 

The analysis requires a feature-table (csv-format) obtained through mzMine, annotation file (tsv-format) obtained through GNPS (https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp) and a metadata table (txt-format). It implies the blank removal from the bacterial single cultures and from the interaction assays, conducted on solid agar plates. Statistical analysis includes the assessment of the most significant features from the interaction in comparison with metabolites found in the single cultures of the bacterial strains, using Student’s t-test. Graphical illustration was carried out using Vulcano plots.
