Rscript /path/to/get_cellfracs_seuratv3.R brca_scRNA_GEP.txt brca_scRNA_celllabels.txt brca_STdata_slide1_GEP.txt brca_cell_fraction_estimates.txt

cytospace -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp brca_STdata_slide1_GEP.txt -cp brca_STdata_slide1_coordinates.txt -ctfep brca_cell_fraction_estimates.txt
