# SARS-CoV-2-ensemble-model
Found here is the R code used to generate the mathematical modelling results found in the paper titled "Ensemble modeling of SARS-CoV-2 immune dynamics in immunologically naïve rhesus macaques predicts that potent, early innate immune responses drive viral elimination". All code was written by Catherine Byrne. All queries can be directed to cbyrne@fredhutch.org.

If wanting to reproduce modelling results and figures found in the article, all data and R files should first be downloaded. 

R files can be run in the following order:
1. fitting_sarscov2_ensemble.r
2. analyze_sarscov2_ensemble.r
3. reinfection_sarscov2_ensemble.r

Data comes from _Nelson et al. Mild SARS-CoV-2 infection in rhesus macaques is associated with viral control prior to antigen-specific T cell responses in tissues. Sci Immunol (2022) 7:eabo0535. doi: 10.1126/sciimmunol.abo0535_. The subset of data used in this analysis can be found in the file titled "all_dat_extrapolated.csv". Other data from this study can be found on NCBI GEO under accession GSE196980 and in the Supplementary Materials of the original paper.  

