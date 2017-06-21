# Joint Modeling and Landmarking

This repository contains the code for the paper entitled: *Dynamic Predictions with 
Time-Dependent Covariates in Survival Analysis using Joint Modeling and Landmarking* by
Dimitris Rizopoulos, Geert Molenberghs and Emmanuel M.E.H. Lesaffre.

The author mainly responsible for writing the code and whom readers should approach with
questions or bug reports is [D. Rizopoulos](mailto:d.rizopoulos@erasmusmc.nl).

The folder `case_study` has the following contents:

- `simulated_AoValve.RData` an R workspace that contains a simulated data set with the same
design as the Aortic Valve dataset presented in the paper.

- `AoValve_analysis.R` an R script file that 
    + loads the simulated Aortic Valve dataset, and
    + does all analysis presented in Section 4 of the paper and Section 1 of the 
    supplementary material in this simulated dataset, producing the corresponding tables
    and figures (namely, Figures 2-3 in the main paper, and Table 1-3 in the 
    supplementary material).

- `AoValve_CV.R` an R script file that produces the cross-validated AUC and PE 
measures in the simulated data set (Table 1 of main paper).

The folder `Simulation` has the following contents:

- The folder `intermediate_results` that contains the corresponding R workspaces with the
simulation results

- The R script files `SimulateI.R`, `SimulateII.R`, and `SimulateIII.R` that were used to
simulate data under the corresponding scenarios described in Section 5.1.

- The R script file `simulation_Funs.R` that contains the supporting functions for 
calculating the AUC and PE measures utilizing the true censoring weights.

- The R script file `plots_SuppMaterial.R` that produces the figures with the simulation 
results shown in the supplementary material.

The results have been produced in [R](https://cran.r-project.org/) under the following 
specification:
```r
R version 3.4.0 (2017-04-21)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] splines   parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] xtable_1.8-2         lattice_0.20-35      JMbayes_0.8-61       rstan_2.15.1         StanHeaders_2.15.0-1
 [6] ggplot2_2.2.1        doParallel_1.0.10    iterators_1.0.8      foreach_1.4.3        survival_2.41-3     
[11] nlme_3.1-131        

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.11        RColorBrewer_1.1-2  compiler_3.4.0      plyr_1.8.4          base64enc_0.1-3    
 [6] tools_3.4.0         rpart_4.1-11        digest_0.6.12       checkmate_1.8.2     htmlTable_1.9      
[11] evaluate_0.10       tibble_1.3.3        gtable_0.2.0        rlang_0.1.1         Matrix_1.2-10      
[16] jagsUI_1.4.4        coda_0.19-1         gridExtra_2.2.1     cluster_2.0.6       stringr_1.2.0      
[21] knitr_1.16          rjags_4-6           htmlwidgets_0.8     nnet_7.3-12         stats4_3.4.0       
[26] rprojroot_1.2       grid_3.4.0          data.table_1.10.4   inline_0.3.14       foreign_0.8-68     
[31] rmarkdown_1.6       latticeExtra_0.6-28 Formula_1.2-1       magrittr_1.5        backports_1.1.0    
[36] scales_0.4.1        Hmisc_4.0-3         codetools_0.2-15    htmltools_0.3.6     MASS_7.3-47        
[41] rsconnect_0.8       colorspace_1.3-2    stringi_1.1.5       acepack_1.4.1       lazyeval_0.2.0     
[46] munsell_0.4.3    
```