Assessment of modelling strategies for drug response prediction in cell lines and xenografts
=================================================================


Abstract
--------

Data from several large high-throughput drug response screens have become available to the scientific community recently. Although many efforts have been made to use this information to predict drug sensitivity, our ability to accurately predict drug response based on genetic data remains limited. In order to systematically examine how different aspects of modelling affect the resulting prediction accuracy, we built a range of models for seven drugs (erlotinib, pacliatxel, lapatinib, PLX4720, sorafenib, nutlin-3 and nilotinib) using data from the largest available cell line and xenograft drug sensitivity screens. We found that the drug response metric, the choice of the molecular data type and the number of training samples have a substantial impact on prediction accuracy. We also compared the tasks of drug response prediction with tissue type prediction and found that, unlike for drug response, tissue type can be predicted with high accuracy. Furthermore, we assessed our ability to predict drug response in four xenograft cohorts (treated either with erlotinib, gemcitabine or paclitaxel) using models trained on cell line data. We could predict response in an erlotinib-treated cohort with a moderate accuracy (correlation≈0.5), but were unable to correctly predict responses in cohorts treated with gemcitabine or paclitaxel.

Citation
--------

To cite this work in publication, use

Kurilov R., Haibe-Kains B., Brors B. 2020. Assessment of modelling strategies for drug response prediction in cell lines and xenografts. _Scientific Reports_


Reproducibility of the Analysis Results
--------------------------------------------
In all scripts, the user needs to change the variable **path** to the repository root folder.

Analysis I
-------------------------------

Scripts should be executed in the following order:

|     | name | description                              |
|-----|------|------------------------------------------|
|1. | data_for_analysis1.R  | obtains data from PharmacoGx and preprocess it for a subsequent analysis |
|2. | analysis.R            | performs the analysis and saves results into results_table.RData |
|3. | plotting_R2.R         | produces figures 2a, 2b and supplementary figures s3, s4, s6, s7, s8, s9 |
|4. | top_features.R        | calculates and saves molecular features associated with drug response vectors and produces figure 2c |
|5. | modelling_obs_pred.R  | produces predictions for each modelling task and saves them into raw_predictions2.RData |
|6. | plotting_obs_vs_pred.R| produces figure 3 and supplementary figures s5-1, s5-2, s5-3 |


Analysis II
-------------------------------
Data: gCSI molecular data, particularly files “gcsi.genomics.rda”, “gcsi.genomics.feature.info.rda”, 
and “gcsi.line.info.rda” should be downloaded from [here](http://research-pub.gene.com/gCSI-cellline-data/compareDrugScreens_current.tar.gz) and extracted in the folder analysis2/data

Scripts should be executed in the following order:

|     | name | description                                   |
|-----|------|-----------------------------------------------|
|1. |nibr_preprocessing.R                  | preprocess molecular and drug response xenograft data |
|2. |main_analysis.R                       | performs the analysis, saves results and produces supplementary figures s11-1 and s11-2 |
|3. |tissue_classification_equal_groups.R  | performs tissue classification for the case with equal number of samples per each tissue and produces figure 4 |
|4. |plotting_obs_vs_pred.R                | produces figure 5a |
|5. |drug_resp_gcsi_to_nibr.R              | creates and tests model for xenograft predictions trained on cell line data, produces figure 5b and supplementary figure s10 |

SessionInfo
------------------------

Analysis has been done in the following session environment. 

#sessionInfo()

R version 3.4.4 (2018-03-15)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] PharmacoGx_1.8.3 ggfortify_0.4.5  gam_1.15         foreach_1.4.4    caret_6.0-78    
 [6] ggplot2_2.2.1    lattice_0.20-35  survcomp_1.28.5  prodlim_1.6.1    survival_2.42-3 

loaded via a namespace (and not attached):
  [1] nlme_3.1-137        lsa_0.73.1          survivalROC_1.0.3   bitops_1.0-6       
  [5] lubridate_1.7.3     dimRed_0.1.0        RColorBrewer_1.1-2  SnowballC_0.5.1    
  [9] tools_3.4.4         R6_2.2.2            rpart_4.1-13        KernSmooth_2.23-15 
 [13] sm_2.2-5.4          lazyeval_0.2.1      BiocGenerics_0.24.0 colorspace_1.3-2   
 [17] rmeta_3.0           nnet_7.3-12         withr_2.1.1         tidyselect_0.2.4   
 [21] gridExtra_2.3       mnormt_1.5-5        compiler_3.4.4      Biobase_2.38.0     
 [25] slam_0.1-42         caTools_1.17.1      scales_0.5.0        sfsmisc_1.1-2      
 [29] DEoptimR_1.0-8      psych_1.7.8         robustbase_0.92-8   randomForest_4.6-12
 [33] relations_0.6-7     stringr_1.3.0       digest_0.6.15       foreign_0.8-70     
 [37] pkgconfig_2.0.1     plotrix_3.7         limma_3.34.9        maps_3.3.0         
 [41] rlang_0.2.0         ddalpha_1.3.1.1     MLmetrics_1.1.1     SuppDists_1.1-9.4  
 [45] bindr_0.1           BiocParallel_1.12.0 gtools_3.5.0        dplyr_0.7.4        
 [49] ModelMetrics_1.1.0  magrittr_1.5        Matrix_1.2-14       Rcpp_0.12.15       
 [53] celestial_1.4.1     munsell_0.4.3       piano_1.18.1        stringi_1.1.6      
 [57] MASS_7.3-50         gplots_3.0.1        plyr_1.8.4          recipes_0.1.2      
 [61] grid_3.4.4          gdata_2.18.0        parallel_3.4.4      mapproj_1.2.6      
 [65] pillar_1.2.1        fgsea_1.4.1         igraph_1.2.1        xgboost_0.6.4.1    
 [69] marray_1.56.0       reshape2_1.4.3      codetools_0.2-15    stats4_3.4.4       
 [73] fastmatch_1.1-0     CVST_0.2-1          NISTunits_1.0.1     glue_1.2.0         
 [77] downloader_0.4      data.table_1.10.4-3 bootstrap_2017.2    gtable_0.2.0       
 [81] RANN_2.5.1          purrr_0.2.4         tidyr_0.8.0         kernlab_0.9-25     
 [85] assertthat_0.2.0    DRR_0.0.3           gower_0.1.2         broom_0.4.3        
 [89] pracma_2.1.4        e1071_1.6-8         class_7.3-14        timeDate_3043.102  
 [93] RcppRoll_0.2.2      tibble_1.4.2        iterators_1.0.9     cluster_2.0.7-1    
 [97] sets_1.0-18         bindrcpp_0.2        lava_1.6            ROCR_1.0-7         
[101] magicaxis_2.0.3     ipred_0.9-6  

