# Fast-latent-correlation


R code with functions that were used for simulations and data analysis in the paper ["Fast computation of latent correlations" by Grace Yoon, Christian L. M&uuml;ller and Irina Gaynanova](https://arxiv.org/abs/2006.13875).

### Prerequisites

The functions for ORG and MLBD methods are available from R package **mixedCCA**

To install the most recent development version:
```{r}
devtools::install_github("irinagain/mixedCCA", ref = "TBcutoff", force = TRUE)
```

To install CRAN version:
```{r}
install.packages("mixedCCA")
```

The comparison of graphical model estimation on QMP data also requires **SPRING** R package
```{r}
devtools::install_github("GraceYoon/SPRING", force = TRUE)
```
### Function scripts

* **Twovariable_simulation_functions.R** + **Functions** folder - Functions to implement ML method (ORG and MLBD are available from mixedCCA R package).

### Figure 1 scripts

* **3dfigure1.R** - 3D-figure of bridge inverse function for TC case (Figure 1 (left) in the manuscript)

* **Delta_vs_Zeroprop.R** - generate a figure showing the relationship between latent cutoff Delta and zero proportion of data. This is the figure of the inverse of normal cdf (Figure 1 (right) in the manuscript)

### Simulation scripts

#### Simulated data with 2 variables

* **Twovariable_simulation.R** - R code to run simulation studies with two variables that compare computation time and approximation errors for 5 variable combinations (TC, TT, BC, BB, and TB) across three methods: ORG (original), ML (linear interpolation) and MLBD (hybrid) (Section 4.1.1, Section 4.2 and Supplement Section S.2 in the manuscript)

* **Twovariable_simulation.Rout** - Saved console output of running **Twovariable_simulation.R.**



* **xtable_figure_for_TwoVariable_simulation.R** - Create latex table of computation times and visualize the error of two variable simulation (Table 1, Figure 3, and Supplement Figures S1-S5 in the manuscript)

#### Comparisons on real data examples

* **Realdata3_tab.R** - Compare accuracy and computation time on three real datasets: QMP, AGP, TCGA-BRCA (Sections 4.1.2 and 4.2 in the manuscript)

* **Realdata3_tab.Rout** - Saved console output of running **Realdata3.R**

* **Realdata-SPRING.R** - Compare accuracy and computation time of graphical model estimation on QMP data (Section 4.1.3 in the manuscript)

* **Realdata-SPRING.Rout** - Saved console output of running **Realdata-SPRING.R**

* **xtable_for_Realdata.R** - Create latex table of results on three real datasets: QMP, AGP, TCGA-BRCA, and SPRING on QMP data. (Table 2 in the manuscript)

* **RunTime_v2.R** - Run time versus dimension in log10 scale: comparing ORG, MLBD and Kendall only.

* **RunTimePlot_v2_loglog.R** - Using saved results, create data.frame format and draw a figure (Figure 4 in the manuscript)



### Precomputed results in **Data** folder:

#### Grid for Figure 1

* **tcInv_equigrid.Rda** - For Figure 1 in manuscript (plot_ly function works with equidistant grid values)

#### Results from two variable simulations 

* **TwoSim_TC_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is truncated continuous and the other is continuous type.

  **TwoSim_TT_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where both variables are truncated continuous type.

  **TwoSim_BC_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is binary and the other is continuous type.

  **TwoSim_BB_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where both variables are binary type.

  **TwoSim_TB_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is truncated continuous and the other is binary type.

#### American Gut Project data

<!-- * **AGP_prunedata.R** - how to load and clean Americun Gut Project data.-->

* **amgutpruned.rdata** - American Gut Project microbiome data that was used for the analyses.
  
#### TCGA data

* **TCGA_dataloading.R** - R script to process gene expression and microRNA data from TCGA-BRCA data portal. See the comments at the top of the script on how to load the original data for subsequent processing.

* **matchedTCGA.Rdata** - processed gene expression and microRNA TCGA-BRCA data used for the analyses.

#### Results from numerical comparisons on real data

* **Realdata3_tab.Rda** - Saved computation results of Kendall's tau and latent correlation of three real data: QMP, AGP, TCGA-BRCA.

*  **Realdata-SPRING_tab.Rda** - Saved graphical model estimation result on QMP data.

* **RunTimePlot_range_v2_100.Rda** - saved results to plot run time versus dimension in log10 scale for the sample size 100.

*  **RunTimePlot_range_v2_6482.Rda** - saved results to plot run time versus dimension in log10 scale for the sample size 6482.






