# Fast-latent-correlation


R code with functions that were used for simulations and data analysis in the paper "Fast computation of latent correlations" by Grace Yoon, Christian L. M\"{u}ller and Irina Gaynanova.


### Simulation scripts

**3dfigure1.R** - 3D-figure of bridge inverse function for TC case.

**Twovariable_simulation.R** - Simulation to save results of computation time and errors between two variables in all different type combinations. Tried 

**Twovariable_simulation_functions.R** - Wrapper functions for the simulation between two synthetic variables.
	
**xtable_figure_for_TwoVariable_simulation.R** - Visualize Check computation time of three real data: QMP, AGP, TCGA-BRCA, and graphical model estimation on QMP data

**RunTimePlot_loglog.R** - Run time versus dimension in log10 scale: comparing ORG, MLBD and Kendall only.

**accuracy_realdata.R** - Check accuracy of three real data: QMP, AGP, TCGA-BRCA, and graphical model estimation on QMP data

**microbenchmark_realdata.R** - Check computation time of three real data: QMP, AGP, TCGA-BRCA, and graphical model estimation on QMP data


### Supporting data and precomputed results are in the folder **Data**:

**tcInv_equigrid.Rda** - For figure 1 in manuscript. (plot_ly function seems to be able to handle only equigrid data values)


**TwoSimVariableTC.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is truncated continuous and the other is continuous type.

**TwoSimVariableTT.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where both variables are truncated continuous type.

**TwoSimVariableBC.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is binary and the other is continuous type.

**TwoSimVariableBB.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where both variables are binary type.

**TwoSimVariableTB.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is truncated continuous and the other is binary type.


**Kendallmatrix_QMP_org.Rda** - Calculated Kendall's tau correlation of QMP data.

**Kendallmatrix_amgut_org.Rda** - Calculated Kendall's tau correlation of AGP data.

**Kendallmatrix_TCGA_org.Rda** - Calculated Kendall's tau correlation of TCGA-BRCA data.


**Kcor_TCGA_fast.Rda** - Calculated latent rank-based correlation of TCGA-BRCA data using MLBD method.

**Kcor_QMP_org.Rda** - Calculated latent rank-based correlation of QMP data using ORG method.

**Kcor_QMP_fast.Rda** - Calculated latent rank-based correlation of QMP data using MLBD method.

**Kcor_amgut_org.Rda** - Calculated latent rank-based correlation of AGP data using ORG method.

**Kcor_amgut_fast.Rda** - Calculated latent rank-based correlation of AGP data using MLBD method.

**Kcor_TCGA_org.Rda** - Calculated latent rank-based correlation of TCGA-BRCA data using ORG method.

**Kcor_TCGA_fast.Rda** - Calculated latent rank-based correlation of TCGA-BRCA data using MLBD method.


**RunTimePlot_range_100.Rda** - saved results to plot run time versus dimension in log10 scale for the sample size 100.

**RunTimePlot_range_6482.Rda** - saved results to plot run time versus dimension in log10 scale for the sample size 6482.


**amgutpruned.rdata** - American Gut Project microbiome data.

**matchedTCGA.Rdata** - saved gene expression and microRNA from TCGA-BRCA.

**DataLoading.R** - how to load and clean gene expression and microRNA data from TCGA-BRCA data portal.



