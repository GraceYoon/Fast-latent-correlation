# Fast-latent-correlation


R code with functions that were used for simulations and data analysis in the paper ["Fast computation of latent correlations" by Grace Yoon, Christian L. M&uuml;ller and Irina Gaynanova](https://arxiv.org/abs/2006.13875).


### Simulation scripts

**3dfigure1.R** - 3D-figure of bridge inverse function for TC case.


**Delta_vs_Zeroprop.R** - generate a figure showing the relationship between latent cutoff Delta and zero proportion of data. Basically inverse of normal cdf.


**Twovariable_simulation.R** - Simulation to save results of computation time and errors between two variables in five different type combinations. Compared three methods: ORG, ML and MLBD.

**Twovariable_simulation.Rout** - Saved outputs of Twovariable_simulation.R.

**Twovariable_simulation_functions.R** and two files in folder **Functions**- Functions for ML method only.

**xtable_figure_for_TwoVariable_simulation.R** - Create latex table of computation times and visualize the error of two variable simulation.


**Realdata3_tab.R** - Check accuracy and computation time of three real data: QMP, AGP, TCGA-BRCA.

**Realdata3_tab.Rout** - Saved outputs of Realdata3.R

**Realdata-SPRING.R** - Check accuracy and computation time of graphical model estimation on QMP data.

**Realdata-SPRING.Rout** - Saved outputs of Realdata-SPRING.R

**xtable_for_Realdata.R** - Creat latex table of three real data results: QMP, AGP, TCGA-BRCA, and SPRING on QMP data.


**RunTime_v2.R** - Run time versus dimension in log10 scale: comparing ORG, MLBD and Kendall only.
**RunTimePlot_v2_loglog.R** - Using saved results, create data.frame format and draw a figure.




### Supporting data and precomputed results are in the folder **Data**:

**tcInv_equigrid.Rda** - For figure 1 in manuscript. (plot_ly function seems to be able to handle only equigrid data values)


**TwoSim_TC_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is truncated continuous and the other is continuous type.

**TwoSim_TT_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where both variables are truncated continuous type.

**TwoSim_BC_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is binary and the other is continuous type.

**TwoSim_BB_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where both variables are binary type.

**TwoSim_TB_rep10.Rda** - Saved computation time of three methods (ORG, ML and MLBD) and absolute errors of ML and MLBD method compared to ORG method where one variable is truncated continuous and the other is binary type.


**AGP_prunedata.R** - how to load and clean Americun Gut Project data.

**amgutpruned.rdata** - American Gut Project microbiome data.

**TCGA_dataloading.R** - how to load and clean gene expression and microRNA data from TCGA-BRCA data portal.

**matchedTCGA.Rdata** - saved gene expression and microRNA from TCGA-BRCA.


**Realdata3_tab.Rda** - Saved computation results of Kendall's tau and latent correlation of three real data: QMP, AGP, TCGA-BRCA.

**Realdata-SPRING_tab.Rda** - Saved graphical model estimation result on QMP data.


**RunTimePlot_range_v2_100.Rda** - saved results to plot run time versus dimension in log10 scale for the sample size 100.

**RunTimePlot_range_v2_6482.Rda** - saved results to plot run time versus dimension in log10 scale for the sample size 6482.






