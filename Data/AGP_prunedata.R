# Import American Gut Data and prune the data.
# Adapted from procAG_metadata.R file (code by Christain M\"{u}ller) in Dropbox Shared folder "AmericanGut".
# Iteratively pruned data (step 1 and 2) instead use all step 1, 2 and 3.


library(phyloseq)

# Sys.setenv('R_MAX_VSIZE'=1000000000000) # not worked. So run in the R1 server.

system.time(ag <- import_biom("8237_analysis.biom")) # due to large file size, it takes around 97 seconds.
map <- read.delim("8237_analysis_mapping.txt", sep="\t", header=TRUE, row.names=1) ## import metadata from mapping file
sample_data(ag) <- map ## assign metadata to phyloseq object

## All fecal data
ag.fe  <- subset_samples(ag, body_habitat=="UBERON:feces") ## only fecal samples 

ag.fe
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 27116 taxa and 8440 samples ]
# sample_data() Sample Data:       [ 8440 samples by 517 sample variables ]
# tax_table()   Taxonomy Table:    [ 27116 taxa by 7 taxonomic ranks ]
# ROW = OTUs, COLUMN = samples

dim(ag.fe@otu_table@.Data) # ROW=OTU(p), COLUMN=Sample(n)
totalsum <- colSums(ag.fe@otu_table@.Data) # total abundance for each sample


#####################################
# Prune step 1 (for each sample): sequencing depths (total abundance) > 10,000 
# Prune step 2 (for taxa/bacteria/species): taxa present in at least 30% of samples
# Iteratively prune step 1 and 2 until convergence OR Prune step 3: check sequencing depths again and prune small 10% quantile.


# STEP 1
## Prune samples (for the purpose of covariance estimation, e.g.)
depths <- colSums(ag.fe@otu_table@.Data) ## calculate sequencing depths
# > summary(depths)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1   14112   22206   30331   34695  547447 

## We usually do this BUT for the purpose of the paper you may want to use all samples
ag.filt1 <- prune_samples(depths > 10000, ag.fe) ## Only samples above 10000 sequencing depth
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 27116 taxa and 7203 samples ]
# sample_data() Sample Data:       [ 7203 samples by 517 sample variables ]
# tax_table()   Taxonomy Table:    [ 27116 taxa by 7 taxonomic ranks ]


# STEP 2
## We usually do this BUT for the purpose of the paper you may want to use all taxa
freq <- rowSums(sign(ag.filt1@otu_table@.Data))
# > summary(freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     2.0   138.7    17.0  7138.0 

# > nsamples(ag.filt1)
# [1] 7203
# > .3*nsamples(ag.filt1)
# [1] 2160.9
ag.filt2 <- prune_taxa(freq > .3*nsamples(ag.filt1), ag.filt1) ## more pruning (taxa present in at least 30% of samples)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 481 taxa and 7203 samples ]
# sample_data() Sample Data:       [ 7203 samples by 517 sample variables ]
# tax_table()   Taxonomy Table:    [ 481 taxa by 7 taxonomic ranks ]


##### We decided to do step 3 instead of doing step 1 and 2 iteratively.
# If we do STEP 3,
## Only samples in 10% quantile
depths   <- colSums(ag.filt2@otu_table@.Data)
ag.filt3 <- prune_samples(depths > quantile(depths, probs=seq(0, 1, .1))[2], ag.filt2)
## Filtered phyloseq object
ag.filt3
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 481 taxa and 6482 samples ]
# sample_data() Sample Data:       [ 6482 samples by 517 sample variables ]
# tax_table()   Taxonomy Table:    [ 481 taxa by 7 taxonomic ranks ]
summary(colSums(ag.filt3@otu_table@.Data))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 10426   15915   21856   30025   33664  512384
summary(rowSums(sign(ag.filt3@otu_table@.Data)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1929    2445    3094    3384    4117    6430
# .3*nsamples(ag.filt3)
# [1] 1944.6

amgutpruned <- t(ag.filt3@otu_table@.Data)
# ROW = samples, COLUMN = OTUs
# save(amgutpruned, file = "Data/amgutpruned.rdata")