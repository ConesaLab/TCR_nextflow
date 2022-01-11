#!/usr/bin/env Rscript


# Packages
library(vegan)
library(microbiome)


# # Define current/input/output directories
# inputDir <- "results"
# dir.create(file.path(inputDir, "TCRanalysis_bookdown"))
# outputDir <- file.path(inputDir, "TCRanalysis_bookdown")

# Load patient data
patients = read.table("sampleslist.csv", sep = ",", header = T)
patients = patients[order(nchar(patients$SampleID), patients$SampleID),]
patients$Group = as.factor(patients$Group)

patients_ordered =  patients[order(patients$SampleID),]
patients_ordered = patients_ordered[order(patients_ordered$Group),]


# Clones from mixcr output

## TCR clonotypes
# List input files
files = list.files(path = inputDir, pattern="*.clonotypes.ALL.txt", full.names = T)
files = files[order(nchar(files), files)]
files = files[patients_ordered$files_order]

# Create a list of clone tables
clonotypes_all = setNames(lapply(files, read.delim), nm = patients_ordered$SampleID)

# exclude everything that is not TCR
clonotypes_tcr = lapply(clonotypes_all, function(i){ 
  i = i[-grep("IG", i$allVHitsWithScore),]
  })

# Test minimum CDR3 aa's sequence length
lapply(names(clonotypes_tcr), function(i){min(nchar(clonotypes_tcr[[i]]$aaSeqCDR3))})

# Summary of filter criteria
summary_table_tcr = do.call(rbind, lapply(names(clonotypes_tcr), function(i){ 
  table = clonotypes_tcr[[i]]
  sample = i
  TCRclones = nrow(table)
  filterby_readcount = nrow(table[table$cloneCount >= 2,])
  min_CDR3aa = min(nchar(clonotypes_tcr[[i]]$aaSeqCDR3))
  filterby_CDR3aa = nrow(table[nchar(table$aaSeqCDR3) >= 4,])
  filterby_both = nrow(table[table$cloneCount >= 2 & nchar(table$aaSeqCDR3) >= 4,])
  cbind(sample,TCRclones,filterby_readcount,min_CDR3aa,filterby_CDR3aa,filterby_both)
  }))


# Filtering and output
clones_tcr_filtered = lapply(clonotypes_tcr, function(i){ 
  #Filter
  i = i[i$cloneCount >= 2,]
  #Recalculate clone Fraction (recalculatedcloneFraction = 1)
  recalculatedcloneFraction = microbiome::transform(i[,"cloneCount",drop=F], "compositional")
  colnames(recalculatedcloneFraction) = "recalculatedcloneFraction"
  cbind(i, recalculatedcloneFraction)
})

save(clones_tcr_filtered, file = "clones_tcr_filtered.Rda")



<!-- ## TRA clonotypes -->
<!-- ```{r, results='hide'} -->
<!-- # List input files -->
<!-- files = list.files(path = inputDir, pattern="*.clonotypes.ALL.txt", full.names = T) -->
<!-- files = files[order(nchar(files), files)] -->
<!-- files = files[patients_ordered$files_order] -->

<!-- # Create a list of clone tables -->
<!-- clonotypes_all = setNames(lapply(files, read.delim), nm = patients_ordered$SampleID) -->

<!-- # exclude everything that is not TRA -->
<!-- clonotypes_TRA = lapply(clonotypes_all, function(i){  -->
<!--   i = i[grep("TRA", i$allVHitsWithScore),] -->
<!--   }) -->

<!-- # Test minimum CDR3 aa's sequence length -->
<!-- lapply(names(clonotypes_TRA), function(i){min(nchar(clonotypes_TRA[[i]]$aaSeqCDR3))}) -->

<!-- # Summary of filter criteria -->
<!-- summary_table_TRA = do.call(rbind, lapply(names(clonotypes_TRA), function(i){  -->
<!--   table = clonotypes_TRA[[i]] -->
<!--   sample = i -->
<!--   TRAclones = nrow(table) -->
<!--   filterby_readcount = nrow(table[table$cloneCount >= 2,]) -->
<!--   min_CDR3aa = min(nchar(clonotypes_TRA[[i]]$aaSeqCDR3)) -->
<!--   filterby_CDR3aa = nrow(table[nchar(table$aaSeqCDR3) >= 4,]) -->
<!--   filterby_both = nrow(table[table$cloneCount >= 2 & nchar(table$aaSeqCDR3) >= 4,]) -->
<!--   cbind(sample,TRAclones,filterby_readcount,min_CDR3aa,filterby_CDR3aa,filterby_both) -->
<!--   })) -->

<!-- # Filtering and output -->
<!-- clones_TRA_filtered = lapply(clonotypes_TRA, function(i){  -->
<!--   #Filter -->
<!--   i = i[i$cloneCount >= 2,] -->
<!--   #Recalculate clone Fraction (recalculatedcloneFraction = 1) -->
<!--   recalculatedcloneFraction = microbiome::transform(i[,"cloneCount",drop=F], "compositional") -->
<!--   colnames(recalculatedcloneFraction) = "recalculatedcloneFraction" -->
<!--   cbind(i, recalculatedcloneFraction) -->
<!-- }) -->

<!-- save(clones_TRA_filtered, file = file.path(currentDir, "clones_TRA_filtered.Rda")) -->

<!-- ``` -->

<!-- ## TRB clonotypes -->
<!-- ```{r, results='hide'} -->
<!-- # List input files -->
<!-- files = list.files(path = inputDir, pattern="*.clonotypes.ALL.txt", full.names = T) -->
<!-- files = files[order(nchar(files), files)] -->
<!-- files = files[patients_ordered$files_order] -->

<!-- # Create a list of clone tables -->
<!-- clonotypes_all = setNames(lapply(files, read.delim), nm = patients_ordered$SampleID) -->

<!-- # exclude everything that is not TRB -->
<!-- clonotypes_TRB = lapply(clonotypes_all, function(i){  -->
<!--   i = i[grep("TRB", i$allVHitsWithScore),] -->
<!--   }) -->

<!-- # Test minimum CDR3 aa's sequence length -->
<!-- lapply(names(clonotypes_TRB), function(i){min(nchar(clonotypes_TRB[[i]]$aaSeqCDR3))}) -->

<!-- # Summary of filter criteria -->
<!-- summary_table_TRB = do.call(rbind, lapply(names(clonotypes_TRB), function(i){  -->
<!--   table = clonotypes_TRB[[i]] -->
<!--   sample = i -->
<!--   TRBclones = nrow(table) -->
<!--   filterby_readcount = nrow(table[table$cloneCount >= 2,]) -->
<!--   min_CDR3aa = min(nchar(clonotypes_TRB[[i]]$aaSeqCDR3)) -->
<!--   filterby_CDR3aa = nrow(table[nchar(table$aaSeqCDR3) >= 4,]) -->
<!--   filterby_both = nrow(table[table$cloneCount >= 2 & nchar(table$aaSeqCDR3) >= 4,]) -->
<!--   cbind(sample,TRBclones,filterby_readcount,min_CDR3aa,filterby_CDR3aa,filterby_both) -->
<!--   })) -->

<!-- # Filtering and output -->
<!-- clones_TRB_filtered = lapply(clonotypes_TRB, function(i){  -->
<!--   #Filter -->
<!--   i = i[i$cloneCount >= 2,] -->
<!--   #Recalculate clone Fraction (recalculatedcloneFraction = 1) -->
<!--   recalculatedcloneFraction = microbiome::transform(i[,"cloneCount",drop=F], "compositional") -->
<!--   colnames(recalculatedcloneFraction) = "recalculatedcloneFraction" -->
<!--   cbind(i, recalculatedcloneFraction) -->
<!-- }) -->

<!-- save(clones_TRB_filtered, file = file.path(currentDir, "clones_TRB_filtered.Rda")) -->
<!-- ``` -->


<!-- # Shannon Evenness (SE) -->
<!-- Common diversity indices are special cases of Rényi diversity: -->
<!-- $$ -->
<!-- H_a = \frac{1}{1-a} log_b \sum_{i = 1}^{n}{f_i^a} -->
<!-- $$ -->
<!-- where $a$ is a scale parameter, $f_i$ is the frequency distribution (proportional abundance of species/clones) and $b$ is the base of the logarithm. It is most popular to use natural logarithms, but some argue for base b = 2 (which makes sense, but no real difference). -->

<!-- Hill (1975) suggested to use so-called ‘Hill numbers’ defined as $N_a = exp(H_a)$: -->
<!-- $$ -->
<!-- N_a = (\sum_{i = 1}^{n}{f_i^a}) ^ \frac{1}{1-a} -->
<!-- $$ -->
<!-- Some Hill numbers are: -->
<!-- - **Richness**: the number of species with $a = 0$ or $exp(H')$  -->
<!-- - inverse Simpson with a = 2 and  -->
<!-- - $1/ max(p_i)$ with $a = ∞$ -->
<!-- - the exponent of Shannon diversity with $a = 1$. -->

<!-- Diversity is not defined for the case alpha = 1. However, we use L'Hospital's rule to fiind that as alpha tends to 1, Diversity tends to the Shannon entropy. Thus, the Shannon entropy is a special case of the Diversity for alpha = 1. Shannon or Shannon–Weaver (or Shannon–Wiener) index is defined as $H' = − \sum_{i} f_i log_b f_i$ -->

<!-- Diversity profile (${}^{a}D$) is defined as: -->
<!-- $$ -->
<!-- {}^{a}D = SR \times {}^{a}E -->
<!-- $$ -->
<!-- where ${}^{a}E$ is the **evenness** $SR$ is the species **richness** ($SR={}^{a=0}D$), that is, the number of unique clones in a repertoire dataset. -->

<!-- Shannon-Evenness is 1 if all clones in a repertoire have the same frequency (an "even" repertoire), or it converges to 0 if very few clones dominate in the repertoire ("polarized" repertoire). -->

<!-- ## Manual calculation -->
<!-- ```{r} -->
<!-- f = clones_tcr_filtered$CR199[,"recalculatedcloneFraction", drop = F] -->

<!-- ## SHANNON ENTROPY (H') = D(a=1) -->
<!-- ###Option1 -->
<!-- -sum(f * log(f)) -->
<!-- ###Option2 -->
<!-- microbiome::diversity(x = f, index = "shannon") -->
<!-- ###Option3 -->
<!-- vegan::renyi(x = f, scales = 1) -->

<!-- ## SPECIES RICHNESS (SR) calculated by Hill formula -->
<!-- ###Option1: number of rows -->
<!-- nrow(f) -->
<!-- ###Option2: package microbiome -->
<!-- microbiome::richness(x = f, index = "observed") -->
<!-- ###Option3: Hill formula with a=0 -->
<!-- alpha = 0 -->
<!-- hill = (sum(f^alpha)) ^ (1/(1-alpha)) -->
<!-- hill -->
<!-- ###Option4: exp(Rényi formula) with a=0 -->
<!-- alpha = 0 -->
<!-- renyi = (1/(1-alpha)) * log(sum(f^alpha)) -->
<!-- exp(renyi) -->

<!-- # SpeciesRichness (SR) calculated by Rényi formula -->
<!-- ###Option1 -->
<!-- vegan::renyi(x = f, scales = 0) -->
<!-- ###Option2 -->
<!-- alpha = 0 -->
<!-- (1/(1-alpha)) * log(sum(f^alpha)) -->

<!-- ## SHANNON EVENESS (S-E) = D / SR -->
<!-- SpeciesRichness = nrow(f) -->
<!-- ShannonEntropy = -sum(f * log(f)) -->
<!-- ShannonEvenness = ShannonEntropy / SpeciesRichness -->
<!-- ShannonEvenness -->

<!-- SpeciesRichness = microbiome::richness(x = f, index = "observed") -->
<!-- ShannonEntropy = microbiome::diversity(x = f, index = "shannon") -->
<!-- ShannonEvenness = ShannonEntropy / SpeciesRichness -->
<!-- ShannonEvenness -->
<!-- ``` -->


<!-- ## TCR - SE -->
<!-- ```{r} -->
<!-- TCR_SE <- do.call(rbind, lapply(clones_tcr_filtered, function(x){ -->
<!--   diversity = exp(vegan::renyi(x = x[,"recalculatedcloneFraction"], scales = 1)) -->
<!--   richness = exp(vegan::renyi(x = x[,"recalculatedcloneFraction"], scales = 0)) -->
<!--   #shannon evenness -->
<!--   SE = diversity/richness -->
<!--   })) -->
<!-- colnames(TCR_SE) <- "shannon" -->

<!-- TCR_SE -->
<!-- ``` -->

<!-- ## TRA SE -->
<!-- ```{r} -->
<!-- TRA_SE <- do.call(rbind, lapply(clones_TRA_filtered, function(x){ -->
<!--   diversity = exp(vegan::renyi(x = x[,"recalculatedcloneFraction"], scales = 1)) -->
<!--   richness = exp(vegan::renyi(x = x[,"recalculatedcloneFraction"], scales = 0)) -->
<!--   #shannon evenness -->
<!--   SE = diversity/richness -->
<!--   })) -->
<!-- colnames(TRA_SE) <- "shannon" -->

<!-- TRA_SE -->
<!-- ``` -->

<!-- ## TRB SE -->
<!-- ```{r} -->
<!-- TRB_SE <- do.call(rbind, lapply(clones_TRB_filtered, function(x){ -->
<!--   diversity = exp(vegan::renyi(x = x[,"recalculatedcloneFraction"], scales = 1)) -->
<!--   richness = exp(vegan::renyi(x = x[,"recalculatedcloneFraction"], scales = 0)) -->
<!--   #shannon evenness -->
<!--   SE = diversity/richness -->
<!--   })) -->
<!-- colnames(TRB_SE) <- "shannon" -->

<!-- TRB_SE -->
<!-- ``` -->


<!-- # Output files -->
<!-- ```{r} -->
<!-- # Summary table -->
<!-- summary_table_tcr = cbind(summary_table_tcr, TCR_SE) -->
<!-- summary_table_TRA = cbind(summary_table_TRA, TRA_SE) -->
<!-- summary_table_TRB = cbind(summary_table_TRB, TRB_SE) -->

<!-- write.table(summary_table_tcr, file = file.path(outputDir, "02_filtering&evenness_summarytable_TCR.csv"), sep = "\t", row.names = F, col.names = T, quote = F) -->
<!-- write.table(summary_table_TRA, file = file.path(outputDir, "02_filtering&evenness_summarytable_TRA.csv"), sep = "\t", row.names = F, col.names = T, quote = F) -->
<!-- write.table(summary_table_TRB, file = file.path(outputDir, "02_filtering&evenness_summarytable_TRB.csv"), sep = "\t", row.names = F, col.names = T, quote = F) -->

<!-- # Filtered files -->
<!-- lapply(names(clones_tcr_filtered), function(x){ -->
<!--   write.table(clones_tcr_filtered[[x]],  -->
<!--               file = file.path(currentDir, "clones_TCR_filtered", paste0(x, ".txt")), -->
<!--               quote = F, sep = "\t", row.names = F) -->
<!-- }) -->

<!-- lapply(names(clones_TRA_filtered), function(x){ -->
<!--   write.table(clones_TRA_filtered[[x]],  -->
<!--               file = file.path(currentDir, "clones_TRA_filtered", paste0(x, ".txt")), -->
<!--               quote = F, sep = "\t", row.names = F) -->
<!-- }) -->

<!-- lapply(names(clones_TRB_filtered), function(x){ -->
<!--   write.table(clones_TRB_filtered[[x]],  -->
<!--               file = file.path(currentDir, "clones_TRB_filtered", paste0(x, ".txt")), -->
<!--               quote = F, sep = "\t", row.names = F) -->
<!-- }) -->
<!-- ``` -->


