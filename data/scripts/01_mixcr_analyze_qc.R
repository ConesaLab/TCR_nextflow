#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

outDir = "TCRanalysis_bookdown"

library(ggplot2)
library(gridExtra)

patients = read.table("SampleInfo.csv", sep = ",", header = T)
patients = patients[order(nchar(patients$CRG.code), patients$CRG.code),]
patients$Group = factor(patients$Group, levels = c("control", "withoutMHE", "withMHE"))


# Split the complete report in 5 subreports and then select number 00 (alignment) and number 04 (assemble)
files = list.files(path = args[1], pattern="*.report", full.names = T)
files = files[order(nchar(files), files)]
files = files[patients_ordered$files_order]
for file in *.report; do
  csplit -f $file $file '/^==/' '{*}'
  rm *.report05
done



# Read Alignment Report 

temp = list.files(pattern="*report00")
temp = temp[order(nchar(temp), temp)]
# Concatenate all reports in one
alignment_reports = setNames(do.call(cbind, lapply(temp, read.table, row.names = 1, sep=":", comment.char="=", skip = 6)), 
                             nm = patients$SampleID)

# Remove (%) values
alignment_reports = apply(alignment_reports, 1, stringr::str_remove, pattern = " [(\\)].+")
# Change to numeric
alignment_reports = apply(alignment_reports, 2, as.numeric)
# Change CRG id by Clinic id
rownames(alignment_reports) = patients$SampleID

# Order alphabetically + by group
alignment_reports = alignment_reports[order(patients$SampleID),]
patients = patients[order(patients$SampleID),]
patients = patients[order(patients$Group),]
alignment_reports = alignment_reports[patients$SampleID,]

# Create alignment report table
alignment_reports = data.frame(alignment_reports, check.names = F)
alignment_reports[,"patientID"] = rownames(alignment_reports)
alignment_reports = alignment_reports[c(24,1:23)]


## Auxiliary functions

color = function(x, palette=0, reverse=1, alpha=1){
  library("viridis")
  library(RColorBrewer)
  if(palette==0)
    return(c(viridis(x, direction=reverse, alpha = alpha)))
  if(palette==1)
    return(c(magma(x, direction=reverse, alpha = alpha)))
  if(palette==2)
    return(c(plasma(x, direction=reverse, alpha = alpha)))
  if(palette==3)
    return(c(inferno(x, direction=reverse, alpha = alpha)))
  if(palette==4)
    return(c(cividis(x, direction=reverse, alpha = alpha)))
  if(palette==5)
    return(c(rainbow(x)))
  if(palette==6)
    return(c(RColorBrewer::brewer.pal(x, "Set3")))
}


## Plots

df_matrix = tibble::as_tibble(alignment_reports)

# Plot
patient_color = c("orange", "olivedrab3", "mediumpurple")[factor(patients$Group)]
df_matrix$patientID = factor(df_matrix$patientID, levels=df_matrix$patientID)


df_gather_alignment = tidyr::gather(data = df_matrix[,c(1,3,5,6,7)], key = "parameter", value = "reads", -1)
df_gather_alignment$parameter = factor(df_gather_alignment$parameter, levels = unique(df_gather_alignment$parameter))

p1 = ggplot(df_gather_alignment, aes(fill = parameter, y = reads, x = patientID)) + 
 geom_bar(position="stack", stat="identity") +
  viridis::scale_fill_viridis(discrete = TRUE, option = "D") +
  ggtitle("MiXCR alignment report") +
  ylab("Number of reads") +
  theme(axis.text.x = element_text(angle = 90, color = patient_color))
p1

df_gather_chains = tidyr::gather(data = df_matrix[,c(1,14:18,20,21)], key = "parameter", value = "reads", -1)
df_gather_chains$parameter = factor(df_gather_chains$parameter, levels = unique(df_gather_chains$parameter))
colnames(df_gather_chains)[1] = "patient"

p2 = ggplot(df_gather_chains, aes(fill = parameter, y = reads, x = patient)) + 
 geom_bar(position="stack", stat="identity") +
  ggtitle("MiXCR alignment report") +
    ylab("Number of reads") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, color = patient_color)) +
    scale_fill_manual(values = color(length(unique(df_gather_chains$parameter)), palette = 6))
p2


# plots
png(file.path(outDir, "01_alignment.png"), width = 30, height = 15, units = "cm", res = 300)
p2
dev.off()



# Assemble report

patients = read.table("SampleInfo.csv", sep = ",", header = T)
patients = patients[order(nchar(patients$CRG.code), patients$CRG.code),]
patients$Group = factor(patients$Group, levels = c("control", "withoutMHE", "withMHE"))

temp = list.files(pattern="*report04")
temp = temp[order(nchar(temp), temp)]
# Concatenate all reports in one
assemble_reports = setNames(do.call(cbind, lapply(temp, read.table, row.names = 1, sep=":", comment.char="=", skip = 7)), 
                             nm = patients$SampleID)

# Remove (%) values
assemble_reports = apply(assemble_reports, 1, stringr::str_remove, pattern = " [(\\)].+")
# Change to numeric
assemble_reports = apply(assemble_reports, 2, as.numeric)
# Change crgID by clinicID
rownames(assemble_reports) = patients$SampleID

# Order alphabetically + by group
assemble_reports = assemble_reports[order(patients$SampleID),]
patients = patients[order(patients$SampleID),]
patients = patients[order(patients$Group),]
assemble_reports = assemble_reports[patients$SampleID,]

rownames(assemble_reports) == patients$SampleID

# Create alignment report table
assemble_reports = data.frame(assemble_reports, check.names = F)
assemble_reports[,"patientID"] = rownames(assemble_reports)
assemble_reports = assemble_reports[c(23,1:22)]


## Plots

df_matrix_assemble = data.frame(assemble_reports, check.names = F)
df_matrix_assemble[,"patient"] = rownames(df_matrix_assemble)
df_matrix_assemble = df_matrix_assemble[c(ncol(df_matrix_assemble),1:ncol(df_matrix_assemble)-1)]
df_matrix_assemble = tibble::as_tibble(df_matrix_assemble)

# Plot
patient_color = c("orange", "olivedrab3", "mediumpurple")[factor(patients$Group)]
df_matrix_assemble$patient = factor(df_matrix_assemble$patient, levels=df_matrix_assemble$patient)


df_gather_chains = tidyr::gather(data = df_matrix_assemble[,c(1,18:24)], key = "parameter", value = "reads", -1)
df_gather_chains$parameter = factor(df_gather_chains$parameter, levels = unique(df_gather_chains$parameter))

p3 = ggplot(df_gather_chains, aes(fill = parameter, y = reads, x = patient)) + 
 geom_bar(position="stack", stat="identity") +
  ggtitle("MiXCR assemble report") +
    ylab("Number of clonotypes") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, color = patient_color)) +
    scale_fill_manual(values = color(length(unique(df_gather_chains$parameter)), palette = 6))
p3


# plots
png(file.path(outDir, "01_assemble.png"), width = 30, height = 15, units = "cm", res = 300)
p3
dev.off()

