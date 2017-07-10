############
#RAW-HAMMER#
############

#Written by Richard Allen White III and Paul D Piehowski
#Created on Dec 1, 2014
#R script for Spectral count Metaproteomics 
#Updated Feb 4, 2016

#Load libraries
library(plyr)
library(reshape2)
library(ggplot2)
library(sqldf)
library(tcltk)

#Load files
inputFiles <- c('.txt')

massSpecOutput <- read.delim(inputFiles[1], header = T)

#create filename for naming output files by replacing '.[stuff]' with ''
filename_s <- gsub('\\..*$', '', inputFiles[1]) 

# Filter by DelMass PPM +/-  ppm and calculate FDR: %
PPM_cutoffs <- #add cut off here
SpecProb_cutoffs <- #add cut off here
  
massSpecOutput_trim <- subset(massSpecOutput, 
                              DelM_PPM >= -1*PPM_cutoffs & DelM_PPM <= PPM_cutoffs & 
                                MSGF_SpecProb <= SpecProb_cutoffs)

ReverseSpec <- 1*sum(grepl('XXX_', massSpecOutput_trim$Protein))
FDR <- ReverseSpec/(length(massSpecOutput_trim$Protein)-ReverseSpec)*100


#As ggplot2 
p1 <- ggplot(massSpecOutput_trim, aes(x=DelM_PPM))
p1 <- p1 + geom_histogram(binwidth=0.1, fill="blue", col="black")
p1 <- p1 + theme_bw() 
p1 <- p1 + ggtitle("")+ 
  theme(plot.title = element_text(lineheight=2, face="bold", size=20))

#Print ggplot2 to png for DelMass Error
png(".png", width = 15, height = 15, units = 'in', res = 300)
p1
dev.off()

#standard graphics 
hist(massSpecOutput_trim$DelM_PPM, breaks = 1000, xlim = c(-10,10),
     main = "Mass Error Histogram 
     "
)

##Remove reverse spectra
massSpecOutput_trim <- massSpecOutput_trim[!grepl('XXX_', massSpecOutput_trim$Protein),]

# Make a table that can be used for calculating summary statistics 
peptide <- sqldf("Select Peptide, Protein, NTT
                 From massSpecOutput_trim
                 Group by Peptide")

unique_peptides_all <- 1*length(peptide$Peptide)

x <- sqldf("Select Protein 
           From peptide 
           Group by Protein")
unique_proteins_all <- 1*length(x$Protein)

x <- subset(peptide, NTT == 2)
unique_tryptic_all <- 1*length(x$NTT)

x <- subset(peptide, NTT == 1)
unique_parttryptic_all <- 1*length(x$NTT)

x <- subset(peptide, NTT == 0)
unique_nontryptic_all <- 1*length(x$NTT)

summary_table <- data.frame(rbind(unique_peptides_all, unique_proteins_all, unique_tryptic_all, 
                                  unique_parttryptic_all, unique_nontryptic_all))

name <- paste(filename_s, "summary_table.txt", sep = "_")
write.table(summary_table, file = name, quote = F, sep = "\t", col.names = F,
            row.names = c("Unique Peptides", "Unique Proteins", "Tryptic", "Partially Tryptic",
                          "Non-Tryptic"))


# Pull unique peptide row for calculating peptide characteristics
peptide <- as.character(peptide$Peptide)

# Export peptide as txt for protein coverage summarizer
name <- paste(filename_s, "peptide_list.txt", sep = "_")
write.table(peptide, file = name, quote = F, sep = "\t", row.names = F, col.names = F)

# Remove flanking amino acids and terminal R,K and *'s
peptide_m <- gsub("\\#", "", peptide)
peptide_m <- gsub("\\*", "", peptide_m)
peptide_m <- substr(peptide_m,3,(nchar(peptide_m)-3))

# Loop through each peptide element by element
missed_cleavages <- data.frame()
peptide_length <- data.frame()
for (i in 1:length(peptide_m)) {
  pepSeq <- unlist(strsplit(peptide_m[i], NULL)) # Turn the string into a list of characters
  plength <- length(pepSeq) + 1
  peptide_length <- rbind(peptide_length,plength)
  missed <- 0
  for (j in 1:length(pepSeq)){ 
    if (pepSeq[j] == "K" & (pepSeq[j + 1] != "P" || is.na(pepSeq[j + 1]))) { 
      missed <- missed + 1
    } else if (pepSeq[j] == "R" & (pepSeq[j + 1] != "P" || is.na(pepSeq[j + 1]))) {
      missed <- missed + 1
    }} 
  missed_cleavages <- rbind(missed_cleavages, missed)} #Add rows to the data frame

# Calculate Kyte-Doolittle hydrophobicity score

# Enter Kyte-Doolittle Hydropathy Scores
I <- 4.5; V <- 4.2; L <- 3.8; F <- 2.8; C <- 2.5; M <- 1.9; A <- 1.8; G <- -0.4
T <- -0.7; W <- -0.9; S <- -0.8; Y <- -1.3; P <- -1.6; H <- -3.2; E <- -3.5
Q <- -3.5; D <- -3.5; N <- -3.5; K <- -3.9; R <- -4.5

# Remove flanking amino acids and #', *'s
peptide_h <- gsub("\\*", "", peptide)
peptide_h <- gsub("\\#", "", peptide_h)
peptide_h <- substr(peptide_h,3,(nchar(peptide_h)-2))

# Loop through the list and calculate the hydrophobicity for each peptide
hydrophobicity <- data.frame()
gravy <- data.frame()
for (i in 1:length(peptide_h)) {
  pepSeq <- unlist(strsplit(peptide_h[i],NULL))
  pepSeq_h <- unlist(mget(pepSeq, envir = as.environment(-1)))
  peptide_gravy <- sum(pepSeq_h)/length(pepSeq_h)
  peptide_gravy <- round(peptide_gravy, 2)
  peptide_hydrophobicity <- sum(pepSeq_h)
  hydrophobicity <- rbind(hydrophobicity, peptide_hydrophobicity)
  hydrophobicity <- round(hydrophobicity, 2)
  gravy <- rbind(gravy, peptide_gravy)
}

# Clear variables
rm(I,V,L,F,C,M,A,G,T,W,S,Y,P,H,E,Q,D,N,K,R)

# Put newly calculated variables into a single data frame to attach to xtab
peptide_characteristics <- cbind(peptide,peptide_length, missed_cleavages, gravy, hydrophobicity)
names(peptide_characteristics)[1:5] <- c("Peptide", "peptide_length", "missed_cleavages", 
                                         "GRAVY", "hydrophobicity"
)

massSpecOutput_trim <- merge(massSpecOutput_trim, peptide_characteristics, by = "Peptide")

# xtab <- sqldf("Select xtab.*, peptide_characteristics.*
#               From xtab
#               Join peptide_characteristics
#               ON xtab.Peptide = peptide_characteristics.Peptide
#               ")

# Subset the data by Job and calculate the summary statistics to create the metric table
metric_table <- data.frame(row.names = 
                             c("Job", "spectra_searched", "PSMs", "unique_peptides", "unique_proteins", "unique_proteins_2peptides",
                               "tryptic", "part_tryptic", "non_tryptic", "keratin_PSMs", "trypsin_PSMs", "albumin_PSMs", 
                               "percent_tryptic","peptides_w_missed", "missed_per_peptide", "missed_per_missed_peptide",
                               "avg_peptide_length", "median_peptide_length", "GRAVY", "hydrophobicity"
                             )
)

job = as.vector(unique(massSpecOutput_trim$Job))
for(i in 1:length(job)) {
  xtab_job <- subset(massSpecOutput_trim, Job == job[i])
  xtab_full_job <- subset(massSpecOutput, Job == job[i])
  spectra_searched <- length(xtab_full_job$Job)
  # Create a unique-peptide-level list for calculating sample metrics
  xtab_unique_peptide <- sqldf("Select Peptide, peptide_length, missed_cleavages, min(Scan) as Scan, 
                               MH, NTT, Protein, GRAVY, hydrophobicity
                               From xtab_job
                               Group by Peptide
                               ")
  
  # Create a protein level list for calculating sample metrics
  xtab_protein <- sqldf("Select Protein, count(Peptide) as Spectral_Counts, 
                        Count(DISTINCT Peptide) as Unique_Peptides
                        From xtab_job
                        Group by Protein
                        ")
  
  # Calculate sample metrics for peptide length and missed cleavages
  a_missed_per_peptide <- sum(xtab_unique_peptide$missed_cleavages / length(xtab_unique_peptide$peptide))
  a_missed_per_peptide <- round(a_missed_per_peptide, digits = 2)
  
  a_peptides_w_missed <- (length(xtab_unique_peptide$missed_cleavages[xtab_unique_peptide$missed_cleavages > 0])/length(xtab_unique_peptide$peptide))*100
  a_peptides_w_missed <- round(a_peptides_w_missed, digits = 2)
  
  a_avg_peptide_length <- mean(xtab_unique_peptide$peptide_length)
  a_avg_peptide_length <- round(a_avg_peptide_length, 2)
  
  a_median_peptide_length <- median(xtab_unique_peptide$peptide_length)
  
  x <- subset(xtab_unique_peptide, missed_cleavages > 0)
  a_missed_per_missed_peptide <- mean(x$missed_cleavages)
  a_missed_per_missed_peptide <- round(a_missed_per_missed_peptide, digits =2)
  rm(x)
  
  # Calculate hydrophobicity metrics
  a_mean_hydrophobicity <- mean(xtab_unique_peptide$hydrophobicity)
  a_mean_hydrophobicity <- round(a_mean_hydrophobicity, 2)
  a_mean_GRAVY <- mean(xtab_unique_peptide$GRAVY)
  a_mean_GRAVY <- round(a_mean_GRAVY, 2)
  
  # Identification metrics
  a_PSMs <- 1*length(xtab_job$peptide)
  a_unique_peptides <- 1*length(xtab_unique_peptide$peptide)
  a_unique_proteins <- 1*length(xtab_protein$Protein)
  
  # Proteins w/ 2 unique peptide IDs
  x <- subset(xtab_protein, Unique_Peptides > 1)
  a_unique_proteins_2peptides <- 1*length(x$Protein)
  rm(x)
  
  # Tryptic Peptides and percent tryptic
  x <- subset(xtab_unique_peptide, NTT == 2)
  a_tryptic <- 1*length(x$NTT)
  a_percent_tryptic <- round((a_tryptic/a_unique_peptides)*100, digits = 2)
  rm(x)
  
  # Partially Tryptic peptides
  x <- subset(xtab_unique_peptide, NTT == 1)
  a_part_tryptic <- 1*length(x$NTT)
  rm(x)
  
  # Non-Tryptic peptides
  x <- subset(xtab_unique_peptide, NTT == 0)
  a_non_tryptic <- 1*length(x$NTT)
  rm(x)
  
  # Count Contaminant PSMs
  x <- sqldf("Select Peptide
             From xtab_job
             Where Protein Like 'Contaminant_k%'
             ")
  a_keratin_PSMs <- 1*length(x$Peptide)
  rm(x)
  
  x <- sqldf("Select Peptide
             From xtab_job
             Where Protein Like 'Contaminant_Tryp%'
             ")
  a_trypsin_PSMs <- 1*length(x$Peptide)
  rm(x)
  
  x <- sqldf("Select Peptide
             From xtab_job
             Where Protein Like 'Contaminant_ALBU%'
             ")
  a_albumin_PSMs <- 1*length(x$Peptide)
  rm(x)
  
  # Create a table with all the sample metrics
  metric_table_temp <- rbind(job[i], spectra_searched, a_PSMs, a_unique_peptides, a_unique_proteins, 
                             a_unique_proteins_2peptides,a_tryptic, a_part_tryptic, a_non_tryptic, 
                             a_keratin_PSMs, a_trypsin_PSMs, a_albumin_PSMs, a_percent_tryptic,a_peptides_w_missed, 
                             a_missed_per_peptide, a_missed_per_missed_peptide,a_avg_peptide_length, 
                             a_median_peptide_length, a_mean_hydrophobicity, a_mean_GRAVY
  )
  
  metric_table <- cbind(metric_table, metric_table_temp)
}

#Write metric table to text
name <- paste(filename_s, "metric_table.txt", sep = "_")
write.table(metric_table, file = name, quote = F, sep = "\t", col.names = F,
            row.names = c("dataset", "spectra_searched", "PSMs", "unique_peptides", "unique_proteins", "unique_proteins_2peptides",
                          "tryptic", "part_tryptic", "non_tryptic", "keratin_PSMs", "trypsin_PSMs", 
                          "albumin_PSMs","percent_tryptic","peptides_w_missed", "missed_per_peptide", 
                          "missed_per_missed_peptide", "avg_peptide_length", 
                          "median_peptide_length", "mean_hydrophobicity", "mean_GRAVY")
)

# Make peptide-level spectral counting crosstab

for_xtab <- sqldf("Select Job, Protein, Peptide,
                  Count(Peptide) As Spectral_Counts
                  From massSpecOutput_trim
                  Group By  Job, Peptide"
)

spectral_count_xtab <- dcast(for_xtab, Peptide ~ Job, value.var = "Spectral_Counts")

protein_list <- sqldf("Select Protein, Peptide from for_xtab Group by Peptide")

spectral_count_xtab <- merge(protein_list, spectral_count_xtab, by = "Peptide", all.y = T)

rm(for_xtab)
write.table(spectral_count_xtab, "peptide.txt", sep = "\t", row.names = F, quote = F)

# Make protein-level spectral counting crosstab

for_xtab <- sqldf("Select Job, Protein,
                  Count(Peptide) As Spectral_Counts
                  From massSpecOutput_trim
                  Group By Job, Protein")

protein_xtab <- dcast(for_xtab, Protein ~ Job, value.var = "Spectral_Counts")

rm(for_xtab)
write.table(protein_xtab, "1pep-pro.txt", sep = "\t", row.names = F, quote = F)


#Make a crosstab with 2 peptide rule
for_xtab <- sqldf("Select count(DISTINCT PEPTIDE) as obs_cnt, Job, Protein,
                  Count(Peptide) As Spectral_Counts
                  From massSpecOutput_trim
                  Group By Job, Protein")

for_xtab <- subset(for_xtab, obs_cnt >= 2)
for_xtab <- for_xtab[,-1]

protein_xtab_2 <- dcast(for_xtab, Protein ~ Job, value.var = "Spectral_Counts")

rm(for_xtab)
write.table(protein_xtab_2, "2pep-pro.txt", sep = "\t", row.names = F, quote = F)
