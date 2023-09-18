#Prepare TCGA data for survival analysis 

############################################################################

#Load packages and data
library(tidyverse)
library(rlang)

#Load METABRIC tables, and select the important columns

##Patients
patients <- read.delim('/Users/IEO5505/Desktop/ML_TDA/Survival_analysis/TCGA_and_METABRIC/data/METABRIC.BRCA/data_clinical_patient.txt')
pat_to_retain <- c("PATIENT_ID", "SEX", "AGE_AT_DIAGNOSIS", "OS_STATUS", "OS_MONTHS",
                   "RFS_STATUS", "RFS_MONTHS")
patients <- patients[, pat_to_retain] %>% na.omit()

##Samples
samples <- read.delim('/Users/IEO5505/Desktop/ML_TDA/Survival_analysis/TCGA_and_METABRIC/data/METABRIC.BRCA/data_clinical_sample.txt')
sam_to_retain <- c('PATIENT_ID', 'SAMPLE_ID', 'ER_STATUS', 'HER2_STATUS', 'PR_STATUS', 'TUMOR_STAGE')
samples <- samples[, sam_to_retain] %>% na.omit()

##Matrix
matrix <- read.delim('/Users/IEO5505/Desktop/ML_TDA/Survival_analysis/TCGA_and_METABRIC/data/METABRIC.BRCA/data_expression_median.txt')

###Reformat colnames (Samples_ID)
s <- colnames(matrix)[3:dim(matrix)[2]]
s <- sapply( 1:length(s), function(i) { x <- paste0('MB-', paste(str_split(s[i], '')[[1]][4:7], collapse = '')) } )
colnames(matrix)[3:dim(matrix)[2]] <- s
###Reformat gene names and transpose matrix
matrix <- matrix %>% select(-c(Entrez_Gene_Id)) %>% na.omit()
index <- c()
unique_genes <- c()
for (i in 1:length(matrix$Hugo_Symbol)) {
    if (!matrix$Hugo_Symbol[i] %in% unique_genes) {
        unique_genes <- c(unique_genes, matrix$Hugo_Symbol[i])
        index <- c(index, i)
    }
}
matrix <- matrix[index, ]
row.names(matrix) <- matrix$Hugo_Symbol
matrix <- matrix %>% select(-c(Hugo_Symbol)) 
matrix <- as.data.frame(t(matrix))


##


#Filter interesting obs and merge: patients and samples
indeces_to_retain <- intersect(intersect(patients$PATIENT_ID, samples$SAMPLE_ID), row.names(matrix))
clinical <- merge(patients[patients$PATIENT_ID %in% indeces_to_retain, ],
                  samples[samples$PATIENT_ID %in% indeces_to_retain, ], by = 'PATIENT_ID')

#Remove stage iv
clinical <- clinical %>% filter(TUMOR_STAGE != 4)
clinical <- clinical %>% select(-c(SEX))

#Create subtype column
clinical <- clinical %>% mutate(subtype = case_when(
        (ER_STATUS == 'Negative') & (HER2_STATUS == 'Negative') & (PR_STATUS == 'Negative') ~ 'TNBC',
        (ER_STATUS == 'Negative') & (HER2_STATUS == 'Positive') & (PR_STATUS == 'Negative') ~ 'Her2',
        (ER_STATUS == 'Positive') & (HER2_STATUS == 'Negative') & (PR_STATUS == 'Positive') ~ 'Luminal'))
clinical$subtype[is.na(clinical$subtype)] <- 'unknown'

#Reorder and rename columns
clinical <- clinical %>% select(SAMPLE_ID, PATIENT_ID, AGE_AT_DIAGNOSIS, TUMOR_STAGE, subtype, OS_STATUS, 
                    OS_MONTHS, RFS_STATUS, RFS_MONTHS)
colnames(clinical) <- c('sample', 'case', 'age_at_diagnosis', 'stage', 
                        'subtype', 'OS_status', 'OS_time', 'Prog_free_status', 
                        'Prog_free_time')                     

#Change OS and PFS indicators, convert times in days
clinical$OS_status[clinical$OS_status == '0:LIVING'] <- FALSE
clinical$OS_status[clinical$OS_status == '1:DECEASED'] <- TRUE
clinical$Prog_free_status[clinical$Prog_free_status == '0:Not Recurred'] <- FALSE
clinical$Prog_free_status[clinical$Prog_free_status == '1:Recurred'] <- TRUE
clinical$OS_time <- clinical$OS_time * 31
clinical$Prog_free_time <- clinical$Prog_free_time * 31


##


#Merge with matrix
matrix$sample <- row.names(matrix)
row.names(matrix) <- NULL
df <- merge(clinical, matrix, by = 'sample')

#Save entire TCGA processed data
columns <- colnames(df)
saveRDS(df, '/Users/IEO5505/Desktop/ML_TDA/Survival_analysis/TCGA_and_METABRIC/data/new/complete_METABRIC.rds')
#write.table(df, '/Users/IEO5505/Desktop/ML_TDA/TCGA_data/TCGA_data/new/TCGA_processed/complete_df.csv', sep = ',')