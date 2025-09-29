setwd("C:/Users/darin/Desktop/AI_Omics_Internship_2025/Module_I")
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","limma","arrayQualityMetrics",
                       "AnnotationDbi","hgu133plus2.db"))

install.packages("dplyr")

library("GEOquery")
library("affy")
library("limma")
library("arrayQualityMetrics")
library("AnnotationDbi")
library("hgu133plus2.db")
library("dplyr")

gse_data <- getGEO("GSE38642", GSEMatrix = TRUE)

expression_data <- exprs(gse_data$GSE38642_series_matrix.txt.gz)

feature_data <- fData(gse_data$GSE38642_series_matrix.txt.gz)

phenotype_data <- pData(gse_data$GSE38642_series_matrix.txt.gz)

sum(is.na(phenotype_data$source_name_ch1))
 
untar("GSE38642/GSE38642_RAW.tar", exdir = "GSE38642/CEL_Files")

raw_data <- ReadAffy(celfile.path ="GSE38642/CEL_Files")
raw_data


arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)

normalized_data <- rma(raw_data)

arrayQualityMetrics(expressionset = normalized_data[, c(3,20,25,27,45,62)],
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)
processed_data <- as.data.frame(exprs(normalized_data))
dim(processed_data)

row_median <- rowMedians(as.matrix(processed_data))
row_median


hist(row_median, breaks = 100, freq = FALSE, main = "Median Intensity Distribution")

threshold <- 4
abline(v = threshold, col ="black", lwd = 2)
indx <- row_median > threshold
filtered_data <- processed_data[indx, ]

colnames(filtered_data) <- rownames(phenotype_data)
processed_data <- filtered_data

class(phenotype_data$characteristics_ch1.5)
phenotype_data$characteristics_ch1.5
groups <- factor(phenotype_data$characteristics_ch1.5,
                 levels = c("status: non-diabetic donors", "status: T2D donors"),
                 label = c("normal", "diabetic"))
class(groups)
levels(groups)




