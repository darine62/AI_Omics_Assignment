setwd("C:/Users/darin/Desktop/AI_Omics_Internship_2025/Module_I")
input_dir <- "raw_data"
output_dir <- "results"

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

result_list <- list()

classify_gene <- function(logFC, padj) {
  if(!is.na(padj) && padj < 0.05 && logFC > 1) {
    return("Upregulated")
  }
  else if(!is.na(padj) && padj < 0.05 && logFC < -1) {
    return("Downregulated")
  }
  else {
    return("Not_significant")
  }
}

for(file_name in files_to_process) {
  cat("\nProcesing:", file_name, "\n")
  
  input_file_path <- file.path(input_dir, file_name)
  data <- read.csv(input_file_path, header = TRUE)
  
  data$padj[is.na(data$padj)] <- 1
  
  data$status <- mapply(classify_gene, data$padj, data$logFC)
  
  output_file_path <- file.path(output_dir, paste0("Processed_", file_name))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
  
  cat("Summary for", file_name, ":\n")
  print(table(data$status))
}
