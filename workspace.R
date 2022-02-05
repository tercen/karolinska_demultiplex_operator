library(tercen)
library(dplyr, warn.conflicts = FALSE)
library(base64enc)
library(stringr)

options("tercen.workflowId" = "b09f25d899d082c8e044aa0e7204d7d5")
options("tercen.stepId"     = "de830731-4bd3-404b-baee-c3cad331b2d2")

getOption("tercen.workflowId")
getOption("tercen.stepId")


serialize.to.string = function(object){
  con = rawConnection(raw(0), "r+")
  saveRDS(object, con)
  str64 = base64enc::base64encode(rawConnectionValue(con))
  close(con)
  return(str64)
}

ctx = tercenCtx()

documentIds <- ctx$cselect()


for (id in documentIds[[1]]) {
  
  res <- try(ctx$client$fileService$get(id),silent = TRUE)
  if (class(res) == "try-error") stop("Supplied column values are not valid documentIds.")
  
  
}

if (length(documentIds) != 1) stop("Should input only one documentId.")

filename <- ctx$client$fileService$get(documentIds[[1]])$name

writeBin(ctx$client$fileService$download(documentIds[[1]]), filename)

system(paste0("unzip ", filename))

system("python3 -m pip install Levenshtein")

r1_file <- list.files(".", "_R1_", recursive = TRUE,
                      full.names = TRUE)

r2_file <- list.files(".", "_R2_", recursive = TRUE,
                      full.names = TRUE)

col_file <- list.files(".", "col.txt", recursive = TRUE,
                       full.names = TRUE)

row_file <- list.files(".", "row.txt", recursive = TRUE,
                       full.names = TRUE)

system("mkdir demultiplexed_files/")
system(paste("python3 demultiplex_TCR_fastqs_by_row_and_column_barcodes.py",
             r1_file, r2_file, "demultiplexed_files/ --gzip_output yes --row_barcodes_file",
             row_file, "--col_barcodes_file", col_file))

output_r1_files <- list.files("demultiplexed_files",
                              "_R1.fastq",
                              full.names = TRUE)

output_table <- c()

for (sample_R1 in output_r1_files) {
  
  
  bytes_R1 <- readBin(file(sample_R1, 'rb'),
                      raw(),
                      n=file.info(sample_R1)$size)
  
  sample_R2 <- str_replace(sample_R1, "_R1.fastq", "_R2.fastq")
  
  bytes_R2 <- readBin(file(sample_R2, 'rb'),
                      raw(),
                      n=file.info(sample_R2)$size)
  
  string_val1 <- serialize.to.string(bytes_R1)
  string_val2 <- serialize.to.string(bytes_R2)
  
  sample_name <- str_split(basename(sample_R1),
                           "_R1.fastq")[[1]][[1]]
  
  output_table <- bind_rows(output_table,
                            tibble(sample = sample_name,
                                   .forward_read_fastq_data = string_val1,
                                   .reverse_read_fastq_data = string_val2))
  
}

output_table %>%
  mutate(.ci = 0) %>%
  ctx$addNamespace() %>%
  ctx$save()
