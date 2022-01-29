library(tercen)
library(dplyr, warn.conflicts = FALSE)

options("tercen.workflowId" = "b09f25d899d082c8e044aa0e7204d7d5")
options("tercen.stepId"     = "de830731-4bd3-404b-baee-c3cad331b2d2")

getOption("tercen.workflowId")
getOption("tercen.stepId")

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

py_install("Levenshtein")

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

ctx %>%
  select(.y, .ci, .ri) %>%
  group_by(.ci, .ri) %>%
  summarise(mean = mean(.y)) %>%
  ctx$addNamespace() %>%
  ctx$save()