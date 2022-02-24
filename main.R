library(tercen)
library(dplyr, warn.conflicts = FALSE)
library(stringr)

ctx <- tercenCtx()

# Define input and output paths
input_path <- "/var/lib/tercen/share/read/files_to_demultiplex"

output_path <- "/var/lib/tercen/share/write/demultiplexed_fastqs/"

system(paste("mkdir", output_path)

r1_file <- list.files(input_path, "_R1_", recursive = TRUE,
                      full.names = TRUE)

r2_file <- list.files(input_path, "_R2_", recursive = TRUE,
                      full.names = TRUE)

col_file <- list.files(input_path, "col.txt", recursive = TRUE,
                       full.names = TRUE)

row_file <- list.files(input_path, "row.txt", recursive = TRUE,
                       full.names = TRUE)


print(paste0("python3 demultiplex_TCR_fastqs_by_row_and_column_barcodes.py ",
              r1_file, " ", r2_file, " ", output_path, " --gzip_output yes --row_barcodes_file ",
              row_file, " --col_barcodes_file ", col_file))

system(paste0("python3 demultiplex_TCR_fastqs_by_row_and_column_barcodes.py ",
              r1_file, " ", r2_file, " ", output_path, " --gzip_output yes --row_barcodes_file ",
              row_file, " --col_barcodes_file ", col_file))

output_r1_files <- list.files(output_path,
                              "_R1.fastq",
                              full.names = TRUE)

output_table <- c()

for (sample_R1 in output_r1_files) {
  
  sample_name <- str_split(basename(sample_R1),
                           "_R1.fastq")[[1]][[1]]
  
  number_of_lines <- as.integer(system(paste("zcat",
                                             sample_R1,
                                             "| wc -l | awk '{print $1}'"),
                                       intern = TRUE))
  
  number_of_reads <- number_of_lines / 4
  
  output_table <- bind_rows(output_table,
                            tibble(sample = sample_name,
                                   read_number = number_of_reads))
  
}

output_table %>%
  mutate(.ci = 0) %>%
  ctx$addNamespace() %>%
  ctx$save()
