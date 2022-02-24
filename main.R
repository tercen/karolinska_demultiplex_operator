library(tercen)
library(dplyr, warn.conflicts = FALSE)
library(stringr)

ctx <- tercenCtx()

# Define input path
input_path <- "/var/lib/tercen/share/read/files_to_demultiplex"

# Check if a "files_to_demultiplex" folder exists and is not empty
if( dir.exists(input_path) == FALSE) {

  stop("ERROR: files_to_demultiplex folder does not exist in project read folder.")

}

if (length(dir(input_path)) == 0) {
  stop("ERROR: files_to_demultiplex folder is empty.")
}

# Define and create output paths

output_path <- "/var/lib/tercen/share/write/demultiplexed_fastqs/"

system(paste("mkdir", output_path))

# Check if individual files are present in the input folder

r1_file <- list.files(input_path, "_R1_", recursive = TRUE,
                      full.names = TRUE)

if (length(r1_file) == 0) stop("ERROR: could not find a forward read file.")

r2_file <- list.files(input_path, "_R2_", recursive = TRUE,
                      full.names = TRUE)

if (length(r2_file) == 0) stop("ERROR: could not find a reverse read file.")

col_file <- list.files(input_path, "col.txt", recursive = TRUE,
                       full.names = TRUE)

if (length(col_file) == 0) stop("ERROR: could not find a col.txt file.")

row_file <- list.files(input_path, "row.txt", recursive = TRUE,
                       full.names = TRUE)

if (length(row_file) == 0) stop("ERROR: could not find a row.txt file.")


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
