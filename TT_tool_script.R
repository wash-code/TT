# Library -----

library(Rsamtools)
library(stringr)
library(stringi)
library(tidyverse)
library(optparse)


# Parsing ------
option_list <- list(
  make_option(c("-i", "--input1"), type = "character",
              help = "Path to input FASTQ file 1", metavar = "character"),
  make_option(c("-j", "--input2"), type = "character",
              help = "Path to input FASTQ file 2", metavar = "character"),
  make_option(c("-b", "--build_file"), type = "character", 
              help = "Path to build file directory", metavar = "character"),
  make_option(c("-g", "--genome_idx"), type = "character",
              help = "Path to genome index", metavar = "character"),
  
  make_option(c("-s", "--sample"), type = "character", default = "test_file",
              help = "Sample name [default = %default]", metavar = "character"),
  make_option(c("-w", "--write_dir"), type = "character", 
              help = "Path to output directory", metavar = "character"),
  make_option(c("-q", "--query_type"), type = "logical", default = FALSE,
              help = "This should be TRUE if the read contains R. Some sequencing platform provides the read with R! Mostly it will not and needs to be added absed on the flag!  [default = %default]", metavar = "integer"),
  make_option(c("-t", "--threads"), type = "integer", default = 4,
              help = "Number of threads [default = %default]", metavar = "integer"),
  make_option(c("-m", "--map_diff"), type = "integer", default = 4,
              help = "Number of matches by which circ alignment has to supercede the genome alignment [default = %default]", metavar = "integer"),
  
  make_option(c("-k", "--keep_temp"), type = "logical", default = FALSE,
              help = "Keep temporary files (TRUE/FALSE) [default = %default]", metavar = "logical"),
  make_option(c("-c", "--conda_path"), type = "character", default = "/Users/sb/miniforge3/etc/profile.d/conda.sh",
              help = "Path to conda.sh [default = %default]", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

ip1 <- opt$input1
ip2 <- opt$input2
build_file <- opt$build_file
genome_idx <- opt$genome_idx
write_dir <- opt$write_dir
threads <- opt$threads
keep_temp <- opt$keep_temp
samp <- opt$sample
query_type <- opt$query_type
map_diff <- opt$map_diff 
conda_path <- opt$conda_path




rev_comp <- function(nucSeq) { return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq))) }
add_mates <- function(read_names) {
  mate1 <- gsub("/2$", "/1", read_names)
  mate2 <- gsub("/1$", "/2", read_names)
  unique(c(read_names, mate1, mate2))
}

count_m_from_cigar <- function(cigar) {
  if (is.na(cigar)) {
    return(0)
  }
  sum(as.numeric(unlist(regmatches(cigar, gregexpr("\\d+(?=M)", cigar, perl = TRUE)))))
}



# Set WD -----------

# Path to input ---

avoid_cigars <- c( 'H', 'P') # c('N', 'H', 'P', '=', 'X')
avoid_pattern <- paste(avoid_cigars, collapse = "|")
tool <- "HISAT2"

# Path to output ---

index_dir <- file.path(build_file , "Index_Creation" , tool )
fasta_file <- paste0(build_file, "junct_fasta.fa")
bsj_fasta <- file.path (build_file , "bsj_strand_seq_fasta.fa" ) 



output_dir <- file.path(write_dir , samp)
dir.create(output_dir)


# Step 2: Align RNA-seq Reads to circRNA Junction Index ----------

denovo_bam <- file.path(output_dir, "denovo_bsj_alignments.sorted.bam")

hisat2_cmd <- sprintf( # Build HISAT2 command for alignment to circRNA index
  "hisat2 -p %d -x %s -1 %s -2 %s | samtools view -@ %d -bS - | samtools sort -@ %d -o %s && samtools index %s",
  threads,
  file.path(index_dir, "junct_hisat2_index"),
  ip1,
  ip2,
  threads,
  threads,
  denovo_bam,
  denovo_bam
)



full_cmd <- sprintf(
  "source %s && conda activate env_create && %s",
  conda_path,
  #"conda run -n env_create && %s", # conda activate
  
  hisat2_cmd
)


cat("Command to run:\n", full_cmd, "\n")

system(full_cmd)


# Checking True circRNA ---------

result_list <- list()

bsj_dict <- readDNAStringSet( bsj_fasta )
bsj_dict <- setNames( as.character(bsj_dict) , names(bsj_dict) )
bsj_df <- data.frame(
  ReferenceName = names(bsj_dict),
  BSJ = as.character(bsj_dict),
  stringsAsFactors = FALSE
)

bsj_df$RC_BSJ = sapply(bsj_df$BSJ, rev_comp)



bam <- BamFile(denovo_bam)
param <- ScanBamParam(what = c("qname", "rname", "seq", "mapq" , "cigar" , "flag" , "mrnm") , 
                      flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE) )

bam_data <- scanBam(bam, param = param)[[1]]
circ_result_df <- as.data.frame(bam_data)


colnames(circ_result_df) <- c("ReadName", "Flag", "ReferenceName", "MAPQ", "cigar", "mate_info", "Sequence")


if (!query_type) {
  circ_result_df <- circ_result_df %>%
    mutate(read_pair = ReadName ) %>%
    mutate(
      ReadName = case_when(
        bitwAnd(Flag, 64) == 0 ~ paste0(ReadName, ".1"),
        bitwAnd(Flag, 128) == 0 ~ paste0(ReadName, ".2"), 
        TRUE ~ ReadName)
    )
} else {
  circ_result_df <- circ_result_df %>% mutate(read_pair = substr(ReadName, 1, nchar(ReadName) - 2))
}


circ_result_df <- circ_result_df %>% left_join(bsj_df, by = "ReferenceName")




circ_result_df$contain_bsj <- str_detect(circ_result_df$Sequence, fixed(circ_result_df$BSJ)) |
  str_detect(circ_result_df$Sequence, fixed(circ_result_df$RC_BSJ))




selected_read_pairs <- circ_result_df %>%
  filter(contain_bsj == 1) %>%
  pull(read_pair) %>% unique()  # This may not work still for the query_name = TRUE check once



circ_result_df <- circ_result_df %>%
  mutate(Pairwise_Retain = ifelse(read_pair %in% selected_read_pairs, 1, 0 )) %>%  # Pairwise_Retain helps retain any read and its mate if they contain the BSJ becuase selected read_pairs are the ones selected due to the contain_bsj functionality
  mutate(Pair_Match = ifelse(mate_info == ReferenceName , 1, 0 ) )


circ_result_df <- circ_result_df %>% filter(Pairwise_Retain == 1) 
circ_result_df <- circ_result_df %>% mutate( matches = sapply(cigar, count_m_from_cigar) )
colnames(circ_result_df) <- paste0("circular_" , colnames(circ_result_df))



# Sub-setting the FASTQ files --------

read_names <- circ_result_df$circular_read_pair %>% unique() # This may not work still for the query_name = TRUE check once
#complete_read_names <- add_mates(read_names)
subset_read_names <- file.path(output_dir , "selected_read_names.txt" ) 
writeLines(read_names, subset_read_names  )



output_fq1 <- file.path(output_dir , "subset_R1.fq" )
output_fq2 <- file.path(output_dir , "subset_R2.fq" )



full_cmd <- sprintf(
  "source %s && conda activate env_create && %s",
  conda_path,
  #"conda run -n env_create && %s", # conda activate env_create
  paste("seqtk subseq", ip1, subset_read_names, ">", output_fq1)
)
system(full_cmd)



full_cmd <- sprintf(
  "source %s && conda activate env_create && %s",
  conda_path,
  #" conda run -n env_create && %s", #conda activate env_create 
  paste("seqtk subseq", ip2, subset_read_names, ">", output_fq2)
)
system(full_cmd)


# Genome align the fastq Files Genome based ------


hisat_bam <- file.path(output_dir, "genome_aligned_reads.bam")



hisat2_genome_cmd <- sprintf( # Build HISAT2 command for alignment to linear index
  "hisat2 -p %d -x %s -1 %s -2 %s | samtools view -@ %d -bS - | samtools sort -@ %d -o %s && samtools index %s",
  threads,
  genome_idx,
  output_fq1,
  output_fq2,
  threads,
  threads,
  hisat_bam,
  hisat_bam
)



full_cmd <- sprintf(
  "source %s && conda activate env_create && %s",
  conda_path,
  #"conda run -n env_create && %s", #conda activate env_create 
  hisat2_genome_cmd
)


cat("Command to run:\n", full_cmd, "\n")

system(full_cmd)



# Filter more FPs



bam <- BamFile(hisat_bam)
param <- ScanBamParam(what =  c("qname", "rname", "seq", "mapq" , "cigar" , "flag" , "mrnm") , 
                      flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE))
bam_data <- scanBam(bam, param = param)[[1]]




genome_result_df <- as.data.frame(bam_data)


colnames(genome_result_df) <- c("ReadName", "Flag", "ReferenceName", "MAPQ", "cigar", "mate_info", "Sequence")


if (!query_type) {
  genome_result_df <- genome_result_df %>%
    mutate(read_pair = ReadName ) %>%
    mutate(
      ReadName = case_when(
        bitwAnd(Flag, 64) == 0 ~ paste0(ReadName, ".1"),
        bitwAnd(Flag, 128) == 0 ~ paste0(ReadName, ".2"), 
        TRUE ~ ReadName)
    )
} else {
  genome_result_df <- genome_result_df %>% mutate(read_pair = substr(ReadName, 1, nchar(ReadName) - 2))
}

genome_result_df <- genome_result_df %>% mutate( matches = sapply(cigar, count_m_from_cigar) )
colnames(genome_result_df) <- paste0("genome_" , colnames(genome_result_df))

combined_result_df <- left_join(
  circ_result_df, 
  genome_result_df, 
  by = c("circular_ReadName" = "genome_ReadName") ) %>% mutate( genome_matches = ifelse(is.na(genome_matches), 0, genome_matches) )



combined_result_df <- combined_result_df %>%
  mutate(contains_avoid_cigar = ifelse(is.na(genome_cigar), TRUE, grepl(avoid_pattern, genome_cigar))) %>%
  mutate(circ_match_better =  ifelse( circular_matches > genome_matches + map_diff , TRUE, FALSE ) )

combined_result_df$True_Cand <- combined_result_df$contains_avoid_cigar | combined_result_df$circ_match_better


true_circ <- combined_result_df  %>% filter(True_Cand == 1) %>% filter(circular_Pair_Match == 1)

count_values <- true_circ %>%
  group_by(circular_ReferenceName) %>%
  summarise(Count = n(), .groups = 'drop')


save(count_values , file = file.path (output_dir ,  "Results.RData") )
save(combined_result_df , true_circ , file = file.path (output_dir ,  "All_candiates.RData") )
write.table(count_values, file = file.path (output_dir , "count_values.txt" ), sep = "\t", row.names = FALSE, col.names = FALSE , quote = FALSE)


if (!keep_temp) {
  unlink(file.path(output_dir, "*.bam*"), recursive = FALSE)
  unlink(file.path(output_dir, "*.fq"), recursive = FALSE)
}
