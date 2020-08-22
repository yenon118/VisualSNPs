# Clear all
rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(lubridate)
library(ggplot2)

library(foreach)
library(doParallel)

library(argparse)



# Usage:
# Rscript genotyping.R \
# -c 10 \
# -fa /scratch/yenc/datasets/Iowa/Wm82.a2.v1.fa \
# -fq /scratch/yenc/datasets/Iowa/fastq \
# -of Iowa \
# -csw Chr \
# -cs 1 \
# -ce 20 \
# -o /scratch/yenc/datasets/Iowa/genotyping_output \
# --fastqc /scratch/yenc/tools/FastQC/fastqc \
# --bwa /scratch/yenc/tools/bwa/bwa \
# --samtools /scratch/yenc/tools/samtools/bin/samtools \
# --picard /scratch/yenc/tools/picard.jar \
# --gatk /scratch/yenc/tools/gatk/gatk-4.1.5.0/gatk \
# --vcf-to-tab /scratch/yenc/tools/vcftools/bin/vcf-to-tab \
# --snpeff /scratch/yenc/tools/snpEff/snpEff.jar \
# --snpeff-v Wm82.a2.v1 \
# --vcfeffoneperline /scratch/yenc/tools/snpEff/scripts/vcfEffOnePerLine.pl \
# --snpsift /scratch/yenc/tools/snpEff/SnpSift.jar


##########################################################################################################
# Parsing inputs
##########################################################################################################
# create parser object
parser = ArgumentParser()

parser$add_argument("-c", "--cores", help = "Computing cores", type="integer", default=1)
parser$add_argument("-fa", "--fasta-input", help = "Fasta input file path", required = TRUE)
parser$add_argument("-fq", "--fasta-input-folder", help = "Fastq input folder path", required = TRUE)
parser$add_argument("-of", "--output-filename", help = "Output filename", required = TRUE)
parser$add_argument("-csw", "--chromosome-start-with", help = "Chromosome start with")
parser$add_argument("-cew", "--chromosome-end-with", help = "Chromosome end with")
parser$add_argument("-cs", "--chromosome-start", help = "Chromosome start", type="integer", required = TRUE)
parser$add_argument("-ce", "--chromosome-end", help = "Chromosome end", type="integer", required = TRUE)
parser$add_argument("-o", "--output", help = "Output folder path", required = TRUE)

parser$add_argument("--fastqc", help = "FastQC", required = TRUE)
parser$add_argument("--bwa", help = "BWA", required = TRUE)
parser$add_argument("--samtools", help = "SAMtools", required = TRUE)
parser$add_argument("--picard", help = "Picard", required = TRUE)
parser$add_argument("--gatk", help = "GATK", required = TRUE)
parser$add_argument("--vcf-to-tab", help = "VCFtools vcf-to-tab", required = TRUE)
parser$add_argument("--snpeff", help = "SNPEff", required = TRUE)
parser$add_argument("--snpeff-v", help = "SNPEff reference data", required = TRUE)
parser$add_argument("--vcfeffoneperline", help = "vcfEffOnePerLine", required = TRUE)
parser$add_argument("--snpsift", help = "SNPSift", required = TRUE)

args = parser$parse_args()


##########################################################################################################
# input assignments
##########################################################################################################
cores = args$cores
fasta_input = args$fasta_input
fasta_input_folder = args$fasta_input_folder
output_filename = args$output_filename
chromosome_start_with = args$chromosome_start_with
chromosome_end_with = args$chromosome_end_with
chromosome_start = args$chromosome_start
chromosome_end = args$chromosome_end
output = args$output

fastqc = args$fastqc
bwa = args$bwa
samtools = args$samtools
picard = args$picard
gatk = args$gatk
vcf_to_tab = args$vcf_to_tab
snpeff = args$snpeff
snpeff_v = args$snpeff_v
vcfeffoneperline = args$vcfeffoneperline
snpsift = args$snpsift


##########################################################################################################
# Register workers for cluster
##########################################################################################################

registerDoParallel(cores = cores)

##########################################################################################################
# input checking section
##########################################################################################################
cat(rep("\n", 2))
print(cores)

if(file.exists(fasta_input)){
    print(fasta_input)
} else{
    print("fasta_input file does not exists!!!")
    quit(status = 1)
}
if(dir.exists(fasta_input_folder)){
    print(fasta_input_folder)
} else{
    print("fasta_input_folder does not exists!!!")
    quit(status = 1)
}

print(output_filename)
print(chromosome_start_with)
print(chromosome_end_with)
print(chromosome_start)
print(chromosome_end)

if(!dir.exists(output)){
    dir.create(output)
    if(dir.exists(output)){
        print(output)
    } else{
        print("output folder cannot be created!!!")
        quit(status = 1)
    }
}

if(file.exists(fastqc)){
    print(fastqc)
} else{
    print("fastqc does not exists!!!")
    quit(status = 1)
}
if(file.exists(bwa)){
    print(bwa)
} else{
    print("bwa does not exists!!!")
    quit(status = 1)
}
if(file.exists(samtools)){
    print(samtools)
} else{
    print("samtools does not exists!!!")
    quit(status = 1)
}
if(file.exists(picard)){
    print(picard)
} else{
    print("picard does not exists!!!")
    quit(status = 1)
}
if(file.exists(gatk)){
    print(gatk)
} else{
    print("gatk does not exists!!!")
    quit(status = 1)
}
if(file.exists(vcf_to_tab)){
    print(vcf_to_tab)
} else{
    print("vcf-to-tab does not exists!!!")
    quit(status = 1)
}
if(file.exists(snpeff)){
    print(snpeff)
} else{
    print("snpeff does not exists!!!")
    quit(status = 1)
}

print(snpeff_v)

if(file.exists(vcfeffoneperline)){
    print(vcfeffoneperline)
} else{
    print("vcfeffoneperline does not exists!!!")
    quit(status = 1)
}
if(file.exists(snpsift)){
    print(snpsift)
} else{
    print("snpsift does not exists!!!")
    quit(status = 1)
}

##########################################################################################################
# Output folder section
##########################################################################################################

folder_tmp = file.path(output, "tmp")
folder_index_fasta_log = file.path(output, "index_fasta_log")
folder_fastqc = file.path(output, "fastqc")
folder_fastq_coverage = file.path(output, "fastq_coverage")
folder_sam = file.path(output, "sam")
folder_sam_log = file.path(output, "sam_log")
folder_samtools_flagstat = file.path(output, "samtools_flagstat")
folder_samtools_stats = file.path(output, "samtools_stats")
folder_picard_sortsam = file.path(output, "picard_sortsam")
folder_picard_sortsam_log = file.path(output, "picard_sortsam_log")
folder_picard_markduplicates = file.path(output, "picard_markduplicates")
folder_picard_markduplicates_log = file.path(output, "picard_markduplicates_log")
folder_picard_markduplicates_metrics = file.path(output, "picard_markduplicates_metrics")
folder_picard_addorreplacereadgroups = file.path(output, "picard_addorreplacereadgroups")
folder_picard_addorreplacereadgroups_log = file.path(output, "picard_addorreplacereadgroups_log")
folder_gatk_haplotypecaller = file.path(output, "gatk_haplotypecaller")
folder_gatk_haplotypecaller_log = file.path(output, "gatk_haplotypecaller_log")
folder_InDel_matrix = file.path(output, "InDel_matrix")
folder_InDel_matrix_missing_percentage_plot = file.path(output, "InDel_matrix_missing_percentage_plot")
folder_InDel_positions = file.path(output, "InDel_positions")
folder_InDel_VCF = file.path(output, "InDel_VCF")
folder_InDel_VCF_EFF = file.path(output, "InDel_VCF_EFF")
folder_InDel_log = file.path(output, "InDel_log")
folder_SNP_matrix = file.path(output, "SNP_matrix")
folder_SNP_matrix_missing_percentage_plot = file.path(output, "SNP_matrix_missing_percentage_plot")
folder_SNP_positions = file.path(output, "SNP_positions")
folder_SNP_VCF = file.path(output, "SNP_VCF")
folder_SNP_VCF_EFF = file.path(output, "SNP_VCF_EFF")
folder_SNP_log = file.path(output, "SNP_log")

folders = c(
    folder_tmp,
    folder_index_fasta_log,
    folder_fastqc,
    folder_fastq_coverage,
    folder_sam,
    folder_sam_log,
    folder_samtools_flagstat,
    folder_samtools_stats,
    folder_picard_sortsam,
    folder_picard_sortsam_log,
    folder_picard_markduplicates,
    folder_picard_markduplicates_log,
    folder_picard_markduplicates_metrics,
    folder_picard_addorreplacereadgroups,
    folder_picard_addorreplacereadgroups_log,
    folder_gatk_haplotypecaller,
    folder_gatk_haplotypecaller_log,
    folder_InDel_matrix,
    folder_InDel_matrix_missing_percentage_plot,
    folder_InDel_positions,
    folder_InDel_VCF,
    folder_InDel_VCF_EFF,
    folder_InDel_log,
    folder_SNP_matrix,
    folder_SNP_matrix_missing_percentage_plot,
    folder_SNP_positions,
    folder_SNP_VCF,
    folder_SNP_VCF_EFF,
    folder_SNP_log
)

if(!dir.exists(output)){
    dir.create(output)
    if(!dir.exists(output)){
        print("Cannot create output directory!!!")
        quit(status = 1)
    }
}
for(i in 1:length(folders)){
    if(!dir.exists(folders[i])){
        dir.create(folders[i])
        if(!dir.exists(folders[i])){
            print(paste("Cannot create folder ", as.character(folders[i]), " !!!"))
            quit(status = 1)
        }
    }
}


##########################################################################################################
# Chromosome array strings
##########################################################################################################
chromosome_array = as.numeric(as.character(chromosome_start)):as.numeric(as.character(chromosome_end))
temp = nchar(chromosome_end)-nchar(chromosome_array)
chromosome_array_strings = ifelse(
    temp > 0,
    paste0(strrep("0",temp), chromosome_array),
    chromosome_array
)
if(!is.null(chromosome_start_with)){
    chromosome_array_strings = paste0(chromosome_start_with, chromosome_array_strings)
}
if(!is.null(chromosome_end_with)){
    chromosome_array_strings = paste0(chromosome_array_strings, chromosome_end_with)
}

cat(rep("\n", 2))
print(chromosome_array_strings)


##########################################################################################################
# Fastq input section
##########################################################################################################
files = list.files(fasta_input_folder)
sample_names = gsub("_.*", "", files)
fastq_number = as.numeric(gsub("[[:alpha:][:blank:][:punct:][:space:]]", "", gsub(".*_", "", files)))

fastq_inputs = data.frame(
    fastq_number = as.numeric(fastq_number),
    sample_names = sample_names,
    files = file.path(fasta_input_folder, files),
    stringsAsFactors = FALSE
)

fastq_inputs_1 = fastq_inputs[fastq_inputs$fastq_number == 1,]
fastq_inputs_1 = fastq_inputs_1[order(fastq_inputs_1$sample_names),]
fastq_inputs_2 = fastq_inputs[fastq_inputs$fastq_number == 2,]
fastq_inputs_2 = fastq_inputs_2[order(fastq_inputs_2$sample_names),]

cat(rep("\n", 2))
print(fastq_inputs_1)
print(fastq_inputs_2)

cat(rep("\n", 2))
if((nrow(fastq_inputs_1) != nrow(fastq_inputs_2)) & (identical(fastq_inputs_1$sample_names, fastq_inputs_2$sample_names))){
    print("Imbalance fastq_inputs!!!")
    quit(status = 1)
} else{
    write.table(
        x = fastq_inputs_1,
        file = file.path(output, "fastq_inputs_1.txt"),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    write.table(
        x = fastq_inputs_2,
        file = file.path(output, "fastq_inputs_2.txt"),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    write.table(
        x = fastq_inputs,
        file = file.path(output, "fastq_inputs.txt"),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
}


# ##########################################################################################################
# # Index fasta
# ##########################################################################################################
# index_fasta_commands = c(
#     paste(bwa, "index", fasta_input, sep = " "),
#     paste(samtools, "faidx", fasta_input, sep = " "),
#     paste(
#         "java", "-jar", picard, "CreateSequenceDictionary",
#         paste0("R=", fasta_input),
#         paste0("O=", sub("\\.fa$", ".dict", fasta_input)),
#         sep = " "
#     )
# )
#
# cat(rep("\n", 2))
# print(index_fasta_commands)
#
# cat(rep("\n", 2), "index fasta start!!!", rep("\n", 1))
# for(i in 1:length(index_fasta_commands)){
#     command = strsplit(index_fasta_commands[i], " ")[[1]]
#     print(file.path(folder_index_fasta_log, paste0("index_fasta_", i, ".log")))
#     system2(
#         command = command[1],
#         args = command[2:length(command)],
#         stdout = file.path(folder_index_fasta_log, paste0("index_fasta_", i, ".log")),
#         stderr = file.path(folder_index_fasta_log, paste0("index_fasta_", i, ".log"))
#     )
# }
# cat("index fasta complete!!!", rep("\n", 1))
#
#
# ##########################################################################################################
# # Check fastq quality with FastQC
# ##########################################################################################################
# fastqc_commands = paste(
#     fastqc, "-o", folder_fastqc, "-f", "fastq", fastq_inputs_1$files, fastq_inputs_2$files,
#     sep = " "
# )
#
# cat(rep("\n", 2))
# print(fastqc_commands)
#
# cat(rep("\n", 2), "FastQC start!!!", rep("\n", 1))
# results = foreach(i=1:length(fastqc_commands), .combine = rbind) %dopar% {
#     command = strsplit(fastqc_commands[i], " ")[[1]]
#     system2(
#         command = command[1],
#         args = command[2:length(command)]
#     )
# }
# cat("FastQC complete!!!", rep("\n", 1))
#
#
# ##########################################################################################################
# # Alignment with BWA
# ##########################################################################################################
# bwa_alignment_commands = paste(
#     bwa, "mem", "-t", "10",
#     "-M", fasta_input, fastq_inputs_1$files, fastq_inputs_2$files, ">",
#     file.path(folder_sam, paste0(fastq_inputs_1$sample_names, ".sam")),
#     sep = " "
# )
#
# cat(rep("\n", 2))
# print(bwa_alignment_commands)
#
# cat(rep("\n", 2), "bwa alignment start!!!", rep("\n", 1))
# results = foreach(i=1:length(bwa_alignment_commands), .combine = rbind) %dopar% {
#     command = strsplit(bwa_alignment_commands[i], " ")[[1]]
#     sample_name = gsub("\\..*", "", gsub(".*/", "", command[length(command)]))
#     print(file.path(folder_sam, paste0(sample_name, ".sam")))
#     print(file.path(folder_sam_log, paste0(sample_name, ".log")))
#     system2(
#         command = command[1],
#         args = command[2:length(command)],
#         stdout = file.path(folder_sam, paste0(sample_name, ".sam")),
#         stderr = file.path(folder_sam_log, paste0(sample_name, ".log"))
#     )
# }
# cat("bwa alignment complete!!!", rep("\n", 1))
#
#
# ##########################################################################################################
# # SAMtools flagstat
# ##########################################################################################################
# input_file_paths = file.path(folder_sam, list.files(folder_sam))
# output_filenames = gsub("\\..*", "", gsub(".*/", "", input_file_paths))
# output_file_paths = file.path(folder_samtools_flagstat, paste0(output_filenames, ".flagstat"))
# samtools_flagstat_commands = paste(
#     samtools, "flagstat", input_file_paths, ">", output_file_paths,
#     sep = " "
# )
#
# cat(rep("\n", 2))
# print(samtools_flagstat_commands)
#
# cat(rep("\n", 2), "samtools flagstat start!!!", rep("\n", 1))
# results = foreach(i=1:length(samtools_flagstat_commands), .combine = rbind) %dopar% {
#     command = strsplit(samtools_flagstat_commands[i], " ")[[1]]
#     log_file = gsub(".*>", "", command[length(command)])
#     print(log_file)
#     system2(
#         command = command[1],
#         args = command[2:length(command)],
#         stdout = log_file,
#         stderr = log_file
#     )
# }
# cat("samtools flagstat complete!!!", rep("\n", 1))
#
#
# ##########################################################################################################
# # SAMtools stats
# ##########################################################################################################
# input_file_paths = file.path(folder_sam, list.files(folder_sam))
# output_filenames = gsub("\\..*", "", gsub(".*/", "", input_file_paths))
# output_file_paths = file.path(folder_samtools_stats, paste0(output_filenames, ".stats"))
# samtools_stats_commands = paste(
#     samtools, "stats", input_file_paths, ">", output_file_paths,
#     sep = " "
# )
#
# cat(rep("\n", 2))
# print(samtools_stats_commands)
#
# cat(rep("\n", 2), "samtools stats start!!!", rep("\n", 1))
# results = foreach(i=1:length(samtools_stats_commands), .combine = rbind) %dopar% {
#     command = strsplit(samtools_stats_commands[i], " ")[[1]]
#     log_file = gsub(".*>", "", command[length(command)])
#     print(log_file)
#     system2(
#         command = command[1],
#         args = command[2:length(command)],
#         stdout = log_file,
#         stderr = log_file
#     )
# }
# cat("samtools stats complete!!!", rep("\n", 1))
#
#
# ##########################################################################################################
# # Sort SAM
# ##########################################################################################################
# input_file_paths = file.path(folder_sam, list.files(folder_sam))
# input_file_paths = input_file_paths[endsWith(input_file_paths, "sam")]
# output_filenames = gsub("\\..*", "", gsub(".*/", "", input_file_paths))
# output_file_paths = file.path(folder_picard_sortsam, paste0(output_filenames, ".bam"))
# sort_sam_commands = paste(
#     "java", "-jar", picard, "SortSam",
#     "CREATE_INDEX=TRUE", "MAX_RECORDS_IN_RAM=5000000", "SORT_ORDER=coordinate",
#     paste0("I=", input_file_paths), paste0("O=", output_file_paths),
#     sep = " "
# )
#
# cat(rep("\n", 2))
# print(sort_sam_commands)
#
# cat(rep("\n", 2), "sort sam start!!!", rep("\n", 1))
# results = foreach(i=1:length(sort_sam_commands), .combine = rbind) %dopar% {
#     command = strsplit(sort_sam_commands[i], " ")[[1]]
#     sample_name = gsub("\\..*", "", gsub(".*/", "", command[length(command)]))
#     print(file.path(folder_picard_sortsam_log, paste0(sample_name, ".log")))
#     system2(
#         command = command[1],
#         args = command[2:length(command)],
#         stdout = file.path(folder_picard_sortsam_log, paste0(sample_name, ".log")),
#         stderr = file.path(folder_picard_sortsam_log, paste0(sample_name, ".log"))
#     )
# }
# cat("sort sam complete!!!", rep("\n", 1))
#
#
# ##########################################################################################################
# # BAM mark duplicates
# ##########################################################################################################
# input_file_paths = file.path(folder_picard_sortsam, list.files(folder_picard_sortsam))
# input_file_paths = input_file_paths[endsWith(input_file_paths, "bam")]
# output_filenames = gsub("\\..*", "", gsub(".*/", "", input_file_paths))
# output_file_paths = file.path(folder_picard_markduplicates, paste0(output_filenames, "_mark_duplicates.bam"))
# output_metrics_file_paths = file.path(folder_picard_markduplicates_metrics, paste0(output_filenames, "_mark_duplicates.metrics"))
# mark_duplicated_commands = paste(
#     "java", "-jar", picard, "MarkDuplicates", "CREATE_INDEX=TRUE", "MAX_RECORDS_IN_RAM=5000000", "VALIDATION_STRINGENCY=LENIENT",
#     paste0("I=", input_file_paths),
#     paste0("METRICS_FILE=", output_metrics_file_paths),
#     paste0("O=", output_file_paths),
#     sep = " "
# )
#
# cat(rep("\n", 2))
# print(mark_duplicated_commands)
#
# cat(rep("\n", 2), "mark duplicates start!!!", rep("\n", 1))
# results = foreach(i=1:length(mark_duplicated_commands), .combine = rbind) %dopar% {
#     command = strsplit(mark_duplicated_commands[i], " ")[[1]]
#     sample_name = gsub("_mark_duplicates.*", "", gsub(".*/", "", command[length(command)]))
#     print(file.path(folder_picard_markduplicates_log, paste0(sample_name, "_mark_duplicates.log")))
#     system2(
#         command = command[1],
#         args = command[2:length(command)],
#         stdout = file.path(folder_picard_markduplicates_log, paste0(sample_name, "_mark_duplicates.log")),
#         stderr = file.path(folder_picard_markduplicates_log, paste0(sample_name, "_mark_duplicates.log"))
#     )
# }
# cat("mark duplicates complete!!!", rep("\n", 1))
#
#
# ##########################################################################################################
# # BAM add or replace reads group
# ##########################################################################################################
# input_file_paths = file.path(folder_picard_markduplicates, list.files(folder_picard_markduplicates))
# input_file_paths = input_file_paths[endsWith(input_file_paths, "bam")]
# output_filenames = gsub("_mark_duplicates.*", "", gsub(".*/", "", input_file_paths))
# output_file_paths = file.path(folder_picard_addorreplacereadgroups, paste0(output_filenames, "_add_replace.bam"))
# add_replace_commands = paste(
#     "java", "-jar", picard, "AddOrReplaceReadGroups",
#     "CREATE_INDEX=TRUE", "MAX_RECORDS_IN_RAM=5000000", "VALIDATION_STRINGENCY=LENIENT", "SORT_ORDER=coordinate",
#     paste0("RGID=", output_filenames), paste0("RGLB=", output_filenames), "RGPL=Illumina",
#     paste0("RGSM=", output_filenames), "RGCN=BGI", paste0("RGPU=", output_filenames),
#     paste0("I=", input_file_paths),
#     paste0("O=", output_file_paths),
#     sep = " "
# )
#
# cat(rep("\n", 2))
# print(add_replace_commands)
#
# cat(rep("\n", 2), "add or replace reads group start!!!", rep("\n", 1))
# results = foreach(i=1:length(add_replace_commands), .combine = rbind) %dopar% {
#     command = strsplit(add_replace_commands[i], " ")[[1]]
#     sample_name = gsub("_add_replace.*", "", gsub(".*/", "", command[length(command)]))
#     print(file.path(folder_picard_addorreplacereadgroups_log, paste0(sample_name, "_add_replace.log")))
#     system2(
#         command = command[1],
#         args = command[2:length(command)],
#         stdout = file.path(folder_picard_addorreplacereadgroups_log, paste0(sample_name, "_add_replace.log")),
#         stderr = file.path(folder_picard_addorreplacereadgroups_log, paste0(sample_name, "_add_replace.log"))
#     )
# }
# cat("add or replace reads group complete!!!", rep("\n", 1))


##########################################################################################################
# Haplotype caller
##########################################################################################################
input_file_paths = file.path(folder_picard_addorreplacereadgroups, list.files(folder_picard_addorreplacereadgroups))
input_file_paths = input_file_paths[endsWith(input_file_paths, "bam")]
output_filenames = gsub("_add_replace.*", "", gsub(".*/", "", input_file_paths))
output_file_paths = file.path(folder_gatk_haplotypecaller, paste0(output_filenames, ".g.vcf"))
haplotype_caller_commands = paste(
    gatk, "HaplotypeCaller", "--tmp-dir", folder_tmp,
    "--emit-ref-confidence", "GVCF",
    "--annotation-group", "AS_StandardAnnotation",
    "-R", fasta_input, "-I", input_file_paths, "-O", output_file_paths,
    sep = " "
)

cat(rep("\n", 2))
print(haplotype_caller_commands)

cat(rep("\n", 2), "haplotype caller start!!!", rep("\n", 1))
results = foreach(i=1:length(haplotype_caller_commands), .combine = rbind) %dopar% {
    command = strsplit(haplotype_caller_commands[i], " ")[[1]]
    sample_name = gsub("\\.g\\.vcf.*", "", gsub(".*/", "", command[length(command)]))
    print(file.path(folder_gatk_haplotypecaller_log, paste0(sample_name, ".log")))
    system2(
        command = command[1],
        args = command[2:length(command)],
        stdout = file.path(folder_gatk_haplotypecaller_log, paste0(sample_name, ".log")),
        stderr = file.path(folder_gatk_haplotypecaller_log, paste0(sample_name, ".log"))
    )
}
cat("haplotype caller complete!!!", rep("\n", 1))


##########################################################################################################
# Combine GVCF
##########################################################################################################
input_file_paths = file.path(folder_gatk_haplotypecaller, list.files(folder_gatk_haplotypecaller))
input_file_paths = input_file_paths[endsWith(input_file_paths, ".g.vcf")]
output_file_paths = file.path(output, paste0(output_filename, ".g.vcf"))
combine_gvcf_command = paste(
    gatk, "CombineGVCFs", "--tmp-dir", folder_tmp,
    "-R", fasta_input,
    paste("-V", input_file_paths, sep = " ", collapse = " "),
    "-O", output_file_paths,
    sep = " "
)

cat(rep("\n", 2))
print(combine_gvcf_command)

cat(rep("\n", 2), "combine gvcf start!!!", rep("\n", 1))
command = strsplit(combine_gvcf_command, " ")[[1]]
print(file.path(output, paste0(output_filename, "_combine_gvcf.log")))
system2(
    command = command[1],
    args = command[2:length(command)],
    stdout = file.path(output, paste0(output_filename, "_combine_gvcf.log")),
    stderr = file.path(output, paste0(output_filename, "_combine_gvcf.log"))
)
cat("combine gvcf complete!!!", rep("\n", 1))



##########################################################################################################
# Stop cluster
##########################################################################################################
cat(rep("\n", 2))
stopImplicitCluster()
