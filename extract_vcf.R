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
# Rscript extract_vcf.R \
# -c 10 \
# -si /scratch/yenc/datasets/Nebraska_Yen/Nebraska_snp.vcf \
# -ii /scratch/yenc/datasets/Nebraska_Yen/Nebraska_indel.vcf \
# -of Nebraska \
# -csw Chr \
# -cs 1 \
# -ce 20 \
# -o /scratch/yenc/datasets/Nebraska_Yen/temp \
# --vcf-to-tab /scratch/yenc/tools/vcftools/bin/vcf-to-tab \
# --snpeff /scratch/yenc/tools/snpEff/snpEff.jar \
# --snpeff-v Wm82.a2.v1 \
# --vcfeffoneperline /scratch/yenc/tools/snpEff/scripts/vcfEffOnePerLine.pl \
# --snpsift /scratch/yenc/tools/snpEff/SnpSift.jar


# create parser object
parser = ArgumentParser()

parser$add_argument("-c", "--cores", help = "Computing cores", type="integer", default=1)
parser$add_argument("-si", "--snp-input", help = "SNP input file path", required = TRUE)
parser$add_argument("-ii", "--indel-input", help = "Indel input file path", required = TRUE)
parser$add_argument("-of", "--output-filename", help = "Output filename", required = TRUE)
parser$add_argument("-csw", "--chromosome-start-with", help = "Chromosome start with")
parser$add_argument("-cew", "--chromosome-end-with", help = "Chromosome end with")
parser$add_argument("-cs", "--chromosome-start", help = "Chromosome start", type="integer", required = TRUE)
parser$add_argument("-ce", "--chromosome-end", help = "Chromosome end", type="integer", required = TRUE)
parser$add_argument("-o", "--output", help = "Output folder path", required = TRUE)
parser$add_argument("--vcf-to-tab", help = "VCFtools vcf-to-tab", required = TRUE)
parser$add_argument("--snpeff", help = "SNPEff", required = TRUE)
parser$add_argument("--snpeff-v", help = "SNPEff reference data", required = TRUE)
parser$add_argument("--vcfeffoneperline", help = "vcfEffOnePerLine", required = TRUE)
parser$add_argument("--snpsift", help = "SNPSift", required = TRUE)

args = parser$parse_args()


cores = args$cores
snp_input = args$snp_input
indel_input = args$indel_input
output_filename = args$output_filename
chromosome_start_with = args$chromosome_start_with
chromosome_end_with = args$chromosome_end_with
chromosome_start = args$chromosome_start
chromosome_end = args$chromosome_end
output = args$output
vcf_to_tab = args$vcf_to_tab
snpeff = args$snpeff
snpeff_v = args$snpeff_v
vcfeffoneperline = args$vcfeffoneperline
snpsift = args$snpsift


cat(rep("\n", 2))
print(cores)
print(snp_input)
print(indel_input)
print(output_filename)
print(chromosome_start_with)
print(chromosome_end_with)
print(chromosome_start)
print(chromosome_end)
print(output)
print(vcf_to_tab)
print(snpeff)
print(snpeff_v)
print(vcfeffoneperline)
print(snpsift)


registerDoParallel(cores = cores)


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


chromosome_array = as.numeric(as.character(chromosome_start)):as.numeric(as.character(chromosome_end))
temp = nchar(chromosome_end)-nchar(chromosome_array)
chromosome_array_strings = ifelse(
    temp > 0,
    paste0(strrep("0",temp), chromosome_array),
    chromosome_array
)
if(!is.null(chromosome_start_with)){
    chromosome_array_strings = paste0(chromosome_start_with, chromosome_array_strings)
} else if(!is.null(chromosome_end_with)){
    chromosome_array_strings = paste0(chromosome_array_strings, chromosome_end_with)
} else{
    print("Chromosome related input error!!!")
    quit(status = 1)
}

cat(rep("\n", 2))
print(chromosome_array_strings)


grep_snp_commands = paste0(
    "grep -e \"^#\" -e \"^",
    chromosome_array_strings,
    "\" ",
    snp_input,
    " > ",
    file.path(folder_SNP_VCF, paste0(output_filename, "_", chromosome_array_strings, ".vcf"))
)

grep_indel_commands = paste0(
    "grep -e \"^#\" -e \"^",
    chromosome_array_strings,
    "\" ",
    indel_input,
    " > ",
    file.path(folder_InDel_VCF, paste0(output_filename, "_", chromosome_array_strings, ".vcf"))
)

vcf_to_tab_snp_commands = paste0(
    vcf_to_tab,
    " < ",
    file.path(folder_SNP_VCF, paste0(output_filename, "_", chromosome_array_strings, ".vcf")),
    " > ",
    file.path(folder_SNP_matrix, paste0(output_filename, "_", chromosome_array_strings, ".txt"))
)

vcf_to_tab_indel_commands = paste0(
    vcf_to_tab,
    " < ",
    file.path(folder_InDel_VCF, paste0(output_filename, "_", chromosome_array_strings, ".vcf")),
    " > ",
    file.path(folder_InDel_matrix, paste0(output_filename, "_", chromosome_array_strings, ".txt"))
)

awk_snp_commands = paste0(
    "awk -F \"\\t\" '{print $1, $2}' FS='\\t' OFS='\\t' ",
    file.path(folder_SNP_matrix, paste0(output_filename, "_", chromosome_array_strings, ".txt")),
    " > ",
    file.path(folder_SNP_positions, paste0(output_filename, "_chr_pos_", chromosome_array_strings, ".txt"))
)

awk_indel_commands = paste0(
    "awk -F \"\\t\" '{print $1, $2}' FS='\\t' OFS='\\t' ",
    file.path(folder_InDel_matrix, paste0(output_filename, "_", chromosome_array_strings, ".txt")),
    " > ",
    file.path(folder_InDel_positions, paste0(output_filename, "_chr_pos_", chromosome_array_strings, ".txt"))
)

snpeff_snp_commands = paste0(
    "java -jar ",
    snpeff,
    " -s ",
    file.path(folder_SNP_VCF_EFF, paste0(output_filename, "_", chromosome_array_strings, ".html")),
    " -v ",
    snpeff_v,
    " ",
    file.path(folder_SNP_VCF, paste0(output_filename, "_", chromosome_array_strings, ".vcf")),
    " > ",
    file.path(folder_SNP_VCF_EFF, paste0(output_filename, "_", chromosome_array_strings, ".eff.vcf"))
)

snpeff_indel_commands = paste0(
    "java -jar ",
    snpeff,
    " -s ",
    file.path(folder_InDel_VCF_EFF, paste0(output_filename, "_", chromosome_array_strings, ".html")),
    " -v ",
    snpeff_v,
    " ",
    file.path(folder_InDel_VCF, paste0(output_filename, "_", chromosome_array_strings, ".vcf")),
    " > ",
    file.path(folder_InDel_VCF_EFF, paste0(output_filename, "_", chromosome_array_strings, ".eff.vcf"))
)

vcfeffoneperline_snp_commands = paste0(
    "cat ",
    file.path(folder_SNP_VCF_EFF, paste0(output_filename, "_", chromosome_array_strings, ".eff.vcf")),
    " | ", vcfeffoneperline, " | java -jar ", snpsift,
    " extractFields - \"CHROM\" \"POS\" \"REF\" \"EFF[*].EFFECT\" \"EFF[*].IMPACT\" \"EFF[*].GENE\" \"EFF[*].FUNCLASS\" \"EFF[*].CODON\" \"EFF[*].AA\" \"EFF[*].AA_LEN\"  \"EFF[*].BIOTYPE\" \"EFF[*].CODING\" \"EFF[*].TRID\" \"EFF[*].RANK\" ",
    " > ",
    file.path(folder_SNP_VCF_EFF, paste0(output_filename, "_", chromosome_array_strings, ".eff.txt"))
)

vcfeffoneperline_indel_commands = paste0(
    "cat ",
    file.path(folder_InDel_VCF_EFF, paste0(output_filename, "_", chromosome_array_strings, ".eff.vcf")),
    " | ", vcfeffoneperline, " | java -jar ", snpsift,
    " extractFields - \"CHROM\" \"POS\" \"REF\" \"EFF[*].EFFECT\" \"EFF[*].IMPACT\" \"EFF[*].GENE\" \"EFF[*].FUNCLASS\" \"EFF[*].CODON\" \"EFF[*].AA\" \"EFF[*].AA_LEN\"  \"EFF[*].BIOTYPE\" \"EFF[*].CODING\" \"EFF[*].TRID\" \"EFF[*].RANK\" ",
    " > ",
    file.path(folder_InDel_VCF_EFF, paste0(output_filename, "_", chromosome_array_strings, ".eff.txt"))
)


cat(rep("\n", 2))
cat(grep_snp_commands, sep = "\n")
cat(grep_indel_commands, sep = "\n")

cat(rep("\n", 2))
cat(vcf_to_tab_snp_commands, sep = "\n")
cat(vcf_to_tab_indel_commands, sep = "\n")

cat(rep("\n", 2))
cat(awk_snp_commands, sep = "\n")
cat(awk_indel_commands, sep = "\n")

cat(rep("\n", 2))
cat(snpeff_snp_commands, sep = "\n")
cat(snpeff_indel_commands, sep = "\n")

cat(rep("\n", 2))
cat(vcfeffoneperline_snp_commands, sep = "\n")
cat(vcfeffoneperline_indel_commands, sep = "\n")


cat(rep("\n", 2), "grep start!!!", rep("\n", 2))
results = foreach(i=1:length(grep_snp_commands), .combine = c) %dopar% {
    log = system(grep_snp_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_SNP_log, "SNP_grep.log"))

results = foreach(i=1:length(grep_indel_commands), .combine = c) %dopar% {
    log = system(grep_indel_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_InDel_log, "InDel_grep.log"))
cat(rep("\n", 2), "grep complete!!!", rep("\n", 2))


cat(rep("\n", 2), "vcf-to-tab start!!!", rep("\n", 2))
results = foreach(i=1:length(vcf_to_tab_snp_commands), .combine = c) %dopar% {
    log = system(vcf_to_tab_snp_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_SNP_log, "SNP_vcf_to_tab.log"))

results = foreach(i=1:length(vcf_to_tab_indel_commands), .combine = c) %dopar% {
    log = system(vcf_to_tab_indel_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_InDel_log, "InDel_vcf_to_tab.log"))
cat(rep("\n", 2), "vcf-to-tab complete!!!", rep("\n", 2))


cat(rep("\n", 2), "awk start!!!", rep("\n", 2))
results = foreach(i=1:length(awk_snp_commands), .combine = c) %dopar% {
    log = system(awk_snp_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_SNP_log, "SNP_awk.log"))

results = foreach(i=1:length(awk_indel_commands), .combine = c) %dopar% {
    log = system(awk_indel_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_InDel_log, "InDel_awk.log"))
cat(rep("\n", 2), "awk complete!!!", rep("\n", 2))


cat(rep("\n", 2), "snpeff start!!!", rep("\n", 2))
results = foreach(i=1:length(snpeff_snp_commands), .combine = c) %dopar% {
    log = system(snpeff_snp_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_SNP_log, "SNP_snpeff.log"))

results = foreach(i=1:length(snpeff_indel_commands), .combine = c) %dopar% {
    log = system(snpeff_indel_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_InDel_log, "InDel_snpeff.log"))
cat(rep("\n", 2), "snpeff complete!!!", rep("\n", 2))


cat(rep("\n", 2), "vcfeffoneperline start!!!", rep("\n", 2))
results = foreach(i=1:length(vcfeffoneperline_snp_commands), .combine = c) %dopar% {
    log = system(vcfeffoneperline_snp_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_SNP_log, "SNP_vcfeffoneperline.log"))

results = foreach(i=1:length(vcfeffoneperline_indel_commands), .combine = c) %dopar% {
    log = system(vcfeffoneperline_indel_commands[i], intern = TRUE)
    log = c(log, "", "")
    return(log)
}
writeLines(text = results, con = file.path(folder_InDel_log, "InDel_vcfeffoneperline.log"))
cat(rep("\n", 2), "vcfeffoneperline complete!!!", rep("\n", 2))





stopImplicitCluster()
