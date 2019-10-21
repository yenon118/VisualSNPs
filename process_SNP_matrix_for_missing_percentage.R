# Clear all
rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

library(argparse)

library(foreach)
library(doParallel)



# create parser object
parser = ArgumentParser()

parser$add_argument("-c", "--cores", help = "Number of processing cores", type = "integer", default = 1)

parser$add_argument("-i", "--input", help = "Input folder path", type = "character", required = TRUE)
parser$add_argument("-ifsw", "--input-file-starts-with", help = "Input file start with", type = "character", default = "Chr")
parser$add_argument("-ife", "--input-file-extension", help = "Input file extension", type = "character", required = TRUE)

parser$add_argument("-o", "--output", help = "Output folder path", type = "character", required = TRUE)

args = parser$parse_args()


cores = args$cores

input_folder_path = file.path(as.character(args$input))
input_file_starts_with = args$input_file_starts_with
input_file_extension = args$input_file_extension

output_folder_path = file.path(as.character(args$output))


cat(rep("\n", 2))
print(cores)
print(input_folder_path)
print(input_file_starts_with)
print(input_file_extension)
print(output_folder_path)


if(!dir.exists(input_folder_path)){
    cat("Input folder does not exist!!!")
    quit(save = "no", status = -1)
} else{
    input_folder_path = normalizePath(input_folder_path)
}
if(!dir.exists(input_folder_path)){
    cat("Input folder does not exist!!!")
    quit(save = "no", status = -1)
}


if(!dir.exists(output_folder_path)){
    dir.create(output_folder_path)
}
if(!dir.exists(output_folder_path)){
    cat("Output directory cannot be created!!!")
    quit(save = "no", status = -1)
} else{
    output_folder_path = normalizePath(output_folder_path)
}
if(!dir.exists(output_folder_path)){
    cat("Output directory cannot be created!!!")
    quit(save = "no", status = -1)
}


filenames = list.files(input_folder_path)

if(length(filenames) == 0){
    cat("There is no file in that directory!!!")
    quit(save = "no", status = -1)
}



filenames = filenames[startsWith(filenames, input_file_starts_with) & endsWith(filenames, input_file_extension)]

file_paths = file.path(input_folder_path, filenames)
chromosomes = as.integer(gsub("\\D", "", gsub(input_file_starts_with, "", gsub(input_folder_path, "", file_paths))))

chromosomes_dat_filepaths = data.frame(Chromosome = chromosomes, Filepaths = file_paths, stringsAsFactors = FALSE)

chromosomes_dat_filepaths = chromosomes_dat_filepaths[order(as.numeric(chromosomes_dat_filepaths[,1])),]
cat(rep("\n", 2));print(chromosomes_dat_filepaths)


registerDoParallel(cores = cores)

result = foreach(i=1:nrow(chromosomes_dat_filepaths), .combine = rbind) %dopar% {

    temp_index = match(i, chromosomes_dat_filepaths[,1])
    dat = read.table(file = file.path(as.character(chromosomes_dat_filepaths[temp_index, 2])), header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")

    colnames(dat)[colnames(dat) == "#CHROM"] = "CHROM"

    missing = as.integer(apply(dat[,4:ncol(dat)], 1, function(x) { return(length(x[x == "./."])) }))

    return(
        data.frame(
            CHROM = dat[,1],
            POS = dat[,2],
            REF = dat[,3],
            Missing = missing,
            Total = ncol(dat)-3,
            stringsAsFactors = FALSE
        )
    )
}

result$Missing_percentage = as.integer(100 * (result$Missing / result$Total))

result = result %>% arrange(desc(Missing_percentage))

write.table(x = result, file = file.path(output_folder_path, "missing_percentage_per_chromosome.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

missing_percentage_per_chromosome = result %>%
                                        group_by(CHROM, Missing_percentage) %>%
                                        summarize(Count = n()) %>%
                                        arrange(CHROM, desc(Missing_percentage))

write.table(x = missing_percentage_per_chromosome, file = file.path(output_folder_path, "missing_percentage_per_chromosome_total.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")



result = foreach(i=1:nrow(chromosomes_dat_filepaths), .combine = rbind) %dopar% {

    temp_index = match(i, chromosomes_dat_filepaths[,1])
    dat = read.table(file = file.path(as.character(chromosomes_dat_filepaths[temp_index, 2])), header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")

    colnames(dat)[colnames(dat) == "#CHROM"] = "CHROM"

    missing = as.integer(apply(dat[,4:ncol(dat)], 2, function(x) { return(length(x[x == "./."])) }))

    return(
        data.frame(
            CHROM = dat[1,1],
            LINE = colnames(dat)[4:ncol(dat)],
            Missing = missing,
            Total = nrow(dat),
            stringsAsFactors = FALSE
        )
    )
}

result$Missing_percentage = as.numeric(100 * (result$Missing / result$Total))

result = result %>% arrange(CHROM, desc(Missing_percentage), LINE)

write.table(x = result, file = file.path(output_folder_path, "missing_percentage_per_line.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

missing_percentage_per_line = result %>%
                                group_by(LINE) %>%
                                summarize(Missing = sum(Missing, na.rm = TRUE), Total = sum(Total, na.rm = TRUE)) %>%
                                mutate(Missing_percentage = 100 * (Missing / Total)) %>%
                                arrange(desc(Missing_percentage))

write.table(x = missing_percentage_per_line, file = file.path(output_folder_path, "missing_percentage_per_line_total.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

stopImplicitCluster()



p <- ggplot(data = missing_percentage_per_chromosome, mapping = aes(x = factor(Missing_percentage), y = Count, fill = factor(CHROM))) +
        geom_bar(stat="identity", color="black") +
        labs(title = "Number of Missing Percentage for Positions in Each Chromosome",
                y = "Count",
                x = "Missing Percentage \n (for each x in integer range [0,100): >= x and < x+1)") +
        scale_fill_discrete(name = "Chromosome") +
        theme_classic() +
        theme(
            plot.title = element_text(size = 40, hjust = 0.5, face = "bold"),
            axis.title = element_text(size = 28),
            axis.text.y = element_text(size = 24),
            axis.text.x = element_text(size = 16),
            legend.title = element_text(size = 24),
            legend.text = element_text(size = 24)
        )

ggsave(filename = file.path("missing_percentage_per_chromosome_total.png"), plot = p, path = output_folder_path, width = 32, height = 18, dpi = 800)



