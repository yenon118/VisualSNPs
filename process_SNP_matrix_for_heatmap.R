# Clear all
rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)

library(argparse)



# create parser object
parser = ArgumentParser()

parser$add_argument("-i", "--input", help = "Input file path", required = TRUE)

parser$add_argument("-o", "--output", help = "Output file path", required = TRUE)

args = parser$parse_args()


input_filepath = file.path(as.character(args$input))
output_filepath = file.path(as.character(args$output))


# cat(rep("\n", 2))
# print(input_filepath)
# print(output_filepath)
# print(dirname(output_filepath))
# print(basename(output_filepath))


if(file.exists(input_filepath)){
    input_filepath = normalizePath(input_filepath)
    if(file.exists(input_filepath)){
        dat = read.table(file = input_filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
    } else{
        cat("Input file does not exist!!!")
        quit(save = "no", status = -1)
    }
} else{
    cat("Input file does not exist!!!")
    quit(save = "no", status = -1)
}


if(!dir.exists(dirname(output_filepath))){
    dir.create(dirname(output_filepath))
}
if(!dir.exists(dirname(output_filepath))){
    cat("Output directory cannot be created!!!")
    quit(save = "no", status = -1)
} else{
    output_filepath = file.path(normalizePath(dirname(output_filepath)), basename(output_filepath))
}
if(!dir.exists(dirname(output_filepath))){
    cat("Output directory cannot be created!!!")
    quit(save = "no", status = -1)
}



colnames(dat)[colnames(dat) == "#CHROM"] = "CHROM"

colnames(dat)[4:ncol(dat)] = gsub("\\.variant.*", "", colnames(dat)[4:ncol(dat)])

results = gather(dat, Line, SNP, -CHROM, -POS, -REF) %>% select(CHROM, POS, REF, SNP, Line)



write.table(x = results, file = file.path(output_filepath), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")


