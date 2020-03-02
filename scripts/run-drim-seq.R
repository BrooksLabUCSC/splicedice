#!/usr/bin/env Rscript

library(DRIMSeq)
library(ggplot2)
library(optparse)

parse_manifest <- function(filename) {
    table <- read.table(filename)
    colnames(table) <- c("sample_id", "path", "type", "group")
    return(table)
}

parser <- OptionParser()
parser <- add_option(parser, c("-m", "--manifest"))
parser <- add_option(parser, c("-d", "--drim-table"))
parser <- add_option(parser, c("-o", "--output"))
parser <- add_option(parser, c("-c", "--ctrl"))
parser <- add_option(parser, c("-t", "--threads"), type="integer", default=1)

args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
manifest_filename <- args$manifest
ctrl <- args$ctrl
drim_table_filename <- args$drim_table
output_filename <- args$output
threads <- args$threads

# Make my plots pretty
options(repr.plot.width = 6, repr.plot.height = 4)

# read in the output from my conversion script.
counts <- read.delim(drim_table_filename, header = TRUE, sep = "\t")
names(counts)[names(counts) == "gene"] <- "gene_id"

head(counts)

manifest <- parse_manifest(manifest_filename)
manifest

# Load the object and check that it loaded the right number of

# "GENES" and sampes.
d <- dmDSdata(counts = counts, samples = manifest)
print(d)

head(counts(d), 3)
transcripts_plot <- plotData(d)
ggsave("transcripts_per_gene.png", plot = transcripts_plot)

design_full <- model.matrix(~ group, data = samples(d))
design_full

set.seed(123)

d <- dmPrecision(d, design = design_full,
                 BPPARAM = BiocParallel::MulticoreParam(threads))

d
head(mean_expression(d), 3)
common_precision(d)
head(genewise_precision(d))

# ggplot2 plotting
precision_plot <- plotPrecision(d) + geom_point(size = 4)
ggsave("precision_plot.png", plot = precision_plot)


d <- dmFit(d, design = design_full, verbose = 1,
           BPPARAM = BiocParallel::MulticoreParam(threads))

head(proportions(d))

head(coefficients(d), level = "feature")

coef_str <- paste("group", ctrl, sep="")
d <- dmTest(d, coef = coef_str, verbose = 1,
            BPPARAM = BiocParallel::MulticoreParam(threads))

head(results(d), 3)

res <- results(d)

# TODO: maybe order by adj_pvalue would be better
res <- res[order(res$pvalue, decreasing = FALSE), ]
top_gene_id <- res$gene_id[1]
head(res, 20)

# TODO: additional plots check DRIMSeq documentation
pvalue_plot <- plotPValues(d)
ggsave("pvalue_plot_by_gene.png", plot = pvalue_plot)

write.table(file = output_filename, x = res, sep = "\t", quote = FALSE,
            row.names = FALSE)
