# Load libraries
library(MSstats)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gplots)

# Set prefix for naming output files
prefix <- "20210217-PA3-ab"

# Set L2FC and pval cutoffs
l2fc_max <- 1
l2fc_min <- -1
pval_max <- 0.05
adjpval_max <- 0.05

# Set database file names
dbfiles <- c("20191010.SwissProt.Hsapiens.fasta", "Uniprot.HIVNL43.fasta")

# Get the number of databases indicated in the params file
numdb <- length(dbfiles)

# Loop through each database file
for (i in 1:numdb) {
  # If it's the first database then read lines into dblines
  if (i == 1) {
    dblines <- readLines(con = dbfiles[i])
  } else {
    # If it's not the first database then read lines into tmp
    # Then concatenate tmp onto the end of dblines
    tmp <- readLines(con = dbfiles[i])
    dblines <- c(dblines, tmp)
  }
}

# Grab header lines starting with > character
headerlines <- dblines[grep(">", dblines)]

# Extract protein accession and protein name from header lines 
protein.accession <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+.*)")[,2]
protein.description <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+.*)")[,3]
protein.name <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+)_.*")[,3]

# Create a named vector to lookup protein descriptions by accession
getProteinDescription <- protein.description
names(getProteinDescription) <- protein.accession
getProteinName <- protein.name
names(getProteinName) <- protein.accession

# Set working directory
setwd("C:/Jeff/20210216-PA3-ab/")

# Read in Spectronaut report
raw <- read.delim(file = "20210216_201224_EC20210127_Report.tsv", quote = "")

# Convert to MSstats format
quant <- SpectronauttoMSstatsFormat(raw)

# Normalize by equalize medians
processed <- dataProcess(quant, normalization = "equalizeMedians", summaryMethod = "TMP", cutoffCensored = "minFeature", censoredInt = "0", MBimpute = TRUE, maxQuantileforCensored = 0.999)

# Read in comparisons file
comparisons <- read.delim(file = "20210216-PA3-ab-comparisons.txt")

# Grab condition columns from comparisons
comparisonmatrix <- comparisons[,2:ncol(comparisons)]

# Reorder conditions in alphabetical order (MSstast will fail if you don't do this)
comparisonmatrix <- as.matrix(comparisonmatrix[,order(names(comparisonmatrix))])

# Add names of comparisons as row names 
row.names(comparisonmatrix) <- comparisons[,1]

print("Running group comparisons.")
compared <- groupComparison(contrast.matrix = comparisonmatrix, data = processed)

# Get subquant from processed
subquant <- quantification(processed)

# Store MSstats results (log2fold-changes and adjusted p-values) in results data frame
results <- compared$ComparisonResult

# Split accession column by semi-colons
split_acc <- str_split(string = results$Protein, pattern = ";")

# Go through all accessions, lookup descriptions, and paste to a new column in evidence named Protein.Description
print("Annotating protein descriptions.")
for (i in 1:length(split_acc)) {
  results[i, "Protein.Description"] <- paste(unname(getProteinDescription[split_acc[[i]]]), collapse = ";")
  results[i, "Protein.Name"] <- paste(unname(getProteinName[split_acc[[i]]]), collapse = ";")
}

# Write results file to text
print("Writing results file.")
write.table(x = results, file = paste(prefix, "-results-ann.txt", sep = '',  collapse = ''), sep = "\t", quote = FALSE, row.names = FALSE)

# Melt and cast data to wide format
results_melted <- melt(data = results, id.vars = c("Protein", "Protein.Description", "Label"), measure.vars = c("log2FC", "pvalue", "adj.pvalue"))
results_wide <- dcast(data = results_melted, formula = Protein + Protein.Description ~ Label + variable, value.var = "value")
write.table(x = results_wide, file = paste(prefix, "-results-wide.txt", collapse = '', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)

# Calculate LogP values for volcano plots
results$LogP <- -1 * log(results$pvalue) / log(10)

# Create VolcanoChange column in results data frame to label up/down regulated 
results$VolcanoChange = "NO CHANGE"
results$VolcanoChange[which((results$log2FC > l2fc_max) & (results$adj.pvalue < adjpval_max) & (results$pvalue < pval_max))] = "UP"
results$VolcanoChange[which((results$log2FC < l2fc_min) & (results$adj.pvalue < adjpval_max) & (results$pvalue < pval_max))] = "DOWN"

# Set color palette to red (down) gray (no change) and blue (up)
cols <- c("UP" = "blue", "DOWN" = "red", "NO CHANGE" = "gray")

# Plot volcano plots and save to pdf
pdf(paste(prefix, "-volcano-plots.pdf", collapse = '', sep = ''))
ggplot(data=results, aes(x=results$log2FC, y=results$LogP)) + geom_point(size=0.25, aes(color = factor(results$VolcanoChange))) + facet_wrap(Label~.) + scale_color_manual(values = cols) + theme_bw() + xlab(label = "Log2Fold-Change") + ylab(label = "-Log10(adjusted p-value)") + theme(legend.position = "none") + coord_equal(ratio=length(levels(results$Label)))
dev.off()

# Extract differentially increased/decreased data
# Must be real numbers (i.e., not Inf, -Inf, or NA)
upregulated <- results[which((results$log2FC > l2fc_max) & (results$pvalue < pval_max) & (results$adj.pvalue < adjpval_max) & (is.finite(results$log2FC)) & (is.finite(results$pvalue))),]
downregulated <- results[which((results$log2FC < l2fc_min) & (results$pvalue < pval_max) & (results$adj.pvalue < adjpval_max) & (is.finite(results$log2FC)) & (is.finite(results$pvalue))),]

# Count number of up/down regulated sites using table function
upregulatedcounts <- as.data.frame(table(upregulated$Label))
downregulatedcounts <- as.data.frame(table(downregulated$Label))

# Label change direction so up and down regulated can be combined in a single data frame
upregulatedcounts$Direction <- "Increased"
downregulatedcounts$Direction <- "Decreased"
significantcounts <- rbind(upregulatedcounts, downregulatedcounts)

# Set colors for bar plots
colors <- c("Increased" = "blue", "Decreased" = "red")

# Plot num diff abundant bar plots and save to pdf
pdf(paste(prefix, "-diff-abun-barplots.pdf", sep = '', collapse = ''))
ggplot(data=significantcounts, aes(x=significantcounts$Direction, y=significantcounts$Freq)) + geom_bar(stat = "identity", aes(color = factor(significantcounts$Direction), fill = factor(significantcounts$Direction))) + facet_wrap(significantcounts$Var1 ~ .) + xlab(label = NULL) + ylab(label = "Num differentially abundant") + theme_bw() + scale_color_manual(values = colors) + scale_fill_manual(values = colors) + theme(legend.position = "none", aspect.ratio = 1) + geom_text(aes(label = significantcounts$Freq),  vjust = -0.5) + scale_y_continuous(expand = c(0,0), limits = c(0, 1.1 * max(significantcounts$Freq)))
dev.off()

# Extract numerical data from subquant matrix for sample correlation heatmap
subquantmatrix <- as.matrix(subquant[,2:ncol(subquant)])

# Calculate pairwise correlations for all samples
cormat <- round(cor(subquantmatrix, use="complete.obs"),2)

# Row names are the same as column names
row.names(cormat) <- colnames(cormat)

# Set color scale from white to blue from 0 to 1 in 0.1 increments
breaks <- seq(0, 1, 0.1)
colorpalette <- colorRampPalette(c("white", "blue"))

# Plot heatmap of correlations and save to pdf
pdf(paste(prefix, "-sample-correlation-heatmap.pdf", collapse = '', sep = ''))
heatmap.2(x=cormat, col = colorpalette(100), trace="none", density.info = "none", margins = c(10,10))
dev.off()

# Extract numerical data from subquant matrix for sample correlation heatmap
subquantmatrix <- as.matrix(subquant[,2:ncol(subquant)])

# Calculate principle components based on subquant matrix
pca <- prcomp(na.omit(subquantmatrix), center = TRUE, scale = TRUE)

# Extract PC1 and PC2 for plotting
pcaplot <- as.data.frame(cbind(pca$rotation[,1], pca$rotation[,2]))

# Plot PC1 vs PC2 and save to pdf
pdf(paste(prefix, "-pca-plot.pdf", collapse = '', sep = ''))
ggplot(data=pcaplot, aes(x=pcaplot[,1], y=pcaplot[,2])) + geom_point(aes(color = row.names(pcaplot))) + xlab(label = "Principle Component 1") + ylab(label = "Principle Component 2") + theme_bw() + theme(legend.title = element_blank(), aspect.ratio = 1)
dev.off()

# Create function to count number of proteins detected per sample
# Counting non-NA values in each column of subquant
proteincountfunc <- function (x) length(which(!(is.na(subquant[,x]))))

# Apply function across each column of subquant and store as a data frame
proteincounts <- data.frame(unlist(lapply(X = 2:ncol(subquant), FUN = proteincountfunc)))

# Rename column of data frame
colnames(proteincounts) <- "Num. Proteins"

# Add column of data frame to label Sample names
proteincounts$Sample <- colnames(subquant[,2:ncol(subquant)])

# Plot protein counts as a bar graph
pdf(paste(prefix, "-protein-counts.pdf", collapse = '', sep = ''))
ggplot(data = proteincounts, aes(x = proteincounts$Sample, y = proteincounts$`Num. Proteins`)) + geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(expand = c(0,0), limits = c(0, 1.1 * max(proteincounts$`Num. Proteins`))) + geom_text(aes(label = proteincounts$`Num. Proteins`, vjust = -0.5)) + xlab(label = NULL) + ylab(label = "Num proteins detected")
dev.off()

# Melt subquant data to plot intensity box plots by sample
subquant_melted <- melt(data = subquant, id.vars = "Protein", measure.vars = colnames(subquant[,2:ncol(subquant)]), variable.name = "Sample", value.name = "Log2Intensity")

# Create a function to calculate the median of subquant log2intensity for each sample
subquantmedianfunc <- function (x) round(median(subquant[,x], na.rm = TRUE), 2)

# Apply function across each column of data in subquant
subquantmedians <- data.frame(unlist(lapply(X = 2:ncol(subquant), FUN = subquantmedianfunc)))

# Rename column header
colnames(subquantmedians) <- "Median Log2Intensity"

subquantmedians$Sample <- colnames(subquant[,2:ncol(subquant)])

# Plot
pdf(paste(prefix, "-sample-intensities.pdf", collapse = '', sep = ''))
ggplot(data = subquant_melted, aes(x = subquant_melted$Sample, y = subquant_melted$Log2Intensity)) + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + xlab(label = NULL) + ylab(label = "Log2Intensity") + geom_text(data = subquantmedians, aes(x = subquantmedians$Sample, y = subquantmedians$`Median Log2Intensity`, label = subquantmedians$`Median Log2Intensity`), vjust = -0.5, size = 2)
dev.off()

