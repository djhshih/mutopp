# Preprocess COSMIC signatures

uri <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt";
output.fn <- "cosmic-signature-spectra.tsv";

# download data table
x <- read.table(uri, sep="\t", header=TRUE, check.names=FALSE);

# remove empty columns
y <- x[, colnames(x) != ""];

# sort rows based on substituion type and trinucleotide context
y <- y[order(y[, 1], y[, 2]), ];

# visually check that mutation channels are in correct order
barplot(y[, "Signature 1"])  # spontaneous deamination at CpG
barplot(y[, "Signature 2"])  # APOBEC
barplot(y[, "Signature 4"])  # Smoking
barplot(y[, "Signature 6"])  # Mismatch repair
barplot(y[, "Signature 7"])  # Ultraviolet radiation

# output to file
write.table(y, output.fn, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
