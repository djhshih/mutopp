#!/usr/bin/env Rscript

library(io, quietly=TRUE, warn=FALSE);
library(argparser, quietly=TRUE, warn=FALSE);

####

pr <- arg_parser("Assess an importance sample of mutations according to signature spectrum");
pr <- add_argument(pr, "signature", help="input signature spectrum");
pr <- add_argument(pr, "sample", help="input importance sample");

argv <- parse_args(pr);

signature.fn <- argv$signature;
sample.fn <- argv$sample;

out.fname <- as.filename(sample.fn);
pdf.fname <- set_fext(out.fname, ext="pdf");

###

nucleotides <- c("A", "C", "G", "T");

channel.levels <- c(
	paste0("C>A_", rep(nucleotides, each = 4), "C", rep(nucleotides, times=4)),
	paste0("C>G_", rep(nucleotides, each = 4), "C", rep(nucleotides, times=4)),
	paste0("C>T_", rep(nucleotides, each = 4), "C", rep(nucleotides, times=4)),
	paste0("T>A_", rep(nucleotides, each = 4), "T", rep(nucleotides, times=4)),
	paste0("T>C_", rep(nucleotides, each = 4), "T", rep(nucleotides, times=4)),
	paste0("T>G_", rep(nucleotides, each = 4), "T", rep(nucleotides, times=4))
);

# Complement a nucleotide.
complement <- function(x) {
	unlist(lapply(x,
		function(z) {
			if (z == "A") {
				"T"
			} else if (z == "C") {
				"G"
			} else if (z == "G") {
				"C"
			} else if (z == "T") {
				"A"
			}
		}
	));
}

# entropy using log2
entropy <- function(v) {
	v <- v[v > 0]
	return(sum(-v * logw(v)))
}

# Cosine similarity
# bounded in [-1, 1]
# higher value => more similar
cos_sim <- function(x, y) {
	sum(x * y) / sqrt(sum(x * x) * sum(y * y))
}

kld <- function(p, q) {
	sum( ifelse(p == 0, 0, p * (log2(p) - log2(q))) )
}

# Jensen-Shannon Divergence using log2
# bounded in [0, 1]
# lower value => more similar
jsd <- function(p, q) {
	# ensure that p and q are proper distributions
	p <- p / sum(p);
	q <- q / sum(q);

	m <- 0.5 * (p + q);

	0.5 * (kld(p, m) + kld(q, m))
}

####

# COSMIC mutation signature spectrum
sig <- qread(signature.fn);
names(sig) <- channel.levels;
sig.df <- data.frame(
	channel = channel.levels,
	probability = sig
);

# read weighted samples
x <- qread(sample.fn);


flip <- x$ref %in% c("C", "T");
ref_ns <- ifelse(flip, as.character(x$ref), complement(x$ref));
alt_ns <- ifelse(flip, as.character(x$alt), complement(x$alt));

mutation <- paste0(ref_ns, ">", alt_ns);

x$channel <- factor(paste0(mutation, "_", x$context), levels=channel.levels);
stopifnot(!is.na(x$channel))

message("Importance sampling")

activities.is <- tapply(exp(x$lweight), x$channel, sum);
stopifnot(names(sig) == names(activities.is))
activities.is[is.na(activities.is)] <- 0;

message("Number of empty channels: ", sum(activities.is == 0));

activities.is <- activities.is / sum(activities.is);


qdraw({barplot(sig, las=2, cex.name=0.5)}, width = 12, height = 4, file=insert(pdf.fname, "spectrum"));
qdraw({barplot(activities.is, las=2, cex.name=0.5)}, width = 12, height = 4, file=insert(pdf.fname, c("spectrum", "is")));

message("Cosine similarity: ", cos_sim(sig, activities.is))
message("Jensen-Shannon Divergence: ", jsd(sig, activities.is))

message("");

####

message("Sampling importance resampling")

# subsample the importance samples
b <- 0.01;
nsubsample <- ceiling(b * nrow(x));
norm.weights <- exp(x$lweight - log(sum(exp(x$lweight))));
y <- x[sample.int(nrow(x), nsubsample, replace=FALSE, prob = norm.weights), ];
qwrite(y[, ! colnames(y) %in% c("weight", "channel")], insert(out.fname, "sir"));

activities.sir <- table(y$channel);
stopifnot(names(sig) == names(activities.sir))
message("Number of empty channels: ", sum(activities.sir == 0));

activities.sir <- activities.sir / sum(activities.sir);

qdraw({barplot(activities.sir, las=2, cex.name=0.5)}, width = 12, height = 4, file=insert(pdf.fname, c("spectrum", "sir")));

message("Cosine similarity: ", cos_sim(sig, activities.sir))
message("Jensen-Shannon Divergence: ", jsd(sig, activities.sir))

message("");

####

message("Idealized direct sampling")

activities.ds <- table(factor(sample.int(length(sig), nsubsample, replace=TRUE, prob=sig), levels=1:length(channel.levels), label=channel.levels));
stopifnot(names(sig) == names(activities.ds))
message("Number of empty channels: ", sum(activities.ds == 0));

activities.ds <- activities.ds / sum(activities.ds);

qdraw({barplot(activities.ds, las=2, cex.name=0.5)}, width = 12, height = 4, file=insert(pdf.fname, c("spectrum", "ds")));

message("Cosine similarity: ", cos_sim(sig, activities.ds))
message("Jensen-Shannon Divergence: ", jsd(sig, activities.ds))

message("");

