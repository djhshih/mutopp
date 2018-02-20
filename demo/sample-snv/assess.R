library(io);
library(dplyr);

####

signature.fn <- "signature1.vtr";
sample.fn <- "sample_signature1.tsv";

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

# Jensen-Shannon Divergence using log2
# bounded in [0, 1]
jsd <- function(p, q) {
	# ensure that p and q are proper distributions
	p <- p / sum(p);
	q <- q / sum(q);

	m <- (p + q) / 2;

	0.5 * (
		sum( p * (log2(p) - log2(m)) ) + 
		sum( q * (log2(q) - log2(m)) )
	)
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

activities <- tapply(x$weight, x$channel, sum);
activities[is.na(activities)] <- 1e-4;

stopifnot(names(sig) == names(activities))

barplot(sig, las=2, cex.name=0.5);
barplot(activities, las=2, cex.name=0.5);

message("Jensen-Shannon Divergence: ", jsd(sig, activities))

