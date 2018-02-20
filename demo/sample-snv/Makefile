
all: signature1.tsv

signature1.tsv:
	cut -f 4 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

test: signature1.tsv
	../../target/debug/mutopp sample-snv ../../data/reference.fa ../../data/genes.gff3 signature1.tsv output.tsv
