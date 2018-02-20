
all: signature1.vtr

signature1.vtr:
	cut -f 4 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

test: signature1.vtr
	../../target/debug/mutopp sample-snv ../../data/reference.fa ../../data/genes.gff3 signature1.vtr output.tsv
