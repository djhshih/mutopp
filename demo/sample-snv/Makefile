
all: signature1.vtr

signature1.vtr:
	cut -f 4 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

test: signature1.vtr
	../../target/debug/mutopp sample-snv -n 100 ../../data/reference.fa ../../data/genes.gff3 signature1.vtr sample_signature1.tsv

