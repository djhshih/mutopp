nsamples = 100
fasta = ../../data/reference.fa
gff3 = ../../data/genes.gff3

all: signature1.vtr signature2.vtr signature3.vtr signature4.vtr signature5.vtr signature6.vtr signature7.vtr

signature1.vtr:
	cut -f 4 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature2.vtr:
	cut -f 5 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature3.vtr:
	cut -f 6 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature4.vtr:
	cut -f 7 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature5.vtr:
	cut -f 8 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature6.vtr:
	cut -f 9 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature7.vtr:
	cut -f 10 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

test: signature1.vtr
	../../target/debug/mutopp sample-snv -n $(nsamples) --output-imp imp-sample_signature1.tsv $(fasta) $(gff3) signature1.vtr sample_signature1.tsv

