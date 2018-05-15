nsamples = 100
fasta = ../../data/reference.fa
gff3 = ../../data/genes.gff3

all: signature1.vtr signature2.vtr signature3.vtr signature4.vtr signature5.vtr signature6.vtr signature7.vtr signature8.vtr signature9.vtr signature10.vtr signature11.vtr signature12.vtr signature13.vtr signature14.vtr signature15.vtr signature16.vtr signature17.vtr signature18.vtr signature19.vtr signature20.vtr signature21.vtr signature22.vtr signature23.vtr signature24.vtr signature25.vtr signature26.vtr signature27.vtr signature28.vtr signature29.vtr signature30.vtr

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

signature8.vtr:
	cut -f 11 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature9.vtr:
	cut -f 12 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature10.vtr:
	cut -f 13 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature11.vtr:
	cut -f 14 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature12.vtr:
	cut -f 15 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature13.vtr:
	cut -f 16 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature14.vtr:
	cut -f 17 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature15.vtr:
	cut -f 18 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature16.vtr:
	cut -f 19 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature17.vtr:
	cut -f 20 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature18.vtr:
	cut -f 21 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature19.vtr:
	cut -f 22 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature20.vtr:
	cut -f 23 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature21.vtr:
	cut -f 24 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature22.vtr:
	cut -f 25 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature23.vtr:
	cut -f 26 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature24.vtr:
	cut -f 27 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature25.vtr:
	cut -f 28 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature26.vtr:
	cut -f 29 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature27.vtr:
	cut -f 30 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature28.vtr:
	cut -f 31 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature29.vtr:
	cut -f 32 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

signature30.vtr:
	cut -f 33 ../../data/cosmic-signature-spectra.tsv | sed 1d > $@

test: signature1.vtr signature4.vtr
	../../target/debug/mutopp sample-snv -n $(nsamples) --output-imp imp-sample_signature1.tsv $(fasta) $(gff3) signature1.vtr sample_signature1.tsv
	../../target/debug/mutopp sample-snv -n $(nsamples) --output-imp imp-sample_signature2.tsv $(fasta) $(gff3) signature1.vtr sample_signature2.tsv
	../../target/debug/mutopp sample-snv -n $(nsamples) --output-imp imp-sample_signature3.tsv $(fasta) $(gff3) signature1.vtr sample_signature3.tsv
	../../target/debug/mutopp sample-snv -n $(nsamples) --output-imp imp-sample_signature4.tsv $(fasta) $(gff3) signature4.vtr sample_signature4.tsv

