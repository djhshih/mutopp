pub mod coding;
pub mod genomic;

use std::str;

use constants::*;

pub type Nucleotide = u8;
pub type Residue = u8;

pub type Codon = [Nucleotide; 3];

pub type DnaSeq = Vec<Nucleotide>;
pub type Peptide = Vec<Residue>;

pub fn complement(x: Nucleotide) -> Nucleotide {
    match x {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'n' => b'n',
        _ => b'.',
    }
}

/// Reverse complement sequnce in place.
pub fn reverse_complement(seq: &mut DnaSeq) {
    let n = seq.len();
    // iteratate from 0 to floor(n / 2) - 1
    for i in 0 .. (n / 2) {
        let j = n - 1 - i;
        // complement Nucleotide and reverse order in-place
        let x = complement(seq[i]);
        seq[i] = complement(seq[j]);
        seq[j] = x;
    }
    if n % 2 == 1 {
        // if sequence length is odd, then the lone middle element has not been touched
        let j = n / 2;
        seq[j] = complement(seq[j]);
    }
}

pub fn print_seq(x: &[Nucleotide]) {
    let s = str::from_utf8(x).expect("Found invalid DNA sequence");
    println!("{}", s);
}

pub fn print_seq_padded(x: &[Nucleotide], padding: usize) {
    let n = x.len();
    if n > padding * 2 {
        let upstream = str::from_utf8(&x[0..padding]).unwrap().to_lowercase();
        let core = str::from_utf8(&x[padding..(n-padding)]).unwrap().to_uppercase();
        let downstream = str::from_utf8(&x[(n-padding)..n]).unwrap().to_lowercase();
        println!("{}{}{}", upstream, core, downstream);
    } else {
        let s = str::from_utf8(x).unwrap().to_lowercase();
        println!("{}", s);
    }
}