extern crate bio;

use std::io;
use std::str;
use std::fmt;
use std::ascii::AsciiExt;
use std::slice::Iter;
use std::env;
use std::process;
use std::collections::HashMap;

use bio::io::fasta;
use bio::io::gff;
use bio::utils::Strand;

#[derive(Debug)]
struct Gene {
    name: String,
    chromosome: String,
    start: u64,
    end: u64,
    strand: Strand,
    transcripts: HashMap<String, Transcript>,
}

#[derive(Debug)]
struct Transcript {
    start: u64,
    end: u64,
    /// Disjoint coding regions sorted by position
    coding_regions: Vec<Region>,
}

#[derive(Debug)]
struct Region {
    start: u64,
    end: u64,
}


fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        println!("usage: mutopp <fasta> <gff3>");
        process::exit(1);
    }

    let fasta = match fasta::Reader::from_file(&args[1]) {
        Err(why) => panic!("{:?}", why),
        Ok(reader) => reader,
    };

    for (_, r) in fasta.records().enumerate() {
        let record = r.expect("Error reading FASTA record");
        let desc = match record.desc() {
            Some(x) => x,
            None  => "",
        };

        println!("{}", record.id());
        println!("{}", desc);
        print_cds(record.seq(), true);
        print_protein(record.seq());

        let opp = accum_opp(&[record.seq()], 0).unwrap();
        print_opp_matrix(&opp);

        println!("");
    }

    // extract gene annotation from gff into struct Gene

    let mut gff = match gff::Reader::from_file(&args[2], gff::GffType::GFF3) {
        Err(why) => panic!("{:?}", why),
        Ok(reader) => reader,
    };

    let mut genes: HashMap<String, Gene> = HashMap::new();

    for r in gff.records() {
        let record = r.expect("Error reading GFF3 record");
        match record.feature_type() {
            "gene" => {
                // assume that frame is always 0
                genes.insert(
                    record.attributes().get("gene_id").unwrap().clone(),
                    Gene {
                        name: record.attributes().get("gene_name").unwrap().clone(),
                        chromosome: record.seqname().to_owned(),
                        start: *record.start(),
                        end: *record.end(),
                        strand: record.strand().unwrap(),
                        transcripts: HashMap::new(),
                    }
                );
            },
            "transcript" => {
                // assume that gene has already been inserted
                let gene = genes.get_mut(record.attributes().get("gene_id").unwrap()).expect("gene is inserted before descendent transcript");
                gene.transcripts.insert(
                    record.attributes().get("transcript_id").unwrap().clone(),
                    Transcript {
                        start: *record.start(),
                        end: *record.end(),
                        coding_regions: Vec::new(),
                    }
                );
            },
            "CDS" => {
                // assume that gene and transcript have already been inserted
                let gene = genes.get_mut(record.attributes().get("gene_id").unwrap()).expect("gene is inserted before descendent CDS");
                let transcript = gene.transcripts.get_mut(record.attributes().get("transcript_id").unwrap()).expect("transcript is inserted before descendent CDS");
                transcript.coding_regions.push(
                    Region { start: *record.start(), end: *record.end() }
                );
            },
            _ => {} // ignore all other fields
        }
    }

    println!("{:?}", genes);

    /*
    print_mutation_types();

    print_mutation_channels();

    let xs = mutation_channels();
    for (i, x) in xs.iter().enumerate() {
        println!("{} {}", i, x);
    }
    */
}

fn print_cds(x: &[u8], format: bool) {
    let len = x.len();
    assert!(len % 3 == 0);

    if format {
        let mut s = String::with_capacity(len + (len/3));
        let mut i = 0;
        let mut j = 0;  // counter for codon within a line
        let width = 80;
        for c in x {
            s.push(*c as char);
            i += 1;
            if i % 3 == 0 {
                s.push(' ');
                j += 1;
                if j * 4 > width {
                    s.push('\n');
                    j = 0;
                }
            }
        }
        println!("{}", s);
    } else {
        let s = str::from_utf8(x).expect("Found invalid DNA sequence");
        println!("{}", s);
    }

}

fn print_protein(dna: &[u8]) {
    let protein = translate(dna);
    let s = String::from_utf8(protein).expect("Found invalid protein sequence");
    println!("{}", s);
}

fn translate(dna: &[u8]) -> Vec<u8> {
    let len = dna.len();
    assert!(len % 3 == 0);

    let mut protein = Vec::with_capacity(len / 3);

    let mut i = 0;
    while i < len {
        if i == 0 {
            // first codon always code for methionine (AUG) or
            // formylmethionine, even for alternative start codons
            // (e.g. UUG, CUG, GUG, ACG, AUU and AGG)
            protein.push(b'M');
        } else {
            let codon: [u8;3] = [dna[i], dna[i+1], dna[i+2]];
            protein.push(codon_to_aa(&codon));
        }
        i += 3;
    }

    protein
}

// Standard genetic code
// Ala / A  GCT GCC GCA GCG
// Arg / R  CGT CGC CGA CGG AGA AGG
// Asn / N  AAT AAC
// Asp / D  GAT GAC
// Cys / C  TGT TGC
// Gln / Q  CAA CAG
// Glu / E  GAA GAG
// Gly / G  GGT GGC
// His / H  CAT CAC
// Ile / I  ATT ATC ATA
// Leu / L  TTA TTG CTT CTC CTA CTG
// Lys / K  AAA AAG
// Met / M  ATG
// Phe / F  TTT TTC
// Pro / P  CCT CCC CCA CCG
// Ser / S  TCT TCC TCA TCG AGT AGC
// Thr / T  ACT ACC ACA ACG
// Trp / W  TGG
// Tyr / Y  TAT TAC
// Val / V  GTT GTC GTA GTG
// Ter / *  TAA TGA TAG
fn codon_to_aa(codon: &[u8; 3]) -> u8 {
    match codon {
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"AAT" | b"AAC" => b'N',
        b"GAT" | b"GAC" => b'D',
        b"TGT" | b"TGC" => b'C',
        b"CAA" | b"CAG" => b'Q',
        b"GAA" | b"GAG" => b'E',
        b"GGT" | b"GGC" => b'G',
        b"CAT" | b"CAC" => b'H',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"AAA" | b"AAG" => b'K',
        b"ATG" => b'M',
        b"TTT" | b"TTC" => b'F',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"TGG" => b'W',
        b"TAT" | b"TAC" => b'Y',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TAA" | b"TGA" | b"TAG" => b'*',
        _ => b'X',
    }
}

const nucleotides: [u8; 4] = [b'A', b'C', b'G', b'T'];
const n_nucleotides: usize = 4;

const n_substitution_types: usize = n_nucleotides * (n_nucleotides - 1);

const n_mutation_channels: usize = n_substitution_types * n_nucleotides * n_nucleotides;

#[derive(Copy,Clone)]
enum MutationClass {
    Synonymous,
    Missense,
    Nonsense,
    SpliceSite,
}
const n_mutation_classes: usize = 4;

impl MutationClass {
    pub fn iter() -> Iter<'static, MutationClass> {
        use self::MutationClass::*;
        static classes: [MutationClass; n_mutation_classes] = [Synonymous, Missense, Nonsense, SpliceSite];
        classes.into_iter()
    }
}

impl fmt::Display for MutationClass {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let c = match *self {
            MutationClass::Synonymous => "syn",
            MutationClass::Missense => "mis",
            MutationClass::Nonsense => "nns",
            MutationClass::SpliceSite => "spl",
        };

        write!(f, "{}", c)
    }
}

const n_mutation_types: usize = (n_mutation_classes * n_mutation_channels) as usize;

// mutation channels are indexed by factors: impact class, stranded substitution, 5' context, and 3' context
// each factor serve as a subscript into the mutation channel array
// consider nucleotides in the order of A, C, G, T
// stranded substitutions (12): A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T, T>A, T>C, T>G
// 5' context (4): A, C, G, T
// 3' context (4): A, C, G, T
// number of channels = 12 * 4 * 4 = 192
// later factors vary faster than early factors
fn mutation_type_index(cl: MutationClass, nt_ref: u8, nt_alt: u8, nt_5p: u8, nt_3p: u8) -> usize {
    assert!(nt_ref != nt_alt);

    let i = cl as usize;

    const J: usize = n_substitution_types;
    let j = match nt_ref {
        b'A' => match nt_alt {
            b'C' => 0,
            b'G' => 1,
            b'T' => 2,
            _ => 0,
        },
        b'C' => match nt_alt {
            b'A' => 3,
            b'G' => 4,
            b'T' => 5,
            _ => 0,
        },
        b'G' => match nt_alt {
            b'A' => 6,
            b'C' => 7,
            b'T' => 8,
            _ => 0,
        },
        b'T' => match nt_alt {
            b'A' => 9,
            b'C' => 10,
            b'G' => 11,
            _ => 0,
        },
        _ => 0,
    };

    const K: usize = n_nucleotides;
    let k = match nt_5p {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    };

    const L: usize = n_nucleotides;
    let l = match nt_3p {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    };

    ((i * J + j) * K + k) * L + l
}

fn mutation_channels() -> Vec<String> {
    let cl = MutationClass::Synonymous;
    let mut names = vec![String::new(); n_mutation_channels as usize];
    for &nt_ref in nucleotides.iter() {
        for &nt_alt in nucleotides.iter() {
            if nt_ref != nt_alt {
                for &nt_5p in nucleotides.iter() {
                    for &nt_3p in nucleotides.iter() {
                        let idx = mutation_type_index(cl, nt_ref, nt_alt, nt_5p, nt_3p);
                        let context_5p = (nt_5p as char).to_ascii_lowercase();
                        let context_3p = (nt_3p as char).to_ascii_lowercase();
                        names[idx] = format!("{}{}{}>{}{}{}",
                            context_5p, nt_ref as char, context_3p,
                            context_5p, nt_alt as char, context_3p);
                    }
                }
            }
        }
    }
    names
}

fn print_mutation_channels() {
    let cl = MutationClass::Synonymous;
    for &nt_ref in nucleotides.iter() {
        for &nt_alt in nucleotides.iter() {
            if nt_ref != nt_alt {
                for &nt_5p in nucleotides.iter() {
                    for &nt_3p in nucleotides.iter() {
                        let idx = mutation_type_index(cl, nt_ref, nt_alt, nt_5p, nt_3p);
                        println!(
                            "{} {}>{} {} {}",
                            idx,
                            nt_ref as char, nt_alt as char,
                            nt_5p as char, nt_3p as char
                        );
                    }
                }
            }
        }
    }
}

fn print_mutation_types() {
    for &cl in MutationClass::iter() {
        for &nt_ref in nucleotides.iter() {
            for &nt_alt in nucleotides.iter() {
                if nt_ref != nt_alt {
                    for &nt_5p in nucleotides.iter() {
                        for &nt_3p in nucleotides.iter() {
                            let idx = mutation_type_index(cl, nt_ref, nt_alt, nt_5p, nt_3p);
                            println!(
                                "{} {} {}>{} {} {}",
                                idx,
                                cl,
                                nt_ref as char, nt_alt as char,
                                nt_5p as char, nt_3p as char
                            );
                        }
                    }
                }
            }
        }
    }
}

fn new_opp_matrix() -> [u32; n_mutation_types] {
    [0; n_mutation_types]
}

// print as n_mutation_channels by n_mutation_classes matrix
// matrix is stored in column-major order
fn print_opp_matrix(x: &[u32; n_mutation_types]) {
    const I: usize = n_mutation_channels;

    // print header
    print!("channel ");
    for &cl in MutationClass::iter() {
        // since the length of cl is not known at compile time, std::fmt does not
        // pad properly; therefore, pad manually
        print!(" {} ", cl);
    }
    println!("");

    let rownames = mutation_channels();
    for i in 0 .. I {
        print!("{} ", rownames[i]);
        for &cl in MutationClass::iter() {
            let j = cl as usize;
            print!("{:>4} ", x[j * I + i])
        }
        println!("");
    }

    // print class totals
    print!("# total ");
    let class_counts: Vec<u32> = MutationClass::iter().map(|&cl| sum_opp_in_class(x, cl)).collect();
    for c in class_counts {
        print!("{:>4} ", c);
    }
    println!("");
}

fn sum_opp_in_class(x: &[u32; n_mutation_types], cl: MutationClass) -> u32 {
    const I: usize = n_mutation_channels;
    let j = cl as usize;
    let mut y = 0u32;
    for i in 0 .. I {
        y += x[j * I + i];
    }
    y
}

/*
struct CodingSeqs {
    seqs: [&[u8]],
    padding: usize
}

impl CodingSeqs {
    fn prev_nt(cds_pos: usize) {

    }
    fn next_nt(cds_pos: usize) {

    }
    fn codon_at(cds_pos: usize) {

    }
}
*/

/// Accmulate mutation opportunities.
/// CDS must be a contiguous sequence in the genome (i.e. exons are not joined together).
/// Multiple CDS's should be analyzed by repeated calling this function.
/// If exons have been joined together to form the CDS, then mutations at exon-exon junctions will
/// be likely be incorrectly counted toward the wrong mutation context channel.
/// CDS may optionally include upstream and downstream sequences, the sizes of which are specified
/// by `padding`. These sequences are required for counting mutation
/// opportunities at the first and last positions, as well as counting mutation opportunities at
/// splice sites (first and last 2 bp of each intron).
fn accum_opp(seqs: &[&[u8]], padding: usize) -> Option<[u32; n_mutation_types]> {

    let nseqs = seqs.len();

    let mut x = new_opp_matrix();
    let mut offset = 0;
    let mut error = false;

    for s in 0 .. nseqs {
        let seq = seqs[s];
        let seq_len = seq.len();
        let cds_end = seq_len - padding;  // substract downstream padding
        let cds_len = cds_end - padding;  // substract upstream padding
        // CDS length might not be a mutiple of 3, because the CDS may span multiple exons, resulting
        // in mutiple disjoint CDS regions, which are processed one at a time

        let mut i = offset + padding;
        // step through the first nucleotide of each codon
        while i < cds_end {
            let codon_ref: [u8; 3] = [seq[i], seq[i+1], seq[i+2]];
            // translate codon
            let aa_ref = codon_to_aa(&codon_ref);
            // iterate through each position of the codon
            for j in 0 .. 3 {
                let ii = i + j;
                // ii-1 and ii+1 need to be valid indices
                if ii > 0 && ii < seq_len - 1 {
                    let context: [u8; 3] = [seq[ii-1], seq[ii], seq[ii+1]];
                    count_opp_at_codon_pos(&mut x, j, &codon_ref, aa_ref, &context);
                }
            }
            // go to next codon
            i += 3;
        }

        // consider the case in which last codon of the current CDS extends into the next CDS:
        // process this last codon by borrowing up to 2 nucleotides from the next CDS
        let remainder = (cds_len - offset) % 3;
        let needed = 3 - remainder;
        if remainder != 0 {
            if s == nseqs - 1 {
                // there are no more CDS... yet the last codon is complete
                // there is an error with the sequences
                error = true;
            } else {
                let seq2 = seqs[s+1];
                let cds2_len = seq2.len() - (2 * padding);
                if cds2_len < needed {
                    // there are not enough nucleotides in the next CDS to complete the codon... 
                    // this is an error with the sequences
                    error = true;
                } else {
                    // complete the last codon, using nucleotides from the next CDS
                    let i = cds_len - remainder;
                    let i2 = padding;

                    // remainder can only possibly be 1 or 2
                    match remainder {
                        1 => {
                            let codon_ref: [u8; 3] = [seq[i], seq2[i2], seq2[i2+1]];
                            let aa_ref = codon_to_aa(&codon_ref);
                            // process position 1 of the codon (which is in the current CDS)
                            {
                                let context: [u8; 3] = [seq[i-1], seq[i], seq[i+1]];
                                count_opp_at_codon_pos(&mut x, 0, &codon_ref, aa_ref, &context);
                            }
                            // process position 2 of the codon (which is in the next CDS)
                            {
                                let context: [u8; 3] = [seq2[i2-1], seq2[i2], seq2[i2+1]];
                                count_opp_at_codon_pos(&mut x, 1, &codon_ref, aa_ref, &context);
                            }
                            // process position 3 of the codon (which is in the next CDS)
                            {
                                let context: [u8; 3] = [seq2[i2], seq2[i2+1], seq2[i2+2]];
                                count_opp_at_codon_pos(&mut x, 2, &codon_ref, aa_ref, &context);
                            }
                        },
                        2 => {
                            let codon_ref: [u8; 3] = [seq[i], seq[i+1], seq2[i2]];
                            let aa_ref = codon_to_aa(&codon_ref);
                            // process position 1 of the codon (which is in the current CDS)
                            {
                                let context: [u8; 3] = [seq[i-1], seq[i], seq[i+1]];
                                count_opp_at_codon_pos(&mut x, 0, &codon_ref, aa_ref, &context);
                            }
                            // process position 2 of the codon (which is in the current CDS)
                            {
                                let context: [u8; 3] = [seq2[i], seq2[i+1], seq2[i+2]];
                                count_opp_at_codon_pos(&mut x, 1, &codon_ref, aa_ref, &context);
                            }
                            // process position 3 of the codon (which is in the next CDS)
                            {
                                let context: [u8; 3] = [seq2[i2-1], seq2[i2], seq2[i2+1]];
                                count_opp_at_codon_pos(&mut x, 2, &codon_ref, aa_ref, &context);
                            }
                        },
                        _ => panic!("remainder can only possibly be either 1 or 2")
                    } // match remainder
                } // cds2_len
            } // if s

            // skip nucleotides already processed for the next CDS
            offset = needed;
        } else {
            // the last codon of the current CDS is complete
            offset = 0;
        } // if remainder

    } // for seq

    if error {
        None
    } else {
        Some(x)
    }
}

/// x: opportunity matrix
/// pos: 0, 1, 2 within the codon
/// codon: codon sequence
/// context: sequence centered around the query position
fn count_opp_at_codon_pos(x: &mut [u32; n_mutation_types], pos: usize, codon: &[u8; 3], aa_ref: u8, context: &[u8; 3]) {
    let nt_ref = context[1];
    // iterate through possible substitutions
    for &nt_alt in nucleotides.iter() {
        if nt_alt != nt_ref {
            let codon_alt = mutate_codon(&codon, pos, nt_alt);
            let aa_alt = codon_to_aa(&codon_alt);
            let cl = if aa_alt == aa_ref {
                MutationClass::Synonymous 
            } else if aa_alt == b'*' {
                MutationClass::Nonsense
            } else {
                MutationClass::Missense
            };
            let idx = mutation_type_index(cl, nt_ref, nt_alt, context[0], context[2]);
            x[idx] += 1;
        }
    }
}

fn mutate_codon(codon: &[u8; 3], pos: usize, nt_alt: u8) -> [u8; 3] {
    assert!(pos <= 2);

    let mut codon_new = *codon;
    codon_new[pos] = nt_alt;

    codon_new
}

