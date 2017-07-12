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

const nucleotides: [u8; 4] = [b'A', b'C', b'G', b'T'];
const n_nucleotides: usize = 4;

const n_substitution_types: usize = n_nucleotides * (n_nucleotides - 1);

const n_mutation_channels: usize = n_substitution_types * n_nucleotides * n_nucleotides;

#[derive(Copy,Clone)]
enum MutationClass {
    Synonymous,
    Missense,
    StartOrStop,
    SpliceSite,
}
const n_mutation_classes: usize = 4;

type Seq = Vec<u8>;

struct CodingSequence {
    seqs: Vec<Seq>,
    padding: usize,
}

type MutOppArray = [u32; n_mutation_types];


fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 4 {
        println!("usage: mutopp <fasta> <coord> <gff3>");
        process::exit(1);
    }

    let fasta_fn = &args[1];
    let coord = &args[2];
    let gff_fn = &args[3];

    // read the entire fasta
    /*
    let fasta = match fasta::Reader::from_file(fasta_fn) {
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
        print_cds(record.seq());

        /*
        print_protein(record.seq());

        let opp = count_opp(&[record.seq()], 0).unwrap();
        print_opp_matrix(&opp);
        */

        println!("");
    }
    */

    // extract gene annotation from gff into struct Gene

    let mut gff = match gff::Reader::from_file(gff_fn, gff::GffType::GFF3) {
        Err(why) => panic!("{:?}", why),
        Ok(reader) => reader,
    };


    let genes = genes_from_gff(&mut gff);

    println!("{:?}", genes);

    // read only part of the fasta
    let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
        Err(why) => panic!("{:?}", why),
        Ok(reader) => reader,
    };

    if let Some((seqname, start, end)) = parse_coordinate(coord) {
        println!("{} {} {}", seqname, start, end);
        let mut seq: Vec<u8> = Vec::new();
        match ifasta.read(&seqname, start, end, &mut seq) {
            Err(why) => panic!("{:?}", why),
            Ok(()) => print_seq(&seq),
        }
    } else {
        panic!("Invalid coordinate: {}", coord);
    }

    /*
    let t = Transcript {
       start: 0,
       end: 1403,
       coding_regions: vec![Region{start: 4, end: 10}, Region{start: 15, end: 18}],
    };
    */

    for (gid, gene) in genes.iter() {
        println!("{}", gid);
        for (tid, transcript) in gene.transcripts.iter() {
            println!("{}", tid);
            let padding = 3;
            let cds = get_transcript_cds_from_fasta(&mut ifasta, transcript, &gene.chromosome, padding);
            for seq in &cds.seqs {
                print_seq(seq);
            }
            let opp = count_opp(&cds.seqs, cds.padding).unwrap();
            print_opp_matrix(&opp);
        }
    }

    /*
    print_mutation_types();

    print_mutation_channels();

    let xs = mutation_channels();
    for (i, x) in xs.iter().enumerate() {
        println!("{} {}", i, x);
    }
    */
}

fn get_transcript_cds_from_fasta(reader: &mut FastaIndexedReader, transcript: &Transcript, chromosome: &str, padding: usize) -> CodingSequence {
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    for region in &transcript.coding_regions {
        let mut seq: Vec<u8> = Vec::new();
        match reader.read(chromosome, region.start - padding as u64, region.end + padding as u64, &mut seq) {
            Err(why) => panic!("{:?}", why),
            Ok(()) => seqs.push(seq),
        }
    }
    
    CodingSequence { seqs: seqs, padding: padding }
}

type FastaIndexedReader = fasta::IndexedReader<std::fs::File>;
type GffReader = gff::Reader<std::fs::File>;

/// Construct genes from GFF reader.
///
/// GFF coordinates are 1-based, closed. Convert the coordinates to 0-based, half-open.
fn genes_from_gff(reader: &mut GffReader) -> HashMap<String, Gene> {
    let mut genes: HashMap<String, Gene> = HashMap::new();

    for r in reader.records() {
        let record = r.expect("Error reading GFF3 record");
        match record.feature_type() {
            "gene" => {
                // assume that frame is always 0
                genes.insert(
                    record.attributes().get("gene_id").unwrap().clone(),
                    Gene {
                        name: record.attributes().get("gene_name").unwrap().clone(),
                        chromosome: record.seqname().to_owned(),
                        start: *record.start() - 1,
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
                        start: *record.start() - 1,
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
                    Region { start: *record.start() - 1, end: *record.end() }
                );
            },
            _ => {} // ignore all other fields
        }
    }

    genes
}

/// Parse genomic coordinate.
///
/// Assume coordinate is specified in 1-based, closed format.
/// Internally, coordinate is converted to 0-based, half-open format.
///
fn parse_coordinate(x: &str) -> Option<(String, u64, u64)> {
    match x.find(':') {
        Some(i) => {
            let seqname = &x[..i];
            let range = &x[i+1..];
            match range.find('-') {
                Some(j) => {
                    let start = &range[..j];
                    let end = &range[j+1..];
                    if let Ok(start) = start.parse() {
                        if let Ok(end) = end.parse() {
                            if start > 0 && end > start {
                                Some((String::from(seqname), start-1, end))
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                },
                None => None,
            }
        },
        None => None,
    }
}

fn print_seq(x: &[u8]) {
    let s = str::from_utf8(x).expect("Found invalid DNA sequence");
    println!("{}", s);
}

fn print_cds(x: &[u8]) {
    let len = x.len();
    assert!(len % 3 == 0);

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
// Ter / $  TAA TGA TAG
// Additionally, '^' will be used to denote the start amino acid in other parts
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
        b"TAA" | b"TGA" | b"TAG" => b'$',
        _ => b'X',
    }
}


impl MutationClass {
    pub fn iter() -> Iter<'static, MutationClass> {
        use self::MutationClass::*;
        static classes: [MutationClass; n_mutation_classes] = [Synonymous, Missense, StartOrStop, SpliceSite];
        classes.into_iter()
    }
}

impl fmt::Display for MutationClass {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let c = match *self {
            MutationClass::Synonymous => "syn",
            MutationClass::Missense => "mis",
            MutationClass::StartOrStop => "sos",
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

fn new_opp_matrix() -> MutOppArray {
    [0; n_mutation_types]
}

// print as n_mutation_channels by n_mutation_classes matrix
// matrix is stored in column-major order
fn print_opp_matrix(x: &MutOppArray) {
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

fn sum_opp_in_class(x: &MutOppArray, cl: MutationClass) -> u32 {
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

fn cds_is_valid(seqs: &Vec<Seq>, padding: usize) -> bool {
    let mut valid = true;

    let mut total_cds_len = 0;

    let nseqs = seqs.len();
    for s in 0 .. nseqs {
        let seq = &seqs[s];
        let seq_len = seq.len();
        if seq_len <= 2 * padding {
            // after removing padding, no CDS sequence is left
            valid = false;
            break;
        }
        // substract upstream and upstream padding
        let cds_len = seq_len - (2 * padding);  
        total_cds_len += cds_len;
    }

    // joined CDS should be a multiple of 3
    if total_cds_len % 3 != 0 {
        valid = false;
    }

    // check for presense of stop codon at the final position using the standard genetic code
    {
        let seq = &seqs[nseqs-1];
        let end = seq.len() - padding;
        let codon = [seq[end-3], seq[end-2], seq[end-1]];
        if codon_to_aa(&codon) != b'$' {
            valid = false;
        }
    }
     
    // not checked: presence of premature stop codon due to computation cost (TODO)
    // not checked: start codon sequence, because there are many alternative start codons

    valid
}


/// Count mutation opportunities.
///
/// Each CDS must be a contiguous sequence in the genome (i.e. exons are not joined together).
/// If exons have been joined together to form the CDS, then mutations at exon-exon junctions will
/// be likely be incorrectly counted toward the wrong mutation context channel.
/// CDS may optionally include upstream and downstream sequences, the sizes of which are specified
/// by `padding`. These sequences are required for counting mutation
/// opportunities at the first and last positions, as well as counting mutation opportunities at
/// splice sites (first and last 2 bp of each intron).
/// Splice sites will not be counted unless `padding >= 3`, because splice sites are the first
/// and last 2 bp of each intron and at least one nucleotide context is needed to count the last
/// nucleotide of the splice donor and the first nucleotide of the splice acceptor.
fn count_opp(seqs: &Vec<Vec<u8>>, padding: usize) -> Option<MutOppArray> {

    let nseqs = seqs.len();

    let mut x = new_opp_matrix();
    let mut offset = 0;
    let mut error = false;

    // count mutations in each CDS

    for s in 0 .. nseqs {
        let seq = &seqs[s];
        let seq_len = seq.len();
        let cds_end = seq_len - padding;  // substract downstream padding
        let cds_len = cds_end - padding;  // substract upstream padding
        // CDS length might not be a mutiple of 3, because the CDS may span multiple exons, resulting
        // in mutiple disjoint CDS regions, which are processed one at a time

        let cds_begin = offset + padding;
        let mut i = cds_begin;
        // step through the first nucleotide of each codon
        while i < cds_end - 2 {
            let codon_ref: [u8; 3] = [seq[i], seq[i+1], seq[i+2]];
            // translate codon
            let aa_ref = if s == 0 && i == cds_begin {
                // first codon of the first CDS is the start codon
                // this is usually ATG but alternative start codons are possible
                b'^'
            } else {
                codon_to_aa(&codon_ref)
            };
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
                break;
            } else {
                let seq2 = &seqs[s+1];
                let cds2_len = seq2.len() - (2 * padding);
                if cds2_len < needed {
                    // there are not enough nucleotides in the next CDS to complete the codon... 
                    // this is an error with the sequences
                    error = true;
                    break;
                } else {
                    // complete the last codon, using nucleotides from the next CDS
                    let i = cds_end - remainder;
                    let i2 = padding;

                    // remainder can only possibly be 1 or 2
                    match remainder {
                        1 => {
                            let codon_ref: [u8; 3] = [seq[i], seq2[i2], seq2[i2+1]];
                            let aa_ref = codon_to_aa(&codon_ref);
                            // process position 1 of the codon (which is in the current CDS)
                            if i < seq_len - 1 {
                                let context: [u8; 3] = [seq[i-1], seq[i], seq[i+1]];
                                count_opp_at_codon_pos(&mut x, 0, &codon_ref, aa_ref, &context);
                            }
                            // process position 2 of the codon (which is in the next CDS)
                            if i2 > 0 {
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
                            if i < seq_len - 1 {
                                let context: [u8; 3] = [seq[i-1], seq[i], seq[i+1]];
                                count_opp_at_codon_pos(&mut x, 0, &codon_ref, aa_ref, &context);
                            }
                            // process position 2 of the codon (which is in the current CDS)
                            if i < seq_len - 2 {
                                let context: [u8; 3] = [seq[i], seq[i+1], seq[i+2]];
                                count_opp_at_codon_pos(&mut x, 1, &codon_ref, aa_ref, &context);
                            }
                            // process position 3 of the codon (which is in the next CDS)
                            if i2 > 0 {
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

    // count mutations in the splice site
    if !error {
        count_opp_at_splice_sites(&mut x, seqs, padding)
    }

    if error {
        None
    } else {
        Some(x)
    }
}

fn count_opp_at_splice_sites(mut x: &mut MutOppArray, seqs: &Vec<Vec<u8>>, padding: usize) {
    if padding >= 3 {
        let nseqs = seqs.len();
        for s in 0 .. nseqs {
            let seq = &seqs[s];

            // process splice acceptor site if it is not the first CDS
            if s != 0 {
                count_opp_at_slice_site(&mut x, &seq, padding - 2);
            }

            // process splice donor site if it is not the last CDS
            if s != nseqs - 1 {
                count_opp_at_slice_site(&mut x, &seq, seq.len() - padding);
            }
        }
    }
}

fn count_opp_at_slice_site(x: &mut MutOppArray, seq: &[u8], begin: usize) {
    let splice =  MutationClass::SpliceSite;
    // iterate through the 2 positions of the splice site
    for j in 0 .. 2 {
        let context_5p = seq[begin-1+j];
        let nt_ref = seq[begin+j];
        let context_3p = seq[begin+1+j];
        for &nt_alt in nucleotides.iter() {
            if nt_alt != nt_ref {
                let idx = mutation_type_index(splice, nt_ref, nt_alt, context_5p, context_3p);
                x[idx] += 1;
            }
        }
    }
}

/// x: opportunity matrix
/// pos: 0, 1, 2 within the codon
/// codon: codon sequence
/// context: sequence centered around the query position
fn count_opp_at_codon_pos(x: &mut MutOppArray, pos: usize, codon: &[u8; 3], aa_ref: u8, context: &[u8; 3]) {
    let nt_ref = context[1];
    // iterate through possible substitutions
    for &nt_alt in nucleotides.iter() {
        if nt_alt != nt_ref {
            let cl = if aa_ref == b'^' {
                // assume that any mutation at start codon causes start codon loss
                // this is true for ATG start codon (no degeneracy) and
                // probably true for alternative start codons as well
                MutationClass::StartOrStop
            } else {
                let codon_alt = mutate_codon(&codon, pos, nt_alt);
                let aa_alt = codon_to_aa(&codon_alt);
                if aa_alt == aa_ref {
                    MutationClass::Synonymous 
                } else if aa_ref == b'$' || aa_alt == b'$' {
                    MutationClass::StartOrStop
                } else {
                    MutationClass::Missense
                }
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_count_opp_pad0() {
        // sequences contain only CDS without any padding
        let test = vec![
            b"ATGATG".to_vec(),
            b"ATGATGATG".to_vec(),
            b"ATGATGATGATGTAG".to_vec(),
        ];
        let out = count_opp(&test, 0).expect("Invalid opportunity array");
        
        // since no 5' or 3' padding is available, the first and last nucleotides
        // of each CDS have no context, and they are therefore not counted

        let mut ans = new_opp_matrix();

        // mutations in position 2 of start codon (ATG) with context ATG
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // mutations in position 3 of start codon (ATG) with context TGA
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'A', b'T', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'C', b'T', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'T', b'T', b'A')] = 1;

        // mutations in TGA context (except involving start codon)
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'A')] = 5;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'A')] = 5;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'A')] = 5;

        // mutations in GAT context
        ans[mutation_type_index(MutationClass::Missense, b'A', b'C', b'G', b'T')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'G', b'G', b'T')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'T', b'G', b'T')] = 6;

        // mutations in ATG context (except involving start codon)
        ans[mutation_type_index(MutationClass::Missense, b'T', b'A', b'A', b'G')] = 8;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'C', b'A', b'G')] = 8;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'G', b'A', b'G')] = 8;

        // mutation in position 3 of penultimate codon (ATG) with context TGT
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'T')] = 1;

        // mutation in position 1 of stop codon (TAG) with context GTA
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'G', b'A')] = 1;

        // mutation in position 2 of stop codon (TAG) with context TAG
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        // discard length information from type s.t. compiler does not complain
        let pout: &[u32] = &out;
        let pans: &[u32] = &ans;

        assert_eq!(pout, pans);
    }

    #[test]
    fn test_count_opp_pad1() {
        // sequences contain CDS with 1 nucleotide padding
        let test = vec![
            b"CATGATGC".to_vec(),
            b"CATGATGATGC".to_vec(),
            b"CATGATGATGATGTAGC".to_vec(),
        ];
        let out = count_opp(&test, 1).expect("Invalid opportunity array");
        
        // the result should be similar to the first test,
        // except now that first and last nucleotides are each CDS have
        // 5' and 3' context and are thus counted

        let mut ans = new_opp_matrix();

        // mutations in position 2 of start codon (ATG) with context ATG
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // mutations in position 3 of start codon (ATG) with context TGA
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'A', b'T', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'C', b'T', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'T', b'T', b'A')] = 1;

        // mutations in TGA context (except involving start codon)
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'A')] = 5;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'A')] = 5;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'A')] = 5;

        // mutations in GAT context
        ans[mutation_type_index(MutationClass::Missense, b'A', b'C', b'G', b'T')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'G', b'G', b'T')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'T', b'G', b'T')] = 6;

        // mutations in ATG context (except involving start codon)
        ans[mutation_type_index(MutationClass::Missense, b'T', b'A', b'A', b'G')] = 8;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'C', b'A', b'G')] = 8;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'G', b'A', b'G')] = 8;

        // mutation in position 3 of penultimate codon (ATG) with context TGT
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'T')] = 1;

        // mutation in position 1 of stop codon (TAG) with context GTA
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'G', b'A')] = 1;

        // mutation in position 2 of stop codon (TAG) with context TAG
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        // additional mutations counted in pad1 but not in pad0

        // mutaiton in first nucleotide of start codon with context CAT
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'C', b'C', b'T')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'G', b'C', b'T')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'T', b'C', b'T')] = 1;

        // mutaiton in first nucleotide of CDS 2 and 3 with context CAT
        ans[mutation_type_index(MutationClass::Missense, b'A', b'C', b'C', b'T')] = 2;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'G', b'C', b'T')] = 2;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'T', b'C', b'T')] = 2;

        // mutation last nucleotide of CDS 1 and 2 with context TGC
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'C')] = 2;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'C')] = 2;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'C')] = 2;

        // mutation in position 3 of stop codon (TAG -> TAA) with context AGC
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'C', b'A', b'C')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'T', b'A', b'C')] = 1;

        // synonymous mutation in position 3 of stop codon (TAG -> TAA) with context AGC
        ans[mutation_type_index(MutationClass::Synonymous, b'G', b'A', b'A', b'C')] = 1;

        // discard length information from type s.t. compiler does not complain
        let pout: &[u32] = &out;
        let pans: &[u32] = &ans;

        assert_eq!(pout, pans);
    }

    #[test]
    fn test_count_opp_pad0_split() {
        // sequences contain only CDS without any padding but with split codons
        let test = vec![
            b"ATGATGA".to_vec(),
            b"TGATGATGAT".to_vec(),
            b"GATGATGATGTAG".to_vec(),
        ];
        let out = count_opp(&test, 0).expect("Invalid opportunity array");
        
        // since no 5' or 3' padding is available, the first and last nucleotides
        // of each CDS have no context, and they are therefore not counted

        let mut ans = new_opp_matrix();

        // mutations in position 2 of start codon (ATG) with context ATG
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // mutations in position 3 of start codon (ATG) with context TGA
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'A', b'T', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'C', b'T', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'T', b'T', b'A')] = 1;

        // mutations in TGA context (except involving start codon)
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'A')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'A')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'A')] = 6;

        // mutations in GAT context
        ans[mutation_type_index(MutationClass::Missense, b'A', b'C', b'G', b'T')] = 7;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'G', b'G', b'T')] = 7;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'T', b'G', b'T')] = 7;

        // mutations in ATG context (except involving start codon)
        ans[mutation_type_index(MutationClass::Missense, b'T', b'A', b'A', b'G')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'C', b'A', b'G')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'G', b'A', b'G')] = 6;

        // mutation in position 3 of penultimate codon (ATG) with context TGT
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'T')] = 1;

        // mutation in position 1 of stop codon (TAG) with context GTA
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'G', b'A')] = 1;

        // mutation in position 2 of stop codon (TAG) with context TAG
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        // discard length information from type s.t. compiler does not complain
        let pout: &[u32] = &out;
        let pans: &[u32] = &ans;

        assert_eq!(pout, pans);
    }

    #[test]
    fn test_count_opp_pad1_split() {
        // sequences contain CDS with 1 nucleotide padding and with split codons
        let test = vec![
            b"CATGATGAC".to_vec(),
            b"CTGATGATGATC".to_vec(),
            b"CGATGATGATGTAGC".to_vec(),
        ];
        let out = count_opp(&test, 1).expect("Invalid opportunity array");
        
        // the result should be similar to the first test,
        // except now that first and last nucleotides are each CDS have
        // 5' and 3' context and are thus counted

        let mut ans = new_opp_matrix();

        // mutations in position 2 of start codon (ATG) with context ATG
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // mutations in position 3 of start codon (ATG) with context TGA
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'A', b'T', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'C', b'T', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'T', b'T', b'A')] = 1;

        // mutations in TGA context (except involving start codon)
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'A')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'A')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'A')] = 6;

        // mutations in GAT context
        ans[mutation_type_index(MutationClass::Missense, b'A', b'C', b'G', b'T')] = 7;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'G', b'G', b'T')] = 7;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'T', b'G', b'T')] = 7;

        // mutations in ATG context (except involving start codon)
        ans[mutation_type_index(MutationClass::Missense, b'T', b'A', b'A', b'G')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'C', b'A', b'G')] = 6;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'G', b'A', b'G')] = 6;

        // mutation in position 3 of penultimate codon (ATG) with context TGT
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'T', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'T', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'T', b'T')] = 1;

        // mutation in position 1 of stop codon (TAG) with context GTA
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'G', b'A')] = 1;

        // mutation in position 2 of stop codon (TAG) with context TAG
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        // additional mutations counted in pad1 but not in pad0

        // mutaiton in first nucleotide of start codon with context CAT
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'C', b'C', b'T')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'G', b'C', b'T')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'T', b'C', b'T')] = 1;

        // mutaiton in first nucleotide of CDS 2 with context CTG
        ans[mutation_type_index(MutationClass::Missense, b'T', b'A', b'C', b'G')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'C', b'C', b'G')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'G', b'C', b'G')] = 1;

        // mutaiton in first nucleotide of CDS 3 with context CGA
        ans[mutation_type_index(MutationClass::Missense, b'G', b'A', b'C', b'A')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'C', b'C', b'A')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'G', b'T', b'C', b'A')] = 1;

        // mutation last nucleotide of CDS 1 with context GAC
        ans[mutation_type_index(MutationClass::Missense, b'A', b'C', b'G', b'C')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'G', b'G', b'C')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'A', b'T', b'G', b'C')] = 1;

        // mutation last nucleotide of CDS 2 with context ATC
        ans[mutation_type_index(MutationClass::Missense, b'T', b'A', b'A', b'C')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'C', b'A', b'C')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'T', b'G', b'A', b'C')] = 1;

        // mutation in position 3 of stop codon (TAG -> TAA) with context AGC
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'C', b'A', b'C')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'T', b'A', b'C')] = 1;

        // synonymous mutation in position 3 of stop codon (TAG -> TAA) with context AGC
        ans[mutation_type_index(MutationClass::Synonymous, b'G', b'A', b'A', b'C')] = 1;

        // discard length information from type s.t. compiler does not complain
        let pout: &[u32] = &out;
        let pans: &[u32] = &ans;

        assert_eq!(pout, pans);
    }

    #[test]
    fn test_count_opp_splice() {
        // sequences contain CDS with 3 nucleotide padding,
        // permitting splice site mutation opportunity counting
        let test = vec![
            b"ACCATGCCCCCCGTA".to_vec(),
            b"CAGCCCCCCGTG".to_vec(),
            b"TAGCCCCCCGTA".to_vec(),
            b"AAGCCCCCCGTG".to_vec(),
            b"GCACCCCCCTAGCCC".to_vec(),
        ];
        let out = count_opp(&test, 3).expect("Invalid opportunity array");

        let mut ans = new_opp_matrix();

        // donor splice site mutations

        // two cGTa and two cGTg donor splice sites, mutation at position 1  with context CGT
        ans[mutation_type_index(MutationClass::SpliceSite, b'G', b'A', b'C', b'T')] = 4;
        ans[mutation_type_index(MutationClass::SpliceSite, b'G', b'C', b'C', b'T')] = 4;
        ans[mutation_type_index(MutationClass::SpliceSite, b'G', b'T', b'C', b'T')] = 4;

        // two cGTa donor splice sites, mutation at position 2 with context GTA
        ans[mutation_type_index(MutationClass::SpliceSite, b'T', b'A', b'G', b'A')] = 2;
        ans[mutation_type_index(MutationClass::SpliceSite, b'T', b'C', b'G', b'A')] = 2;
        ans[mutation_type_index(MutationClass::SpliceSite, b'T', b'G', b'G', b'A')] = 2;

        // two cGTg donor splice sites, mutation at position 2 with context GTG
        ans[mutation_type_index(MutationClass::SpliceSite, b'T', b'A', b'G', b'G')] = 2;
        ans[mutation_type_index(MutationClass::SpliceSite, b'T', b'C', b'G', b'G')] = 2;
        ans[mutation_type_index(MutationClass::SpliceSite, b'T', b'G', b'G', b'G')] = 2;

        // acceptor splice site mutations

        // cAGc acceptor site, mutation at position 1 with context CAG
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'C', b'C', b'G')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'T', b'C', b'G')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'G', b'C', b'G')] = 1;

        // tAGc acceptor site, mutation at position 1 with context TAG
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'C', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'T', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'G', b'T', b'G')] = 1;

        // aAGc acceptor site, mutation at position 1 with context AAG
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'C', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'T', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'G', b'A', b'G')] = 1;

        // cAGc, tAGc, and aAGc acceptor sites, mutation at position 2 with context AGC
        ans[mutation_type_index(MutationClass::SpliceSite, b'G', b'A', b'A', b'C')] = 3;
        ans[mutation_type_index(MutationClass::SpliceSite, b'G', b'C', b'A', b'C')] = 3;
        ans[mutation_type_index(MutationClass::SpliceSite, b'G', b'T', b'A', b'C')] = 3;

        // gCAc acceptor site, mutation at position 1 with context GCA
        ans[mutation_type_index(MutationClass::SpliceSite, b'C', b'A', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'C', b'T', b'G', b'A')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'C', b'G', b'G', b'A')] = 1;

        // gCAc acceptor site, mutation at position 2 with context CAC
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'C', b'C', b'C')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'G', b'C', b'C')] = 1;
        ans[mutation_type_index(MutationClass::SpliceSite, b'A', b'T', b'C', b'C')] = 1;

        // start codon mutations

        // start codon mutations, position 1 with context CAT
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'C', b'C', b'T')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'G', b'C', b'T')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'T', b'C', b'T')] = 1;

        // start codon mutations, position 2 with context ATG
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // start codon mutations, position 3 with context TGC
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'A', b'T', b'C')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'C', b'T', b'C')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'T', b'T', b'C')] = 1;

        // stop codon mutations, position 1 with context CTA
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'A', b'C', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'C', b'C', b'A')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'T', b'G', b'C', b'A')] = 1;

        // stop codon mutations, position 2 with context TAG
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        // stop codon mutations, position 3 with context AGC
        // TAG -> TAA is synonymous stop codon mutation
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'C', b'A', b'C')] = 1;
        ans[mutation_type_index(MutationClass::StartOrStop, b'G', b'T', b'A', b'C')] = 1;

        // synonymous mutations
        
        // stop codon mutations, position 3 with context AGC
        // TAG -> TAA is synonymous stop codon mutation
        ans[mutation_type_index(MutationClass::Synonymous, b'G', b'A', b'A', b'C')] = 1;
        
        // position 3 mutations in codon CCC are synonymous
        // CCC -> CCA, CCG, CCT (all code for Proline)

        // context CCC (first CCC codon of ech of 5 CDS's)
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'A', b'C', b'C')] = 5;
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'G', b'C', b'C')] = 5;
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'T', b'C', b'C')] = 5;
        
        // context CCG (last codon of each of 4 CDS's)
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'A', b'C', b'G')] = 4;
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'G', b'C', b'G')] = 4;
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'T', b'C', b'G')] = 4;

        // context CCT (codon before stop codon TAG)
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'A', b'C', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'G', b'C', b'T')] = 1;
        ans[mutation_type_index(MutationClass::Synonymous, b'C', b'T', b'C', b'T')] = 1;
        
        // missense mutations
        
        // position 1 mutations in codon CCC: 5 codons with context CCC
        // (second codon of each of 5 CDS's)
        // position 2 mutations in codon CCC: 10 CCC codons with context CCC
        ans[mutation_type_index(MutationClass::Missense, b'C', b'A', b'C', b'C')] = 15;
        ans[mutation_type_index(MutationClass::Missense, b'C', b'G', b'C', b'C')] = 15;
        ans[mutation_type_index(MutationClass::Missense, b'C', b'T', b'C', b'C')] = 15;

        // position 1 mutation in codon CCC with context GCC (first codon of each of CDS 1-4)
        ans[mutation_type_index(MutationClass::Missense, b'C', b'A', b'G', b'C')] = 4;
        ans[mutation_type_index(MutationClass::Missense, b'C', b'G', b'G', b'C')] = 4;
        ans[mutation_type_index(MutationClass::Missense, b'C', b'T', b'G', b'C')] = 4;

        // position 1 mutation in codon CCC with context ACC (first codon of CDS 5)
        ans[mutation_type_index(MutationClass::Missense, b'C', b'A', b'A', b'C')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'C', b'G', b'A', b'C')] = 1;
        ans[mutation_type_index(MutationClass::Missense, b'C', b'T', b'A', b'C')] = 1;

        // discard length information from type s.t. compiler does not complain
        let pout: &[u32] = &out;
        let pans: &[u32] = &ans;

        assert_eq!(pout, pans);
    }

}

