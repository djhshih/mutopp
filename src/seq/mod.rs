use std::str;

use bio::utils::Strand;

use gene::{Gene,Transcript};
use mutation::{MutImpact, MutOpps, MutEffect};
use constants::*;

pub type Nucleotide = u8;
pub type Residue = u8;

pub type Codon = [Nucleotide; 3];

pub type DnaSeq = Vec<Nucleotide>;
pub type Peptide = Vec<Residue>;

pub type DnaSeqs = Vec<DnaSeq>;

pub fn complement(x: Nucleotide) -> Nucleotide {
    match x {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'X',
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

pub fn print_cds(x: &[Nucleotide]) {
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

pub fn translate(dna: &[Nucleotide]) -> Peptide {
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
            let codon: Codon = [dna[i], dna[i+1], dna[i+2]];
            protein.push(codon_to_aa(&codon));
        }
        i += 3;
    }

    protein
}

fn print_translation(dna: &[Nucleotide]) {
    let protein = translate(dna);
    let s = String::from_utf8(protein).expect("Found invalid protein sequence");
    println!("{}", s);
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
fn codon_to_aa(codon: &Codon) -> Residue {
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
        _ => b'.',
    }
}

pub struct CodingSequence {
    /// Coding sequences, each of which is contiguous in the genome
    /// If a feature is on the reverse strand, the sequence should have already been reverse-complemented.
    pub seqs: Vec<DnaSeq>,
    /// Number of nucleotides padded to the 5' and 3' 
    pub padding: usize,
}

impl CodingSequence {
    /// Assign type to observed mutation.
    /// 
    /// Assign type index to observed mutation at `cds_pos` with a 
    /// possibly non-zero `offset` for mutations close to coding region 
    /// (e.g. splicing site). Coding position 0 is the first nucleotide of 
    /// the start codon. The query position, as well as both the 5' and 3'
    /// positions, must be within the available padded coding sequence.
    /// The `offset` cannot shift the query position outside the padded coding
    /// sequence region containing `cds_pos`. 
    pub fn assign_mutation(&self, gene: &Gene, transcript: &Transcript, g_pos: u64, g_nt_ref: Nucleotide, g_nt_alt: Nucleotide) -> Option<MutEffect> {
        let padding = self.padding;
        let seqs = &self.seqs;
        let strand = gene.strand;
        let regions = &transcript.coding_regions;
        
        let splice_site_len: u64 = 2;
        
        // find CDS region i (including splice site) that contains g_pos
        // convert genomic position to CDS position with possibly non-zero offset
        let mut c_pos: u64 = 0;
        let mut offset: i64 = 0;
        let nt_ref: Nucleotide;
        let nt_alt: Nucleotide;
        match strand {
            Strand::Forward => {
                nt_ref = g_nt_ref;
                nt_alt = g_nt_alt;
                let mut i = 0;
                for region in regions {
                    // pad region to include splice sites,
                    // but splice site padding do not apply to first and last CDS
                    let start_pad: u64 = if i == 0 { 0 } else { splice_site_len };
                    let end_pad: u64 = if i == regions.len() - 1 { 0 } else { splice_site_len };
                    if g_pos < region.end + end_pad {
                        let left = if region.start > start_pad {
                            region.start - start_pad
                        } else { 0 };
                        if g_pos >= left {
                            // found region
                            if g_pos < region.start {
                                // acceptor splice site (5' of CDS)
                                // leave c_pos at first nucleotide of current CDS but update offset
                                offset = -((region.start - g_pos) as i64);
                            } else if g_pos >= region.end {
                                // donor splice site (3' of CDS)
                                // set c_pos to last nucleotide of current CDS and update offset
                                c_pos += region.end - region.start;
                                offset = (g_pos - (region.end - 1)) as i64;
                            } else {
                                // in coding region
                                c_pos += g_pos - region.start;
                                offset = 0;
                            }
                            break;
                        } else {
                            // since regions are sorted, g_pos is not in a region
                            return None;
                        }
                    }
                    c_pos += region.end - region.start;
                    i += 1;
                }
            },
            Strand::Reverse => {
                // complement the nucleotides to convert them from genomic space to coding space
                nt_ref = complement(g_nt_ref);
                nt_alt = complement(g_nt_alt);
                let mut i = 0;
                for region in regions {
                    // pad region to include splice sites,
                    // but splice site padding do not apply to first and last CDS
                    let start_pad: u64 = if i == 0 { 0 } else { splice_site_len };
                    let end_pad: u64 = if i == regions.len() - 1 { 0 } else { splice_site_len };
                    // regions are in genomic coordinates but listed in 5' to 3' order of the gene
                    // therefore, the start coordinate of each CDS correspond to the end of the region,
                    // and the end coordinate correspond to the start of the region
                    // futher, regions are in [region.start, region.end)
                    let right = region.end + start_pad;
                    if g_pos < right {
                        let left = if region.start > end_pad {
                            region.start - end_pad
                        } else { 0 };
                        if g_pos >= left {
                            // found region
                            if g_pos < region.start {
                                // donor splice site (3' of CDS)
                                // set c_pos to last nucleotide of current CDS and update offset
                                c_pos += region.end - region.start;
                                offset = (region.start - g_pos) as i64;
                            } else if g_pos >= region.end {
                                // acceptor splice site (5' of CDS)
                                // leave c_pos at first nucleotide of current CDS but update offset
                                offset = -((g_pos - (region.end - 1)) as i64);
                            } else {
                                // in coding region
                                c_pos += (region.end - 1) - g_pos;
                                offset = 0;
                            }
                            break;
                        } else {
                            // since regions are sorted, g_pos is not in a region
                            return None;
                        }
                    }
                    c_pos += region.end - region.start;
                    i += 1;
                }
            },
            Strand::Unknown => {
                // cannot classify mutation if gene strand is unknown
                return None;
            },
        }
        
        // TODO merge this into previous code block
        // find segment i containing coding_pos
        // cds_begin and cds_end are coordinates in coding space
        let cds_pos = c_pos as usize;
        let mut cds_end: usize = 0;
        let mut cds_len: usize = 0;
        let mut i = 0;
        // number of nucleotides to remove from the current CDS until the next codon
        let mut phase = 0;
        for seq in seqs {
            cds_len = seq.len() - (2 * padding);
            cds_end += cds_len;
            if cds_pos < cds_end {
                break;
            }
            i += 1;
            // compute phase for the next CDS
            phase = (3 - (cds_len % 3)) % 3;
        }
        let cds_begin = cds_end - cds_len;
        // TODO check if calculated phase matches annotated phase
        //let phase = self.regions[i].phase;
        //assert!(phase >= 0 && phase <= 2);
        let seq = &seqs[i];
        
        // compute local index of the query position within CDS i
        let local_idx = (cds_pos - cds_begin + padding) as i64 + offset;
        
        // Check if mutation occur outside the bounds of legal region
        
        // mutation cannot occur at (or before) first local index (0), because no 5' context is available
        // mutation cannot occur at (or after) last local index (n - 1), because no 3' context is available
        if local_idx <= 0 || local_idx >= seq.len() as i64 - 1 {
            // NB  When seq.len() == 1, the mutation will not be assigned a type.
            //     This can happen if the last nuleoide of the first exon is first position of the start codon
            //     and the remaining two nucleotides of the start codon occur in the next exon.
            //     But `padding` must have been set to 0, which provides no context to classify mutation.
            //     If `padding > 1, then seq.len() >= 3, and this edge case will be handled later.
            //     Therefore, it is appropriate to not classify mutation here.
            return None;
        }
        
        let local_idx = local_idx as usize;
        // check that mutation reference allele matches reference sequence
        if seq[local_idx] != nt_ref {
            return None;
        }
        // previous local index bound check ensures that nt_5p and nt_3p are well defined
        let nt_5p = seq[local_idx - 1];
        let nt_3p = seq[local_idx + 1];
        
        // Check whether mutation occurs at a splice site
        let splice_site_len = splice_site_len as usize;
        
        // padding must be at least 3 for splice site mutation to be considered
        if padding >= 3 {
            // acceptor splice site is the last two nucleotides of each intron,
            // i.e. the two intronic nucleotides 5' of a CDS (but not the first CDS)
            if (i > 0 && local_idx >= padding - splice_site_len && local_idx < padding) ||
            // donor splice site is the first two nucleotides of each intron,
            // i.e. the two intronic nucleotides 3' of a CDS (but not the last CDS)
            (i < seqs.len() - 1 && local_idx >= padding + cds_len && local_idx < padding + cds_len + splice_site_len) {
                let cl = MutImpact::SpliceSite;
                return Some( MutEffect {
                    c_pos: cds_pos as u64,
                    c_offset: offset as i32,
                    nt_ref: nt_ref,
                    nt_alt: nt_alt,
                    aa_ref: b'.',
                    aa_alt: b'.',
                    impact: cl,
                    type_index: MutOpps::index(cl, nt_ref, nt_alt, nt_5p, nt_3p),
                } );
            }
        }
        
        // Check whether mutation occurs in coding region
        
        let local_cds_begin = padding;
        let local_cds_end = padding + cds_len;
        if local_idx >= local_cds_begin && local_idx < local_cds_end {
            // find the codon that encompasses the query position; henceforth, query codon
            let codon_ref: [Nucleotide; 3];
            // codon position of the query {0, 1, 2}
            let codon_pos: usize;
            if local_idx < phase {
                // query codon begins in the previous CDS
                if i == 0 {
                    // Invalid input arguments: first CDS should begin with phase = 0
                    return None;
                }
                let seq_prev = &seqs[i-1];
                let seq_prev_cds_end = seq_prev.len() - padding;
                let seq_prev_cds_len = seq_prev_cds_end - padding;
                if seq_prev_cds_len < phase {
                    // Invalid input arguments:
                    // there is not enough CDS nucleotides in the previous CDS
                    return None;
                }
                codon_ref = match phase {
                    1 => [seq_prev[seq_prev_cds_end - 2], seq_prev[seq_prev_cds_end - 1], seq[local_cds_begin]],
                    2 => [seq_prev[seq_prev_cds_end - 1], seq[local_cds_begin], seq[local_cds_begin+1]],
                    _ => panic!("Invalid phase value"),
                };
                codon_pos = 3 - phase;
            } else {
                // query codon begins in current CDS
                codon_pos = (local_idx - phase) % 3;
                
                // local index of the start of the query codon
                let local_codon_begin = local_idx - codon_pos;
                
                if local_codon_begin + 2 >= local_cds_end - 1 {
                    // codon extends into the next CDS
                    if i == seqs.len() - 1 {
                        // Invalid input arguments: codon is incomplete and no more CDS available
                        return None;
                    }
                    let seq_next = &seqs[i+1];
                    if local_codon_begin + 1 < local_cds_end {
                        // codon extends 1 nt into the next CDS
                        codon_ref = [seq[local_codon_begin], seq[local_codon_begin+1], seq_next[padding]];
                    } else {
                        // codon extends 2 nt into the next CDS
                        codon_ref = [seq[local_codon_begin], seq_next[padding], seq_next[padding+1]];
                    }
                } else {
                    // codon is within the current CDS
                    codon_ref = [seq[local_codon_begin], seq[local_codon_begin+1], seq[local_codon_begin+2]];
                }
            }
            
            let context = [nt_5p, nt_ref, nt_3p];
            
            let aa_ref = if i == 0 && local_idx == local_cds_begin {
                // first codon of the first CDS is the start codon
                // this is usually ATG but alternative start codons are possible
                b'^'
            } else {
                codon_to_aa(&codon_ref)
            };
            
            let codon_alt = mutate_codon(&codon_ref, codon_pos, nt_alt);
            let aa_alt = codon_to_aa(&codon_alt);
            let cl = mutation_impact(aa_ref, aa_alt);
            return Some( MutEffect {
                c_pos: cds_pos as u64,
                c_offset: offset as i32,
                nt_ref: nt_ref,
                nt_alt: nt_alt,
                aa_ref: aa_ref,
                aa_alt: aa_alt,
                impact: cl,
                type_index: MutOpps::index(cl, nt_ref, nt_alt, nt_5p, nt_3p),
            } );
        }

        // mutation cannot be assigned a type
        return None;
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
    pub fn count_opp(&self) -> Option<MutOpps> {
        if self.seqs.is_empty() {
            return None;
        }

        let seqs = &self.seqs;
        let padding = self.padding;
        let nseqs = seqs.len();

        let mut x = MutOpps::new();
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
            while i + 2 < cds_end {
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
                        // there is probably a problem with the sequence...
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
            count_opp_at_splice_sites(&mut x, self);
        }

        if error {
            None
        } else {
            Some(x)
        }
    }
    
    pub fn is_valid(&self) -> bool {
        let seqs = &self.seqs;
        let padding = self.padding;
        
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

        // check for presence of stop codon at the final position using the standard genetic code
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
    
/*
    fn prev_nt(cds_pos: usize) {

    }
    fn next_nt(cds_pos: usize) {

    }
    fn codon_at(cds_pos: usize) {

    }
*/
}

fn count_opp_at_splice_sites(mut x: &mut MutOpps, cds: &CodingSequence) {
    let padding = cds.padding;
    if padding >= 3 {
        let nseqs = cds.seqs.len();
        for s in 0 .. nseqs {
            let seq = &cds.seqs[s];

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

fn count_opp_at_slice_site(x: &mut MutOpps, seq: &[u8], begin: usize) {
    let splice =  MutImpact::SpliceSite;
    // iterate through the 2 positions of the splice site
    for j in 0 .. 2 {
        let context_5p = seq[begin-1+j];
        let nt_ref = seq[begin+j];
        let context_3p = seq[begin+1+j];
        for &nt_alt in nucleotides.iter() {
            if nt_alt != nt_ref {
                let idx = MutOpps::index(splice, nt_ref, nt_alt, context_5p, context_3p);
                x[idx] += 1;
            }
        }
    }
}

/// Count mutation opportunities at codon position.
/// x: opportunity matrix
/// pos: 0, 1, 2 within the codon
/// codon: codon sequence
/// context: sequence centered around the query position
fn count_opp_at_codon_pos(x: &mut MutOpps, pos: usize, codon: &[u8; 3], aa_ref: u8, context: &[u8; 3]) {
    let nt_ref = context[1];
    // iterate through possible substitutions
    for &nt_alt in nucleotides.iter() {
        if nt_alt != nt_ref {
            let codon_alt = mutate_codon(&codon, pos, nt_alt);
            let aa_alt = codon_to_aa(&codon_alt);
            let cl = mutation_impact(aa_ref, aa_alt);
            let idx = MutOpps::index(cl, nt_ref, nt_alt, context[0], context[2]);
            x[idx] += 1;
        }
    }
}

/// Assign mutation impact for observed mutation at codon position.
fn mutation_impact(aa_ref: u8, aa_alt: u8) -> MutImpact {
    if aa_ref == b'^' {
        // assume that any mutation at start codon causes start codon loss
        // this is true for ATG start codon (no degeneracy) and
        // probably true for alternative start codons as well
        MutImpact::StartOrStop
    } else {
        if aa_alt == aa_ref {
            MutImpact::Synonymous 
        } else if aa_ref == b'$' || aa_alt == b'$' {
            MutImpact::StartOrStop
        } else {
            MutImpact::Missense
        }
    }
}

/// Mutate codon.
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
        let test = CodingSequence {
            seqs: vec![
                b"ATGATG".to_vec(),
                b"ATGATGATG".to_vec(),
                b"ATGATGATGATGTAG".to_vec(),
            ],
            padding: 0
        };
        let out = test.count_opp().expect("Invalid opportunity array");
        
        // since no 5' or 3' padding is available, the first and last nucleotides
        // of each CDS have no context, and they are therefore not counted

        let mut ans = MutOpps::new();

        // mutations in position 2 of start codon (ATG) with context ATG
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // mutations in position 3 of start codon (ATG) with context TGA
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'A', b'T', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'C', b'T', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'T', b'T', b'A')] = 1;

        // mutations in TGA context (except involving start codon)
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'A')] = 5;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'A')] = 5;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'A')] = 5;

        // mutations in GAT context
        ans[MutOpps::index(MutImpact::Missense, b'A', b'C', b'G', b'T')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'G', b'G', b'T')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'T', b'G', b'T')] = 6;

        // mutations in ATG context (except involving start codon)
        ans[MutOpps::index(MutImpact::Missense, b'T', b'A', b'A', b'G')] = 8;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'C', b'A', b'G')] = 8;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'G', b'A', b'G')] = 8;

        // mutation in position 3 of penultimate codon (ATG) with context TGT
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'T')] = 1;

        // mutation in position 1 of stop codon (TAG) with context GTA
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'G', b'A')] = 1;

        // mutation in position 2 of stop codon (TAG) with context TAG
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        assert_eq!(&out, &ans);
    }

    #[test]
    fn test_count_opp_pad1() {
        // sequences contain CDS with 1 nucleotide padding
        let test = CodingSequence {
            seqs: vec![
                b"CATGATGC".to_vec(),
                b"CATGATGATGC".to_vec(),
                b"CATGATGATGATGTAGC".to_vec(),
            ],
            padding: 1
        };
        let out = test.count_opp().expect("Invalid opportunity array");
        
        // the result should be similar to the first test,
        // except now that first and last nucleotides are each CDS have
        // 5' and 3' context and are thus counted

        let mut ans = MutOpps::new();

        // mutations in position 2 of start codon (ATG) with context ATG
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // mutations in position 3 of start codon (ATG) with context TGA
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'A', b'T', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'C', b'T', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'T', b'T', b'A')] = 1;

        // mutations in TGA context (except involving start codon)
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'A')] = 5;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'A')] = 5;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'A')] = 5;

        // mutations in GAT context
        ans[MutOpps::index(MutImpact::Missense, b'A', b'C', b'G', b'T')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'G', b'G', b'T')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'T', b'G', b'T')] = 6;

        // mutations in ATG context (except involving start codon)
        ans[MutOpps::index(MutImpact::Missense, b'T', b'A', b'A', b'G')] = 8;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'C', b'A', b'G')] = 8;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'G', b'A', b'G')] = 8;

        // mutation in position 3 of penultimate codon (ATG) with context TGT
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'T')] = 1;

        // mutation in position 1 of stop codon (TAG) with context GTA
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'G', b'A')] = 1;

        // mutation in position 2 of stop codon (TAG) with context TAG
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        // additional mutations counted in pad1 but not in pad0

        // mutaiton in first nucleotide of start codon with context CAT
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'C', b'C', b'T')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'G', b'C', b'T')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'T', b'C', b'T')] = 1;

        // mutaiton in first nucleotide of CDS 2 and 3 with context CAT
        ans[MutOpps::index(MutImpact::Missense, b'A', b'C', b'C', b'T')] = 2;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'G', b'C', b'T')] = 2;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'T', b'C', b'T')] = 2;

        // mutation last nucleotide of CDS 1 and 2 with context TGC
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'C')] = 2;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'C')] = 2;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'C')] = 2;

        // mutation in position 3 of stop codon (TAG -> TAA) with context AGC
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'C', b'A', b'C')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'T', b'A', b'C')] = 1;

        // synonymous mutation in position 3 of stop codon (TAG -> TAA) with context AGC
        ans[MutOpps::index(MutImpact::Synonymous, b'G', b'A', b'A', b'C')] = 1;

        assert_eq!(&out, &ans);
    }

    #[test]
    fn test_count_opp_pad0_split() {
        // sequences contain only CDS without any padding but with split codons
        let test = CodingSequence {
            seqs: vec![
                b"ATGATGA".to_vec(),
                b"TGATGATGAT".to_vec(),
                b"GATGATGATGTAG".to_vec(),
            ],
            padding: 0
        };
        let out = test.count_opp().expect("Invalid opportunity array");
        
        // since no 5' or 3' padding is available, the first and last nucleotides
        // of each CDS have no context, and they are therefore not counted

        let mut ans = MutOpps::new();

        // mutations in position 2 of start codon (ATG) with context ATG
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // mutations in position 3 of start codon (ATG) with context TGA
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'A', b'T', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'C', b'T', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'T', b'T', b'A')] = 1;

        // mutations in TGA context (except involving start codon)
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'A')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'A')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'A')] = 6;

        // mutations in GAT context
        ans[MutOpps::index(MutImpact::Missense, b'A', b'C', b'G', b'T')] = 7;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'G', b'G', b'T')] = 7;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'T', b'G', b'T')] = 7;

        // mutations in ATG context (except involving start codon)
        ans[MutOpps::index(MutImpact::Missense, b'T', b'A', b'A', b'G')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'C', b'A', b'G')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'G', b'A', b'G')] = 6;

        // mutation in position 3 of penultimate codon (ATG) with context TGT
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'T')] = 1;

        // mutation in position 1 of stop codon (TAG) with context GTA
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'G', b'A')] = 1;

        // mutation in position 2 of stop codon (TAG) with context TAG
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        assert_eq!(&out, &ans);
    }

    #[test]
    fn test_count_opp_pad1_split() {
        // sequences contain CDS with 1 nucleotide padding and with split codons
        let test = CodingSequence {
            seqs: vec![
                b"CATGATGAC".to_vec(),
                b"CTGATGATGATC".to_vec(),
                b"CGATGATGATGTAGC".to_vec(),
            ],
            padding: 1
        };
        let out = test.count_opp().expect("Invalid opportunity array");
        
        // the result should be similar to the first test,
        // except now that first and last nucleotides are each CDS have
        // 5' and 3' context and are thus counted

        let mut ans = MutOpps::new();

        // mutations in position 2 of start codon (ATG) with context ATG
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // mutations in position 3 of start codon (ATG) with context TGA
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'A', b'T', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'C', b'T', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'T', b'T', b'A')] = 1;

        // mutations in TGA context (except involving start codon)
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'A')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'A')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'A')] = 6;

        // mutations in GAT context
        ans[MutOpps::index(MutImpact::Missense, b'A', b'C', b'G', b'T')] = 7;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'G', b'G', b'T')] = 7;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'T', b'G', b'T')] = 7;

        // mutations in ATG context (except involving start codon)
        ans[MutOpps::index(MutImpact::Missense, b'T', b'A', b'A', b'G')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'C', b'A', b'G')] = 6;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'G', b'A', b'G')] = 6;

        // mutation in position 3 of penultimate codon (ATG) with context TGT
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'T', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'T', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'T', b'T')] = 1;

        // mutation in position 1 of stop codon (TAG) with context GTA
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'G', b'A')] = 1;

        // mutation in position 2 of stop codon (TAG) with context TAG
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        // additional mutations counted in pad1 but not in pad0

        // mutaiton in first nucleotide of start codon with context CAT
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'C', b'C', b'T')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'G', b'C', b'T')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'T', b'C', b'T')] = 1;

        // mutaiton in first nucleotide of CDS 2 with context CTG
        ans[MutOpps::index(MutImpact::Missense, b'T', b'A', b'C', b'G')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'C', b'C', b'G')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'G', b'C', b'G')] = 1;

        // mutaiton in first nucleotide of CDS 3 with context CGA
        ans[MutOpps::index(MutImpact::Missense, b'G', b'A', b'C', b'A')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'C', b'C', b'A')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'G', b'T', b'C', b'A')] = 1;

        // mutation last nucleotide of CDS 1 with context GAC
        ans[MutOpps::index(MutImpact::Missense, b'A', b'C', b'G', b'C')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'G', b'G', b'C')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'A', b'T', b'G', b'C')] = 1;

        // mutation last nucleotide of CDS 2 with context ATC
        ans[MutOpps::index(MutImpact::Missense, b'T', b'A', b'A', b'C')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'C', b'A', b'C')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'T', b'G', b'A', b'C')] = 1;

        // mutation in position 3 of stop codon (TAG -> TAA) with context AGC
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'C', b'A', b'C')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'T', b'A', b'C')] = 1;

        // synonymous mutation in position 3 of stop codon (TAG -> TAA) with context AGC
        ans[MutOpps::index(MutImpact::Synonymous, b'G', b'A', b'A', b'C')] = 1;

        assert_eq!(&out, &ans);
    }

    #[test]
    fn test_count_opp_splice() {
        // sequences contain CDS with 3 nucleotide padding,
        // permitting splice site mutation opportunity counting
        let test = CodingSequence {
            seqs: vec![
                b"ACCATGCCCCCCGTA".to_vec(),
                b"CAGCCCCCCGTG".to_vec(),
                b"TAGCCCCCCGTA".to_vec(),
                b"AAGCCCCCCGTG".to_vec(),
                b"GCACCCCCCTAGCCC".to_vec(),
            ],
            padding: 3
        };
        let out = test.count_opp().expect("Invalid opportunity array");

        let mut ans = MutOpps::new();

        // donor splice site mutations

        // two cGTa and two cGTg donor splice sites, mutation at position 1  with context CGT
        ans[MutOpps::index(MutImpact::SpliceSite, b'G', b'A', b'C', b'T')] = 4;
        ans[MutOpps::index(MutImpact::SpliceSite, b'G', b'C', b'C', b'T')] = 4;
        ans[MutOpps::index(MutImpact::SpliceSite, b'G', b'T', b'C', b'T')] = 4;

        // two cGTa donor splice sites, mutation at position 2 with context GTA
        ans[MutOpps::index(MutImpact::SpliceSite, b'T', b'A', b'G', b'A')] = 2;
        ans[MutOpps::index(MutImpact::SpliceSite, b'T', b'C', b'G', b'A')] = 2;
        ans[MutOpps::index(MutImpact::SpliceSite, b'T', b'G', b'G', b'A')] = 2;

        // two cGTg donor splice sites, mutation at position 2 with context GTG
        ans[MutOpps::index(MutImpact::SpliceSite, b'T', b'A', b'G', b'G')] = 2;
        ans[MutOpps::index(MutImpact::SpliceSite, b'T', b'C', b'G', b'G')] = 2;
        ans[MutOpps::index(MutImpact::SpliceSite, b'T', b'G', b'G', b'G')] = 2;

        // acceptor splice site mutations

        // cAGc acceptor site, mutation at position 1 with context CAG
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'C', b'C', b'G')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'T', b'C', b'G')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'G', b'C', b'G')] = 1;

        // tAGc acceptor site, mutation at position 1 with context TAG
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'C', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'T', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'G', b'T', b'G')] = 1;

        // aAGc acceptor site, mutation at position 1 with context AAG
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'C', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'T', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'G', b'A', b'G')] = 1;

        // cAGc, tAGc, and aAGc acceptor sites, mutation at position 2 with context AGC
        ans[MutOpps::index(MutImpact::SpliceSite, b'G', b'A', b'A', b'C')] = 3;
        ans[MutOpps::index(MutImpact::SpliceSite, b'G', b'C', b'A', b'C')] = 3;
        ans[MutOpps::index(MutImpact::SpliceSite, b'G', b'T', b'A', b'C')] = 3;

        // gCAc acceptor site, mutation at position 1 with context GCA
        ans[MutOpps::index(MutImpact::SpliceSite, b'C', b'A', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'C', b'T', b'G', b'A')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'C', b'G', b'G', b'A')] = 1;

        // gCAc acceptor site, mutation at position 2 with context CAC
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'C', b'C', b'C')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'G', b'C', b'C')] = 1;
        ans[MutOpps::index(MutImpact::SpliceSite, b'A', b'T', b'C', b'C')] = 1;

        // start codon mutations

        // start codon mutations, position 1 with context CAT
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'C', b'C', b'T')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'G', b'C', b'T')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'T', b'C', b'T')] = 1;

        // start codon mutations, position 2 with context ATG
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'A', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'A', b'G')] = 1;

        // start codon mutations, position 3 with context TGC
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'A', b'T', b'C')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'C', b'T', b'C')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'T', b'T', b'C')] = 1;

        // stop codon mutations, position 1 with context CTA
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'A', b'C', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'C', b'C', b'A')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'T', b'G', b'C', b'A')] = 1;

        // stop codon mutations, position 2 with context TAG
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'C', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'G', b'T', b'G')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'A', b'T', b'T', b'G')] = 1;

        // stop codon mutations, position 3 with context AGC
        // TAG -> TAA is synonymous stop codon mutation
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'C', b'A', b'C')] = 1;
        ans[MutOpps::index(MutImpact::StartOrStop, b'G', b'T', b'A', b'C')] = 1;

        // synonymous mutations
        
        // stop codon mutations, position 3 with context AGC
        // TAG -> TAA is synonymous stop codon mutation
        ans[MutOpps::index(MutImpact::Synonymous, b'G', b'A', b'A', b'C')] = 1;
        
        // position 3 mutations in codon CCC are synonymous
        // CCC -> CCA, CCG, CCT (all code for Proline)

        // context CCC (first CCC codon of ech of 5 CDS's)
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'A', b'C', b'C')] = 5;
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'G', b'C', b'C')] = 5;
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'T', b'C', b'C')] = 5;
        
        // context CCG (last codon of each of 4 CDS's)
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'A', b'C', b'G')] = 4;
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'G', b'C', b'G')] = 4;
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'T', b'C', b'G')] = 4;

        // context CCT (codon before stop codon TAG)
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'A', b'C', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'G', b'C', b'T')] = 1;
        ans[MutOpps::index(MutImpact::Synonymous, b'C', b'T', b'C', b'T')] = 1;
        
        // missense mutations
        
        // position 1 mutations in codon CCC: 5 codons with context CCC
        // (second codon of each of 5 CDS's)
        // position 2 mutations in codon CCC: 10 CCC codons with context CCC
        ans[MutOpps::index(MutImpact::Missense, b'C', b'A', b'C', b'C')] = 15;
        ans[MutOpps::index(MutImpact::Missense, b'C', b'G', b'C', b'C')] = 15;
        ans[MutOpps::index(MutImpact::Missense, b'C', b'T', b'C', b'C')] = 15;

        // position 1 mutation in codon CCC with context GCC (first codon of each of CDS 1-4)
        ans[MutOpps::index(MutImpact::Missense, b'C', b'A', b'G', b'C')] = 4;
        ans[MutOpps::index(MutImpact::Missense, b'C', b'G', b'G', b'C')] = 4;
        ans[MutOpps::index(MutImpact::Missense, b'C', b'T', b'G', b'C')] = 4;

        // position 1 mutation in codon CCC with context ACC (first codon of CDS 5)
        ans[MutOpps::index(MutImpact::Missense, b'C', b'A', b'A', b'C')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'C', b'G', b'A', b'C')] = 1;
        ans[MutOpps::index(MutImpact::Missense, b'C', b'T', b'A', b'C')] = 1;

        assert_eq!(&out, &ans);
    }

    #[test]
    fn test_reverse_complement() {
        {
            let mut x = b"ATAGAC".to_vec();
            reverse_complement(&mut x);
            assert_eq!(b"GTCTAT".to_vec(), x);
        }
        {
            let mut x = b"GTATCGAAC".to_vec();
            reverse_complement(&mut x);
            assert_eq!(b"GTTCGATAC".to_vec(), x);
        }
        {
            let mut x = b"TTAACTTXTTATGC".to_vec();
            reverse_complement(&mut x);
            assert_eq!(b"GCATAAXAAGTTAA".to_vec(), x);
        }
    }

}

