use std::fmt;
use std::ops;
use std::ascii::AsciiExt;
use std::slice;

use gene::Pos;
use seq::{Nucleotide,Residue};
use constants::*;

#[derive(Copy, Clone)]
pub enum SeqOntology {
    //TranscriptAblation,
    SpliceAcceptor,
    SpliceDonor,
    StopGained,
    //Frameshift,
    StopLost,
    StartLost,
    //TranscriptAmplification,
    //InframeInsertion,
    //InframeDeletion,
    Missense,
    //?ProteinAltering,
    //?SpliceRegion,
    //IncompleteTerminalCodon,
    StopRetained,
    Synonymous,
    //?CodingSequence,
    //MatureMiRNA,
    FivePrimeUTR,
    ThreePrimeUTR,
    //NoncodingTranscriptExon,
    Intron,
    //NMDTranscript,
    //NoncodingTranscript,
    Upstream,
    Downstream,
    //TFBSAblation,
    //TFBSAmplification,
    //TFBindingSite,
    //RegulatoryRegionAblation,
    //RegulatoryRegionAmplification,
    //FeatureElongation,
    //RegulatoryRegion,
    //FeatureTruncation,
    Intergenic,
}

impl fmt::Display for SeqOntology {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let c = match *self {
            //SeqOntology::TranscriptAblation => "transcript_ablation",
            SeqOntology::SpliceAcceptor => "splice_acceptor_variant",
            SeqOntology::SpliceDonor => "splice_donor_variant",
            SeqOntology::StopGained => "stop_gained",
            //SeqOntology::Frameshift => "frameshift_variant",
            SeqOntology::StopLost => "stop_lost",
            SeqOntology::StartLost => "start_lost",
            //SeqOntology::TranscriptAmplification => "transcript_amplification",
            //SeqOntology::InframeInsertion => "inframe_insertion",
            //SeqOntology::InframeDeletion => "inframe_deletion",
            SeqOntology::Missense => "missense_variant",
            //SeqOntology::?ProteinAltering => "protein_altering_variant",
            //SeqOntology::?SpliceRegion => "splice_region_variant",
            //SeqOntology::IncompleteTerminalCodon => "incomplete_terminal_codon_variant",
            SeqOntology::StopRetained => "stop_retained_variant",
            SeqOntology::Synonymous => "synonymous_variant",
            //SeqOntology::?CodingSequence => "coding_sequence_variant",
            //SeqOntology::MatureMiRNA => "mature_mirRNA_variant",
            SeqOntology::FivePrimeUTR => "5_prime_UTR_variant",
            SeqOntology::ThreePrimeUTR => "3_prime_UTR_variant",
            //SeqOntology::NoncodingTranscriptExon => "non_coding_transcript_exon_variant",
            SeqOntology::Intron => "intron_variant",
            //SeqOntology::NMDTranscript => "NMD_transcript_variant",
            //SeqOntology::NoncodingTranscript => "non_coding_transcript_exon_variant",
            SeqOntology::Upstream => "upstream_gene_variant",
            SeqOntology::Downstream => "downstream_gene_variant",
            //SeqOntology::TFBSAblation => "TFBS_ablation",
            //SeqOntology::TFBSAmplification => "TFBS_amplification",
            //SeqOntology::TFBindingSite => "TF_binding_site_variant",
            //SeqOntology::RegulatoryRegionAblation => "regulatory_region_ablation",
            //SeqOntology::RegulatoryRegionAmplification => "regulatory_region_amplification",
            //SeqOntology::FeatureElongation => "feature_elongation",
            //SeqOntology::RegulatoryRegion => "regulatory_region_variant",
            //SeqOntology::FeatureTruncation => "feature_truncation",
            SeqOntology::Intergenic => "intergenic_variant",
        };

        write!(f, "{}", c)
    }
}

#[derive(Copy,Clone)]
pub enum MutImpact {
    Synonymous,
    Missense,
    StartOrStop,
    SpliceSite,
    // non-coding variants are not currently considered
    Noncoding,
}

impl From<SeqOntology> for MutImpact {
    fn from(x: SeqOntology) -> Self {
        match x {
            SpliceAcceptor => MutImpact::SpliceSite,
            SpliceDonor => MutImpact::SpliceSite,
            StopGained => MutImpact::StartOrStop,
            //Frameshift,
            StopLost => MutImpact::StartOrStop,
            StartLost => MutImpact::StartOrStop,
            //TranscriptAmplification,
            //InframeInsertion,
            //InframeDeletion,
            Missense => MutImpact::Missense,
            //?ProteinAltering,
            //?SpliceRegion,
            //IncompleteTerminalCodon,
            StopRetained => MutImpact::Synonymous,
            Synonymous => MutImpact::Synonymous,
            //?CodingSequence,
            //MatureMiRNA,
            FivePrimeUTR => MutImpact::Noncoding,
            ThreePrimeUTR => MutImpact::Noncoding,
            //NoncodingTranscriptExon,
            Intron => MutImpact::Noncoding,
            //NMDTranscript,
            //NoncodingTranscript,
            Upstream => MutImpact::Noncoding,
            Downstream => MutImpact::Noncoding,
            //TFBSAblation,
            //TFBSAmplification,
            //TFBindingSite,
            //RegulatoryRegionAblation,
            //RegulatoryRegionAmplification,
            //FeatureElongation,
            //RegulatoryRegion,
            //FeatureTruncation,
            Intergenic => MutImpact::Noncoding,
        }
    }
}

impl MutImpact {
    pub fn iter() -> slice::Iter<'static, MutImpact> {
        use self::MutImpact::*;
        static classes: [MutImpact; n_mutation_classes] = [Synonymous, Missense, StartOrStop, SpliceSite];
        classes.into_iter()
    }
}

impl fmt::Display for MutImpact {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let c = match *self {
            MutImpact::Synonymous => "syn",
            MutImpact::Missense => "mis",
            MutImpact::StartOrStop => "sos",
            MutImpact::SpliceSite => "spl",
            MutImpact::Noncoding => "ncd",
        };

        write!(f, "{}", c)
    }
}

/// Mutation effect on a transcript from a protein-coding gene.
/// All positions are 0-based
pub struct MutEffect {
    /// HGNC symbol of affected gene
    //gene_symbol: String,
    /// ID of affected transcript
    //transcript_id: String,
    /// Position of mutation in coding DNA sequence
    pub c_pos: Pos,
    /// Offset from the coding sequence position to indicate intronic positions
    pub c_offset: i32,
    /// Reference nucleotide of the coding sequence
    pub nt_ref: Nucleotide,
    /// Alternate nucleotide of the coding sequence
    pub nt_alt: Nucleotide,
    /// Reference amino acid residue
    pub aa_ref: Residue,
    /// Alternate amino acid residue
    pub aa_alt: Residue,
    /// Mutation impact
    pub impact: MutImpact,
    /// Mutation type index
    pub type_index: usize,
}

impl MutEffect {
    #[inline]
    pub fn cdna_pos(&self) -> Pos {
        self.c_pos + 1
    }
    
    #[inline]
    pub fn aa_pos(&self) -> Pos {
        (self.c_pos / 3) + 1
    }
    
    pub fn cdna_change(&self) -> String {
        match self.impact {
            MutImpact::Synonymous | MutImpact::Missense | MutImpact::StartOrStop => {
                format!("c.{}{}>{}", self.cdna_pos(), self.nt_ref as char, self.nt_alt as char)
            },
            MutImpact::SpliceSite => {
                let op = if self.c_offset >= 0 { '+' } else { '-' };
                format!("c.{}{}{}{}>{}", self.cdna_pos(), op, self.c_offset.abs(), self.nt_ref as char, self.nt_alt as char)
            },
            MutImpact::Noncoding => {
                format!("c.?")
            }
        }
    }
    
    pub fn protein_change(&self) -> String {
        match self.impact {
            MutImpact::Synonymous => {
                let aa_ref = if self.aa_ref == b'$' { b'*' } else { self.aa_ref };
                format!("p.{}{}=", aa_ref as char, self.aa_pos())
            },
            MutImpact::Missense => {
                format!("p.{}{}{}", self.aa_ref as char, self.aa_pos(), self.aa_alt as char)
            },
            MutImpact::StartOrStop => {
                if self.aa_ref == b'^' {
                    // start-lost
                    String::from("p.M1?")
                } else if self.aa_ref == b'$' {
                    // stop-lost or no-stop
                    // NB stop-lost and no-stop are currently not distinguished;
                    //    we do not track where the new stop codon would eventually occur
                    format!("p.*{}{}ext*?", self.aa_pos(), self.aa_alt as char)
                } else if self.aa_alt == b'$' {
                    // stop-gained
                    format!("p.{}{}*", self.aa_ref as char, self.aa_pos())
                } else {
                    // NB start-gained is currently not annotated!
                    String::from("p.?")
                }
            },
            MutImpact::SpliceSite => String::from("p.?"),
            MutImpact::Noncoding => String::from("p.?"),
        }
    }
}

pub struct MutOpps([u32; n_mutation_types]);

impl MutOpps {
    #[inline]    
    pub fn iter(&self) -> slice::Iter<u32> {
        self.0.iter()
    }
}

impl PartialEq for MutOpps {
    fn eq(&self, other: &MutOpps) -> bool {
        &self.0 as &[u32] == &other.0 as &[u32]
    }
}

impl fmt::Debug for MutOpps {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(&self.0 as &[u32], f)
    }
}

impl ops::Index<usize> for MutOpps {
    type Output = u32;

    #[inline]
    fn index(&self, index: usize) -> &u32 {
        &self.0[index]
    }
}

impl ops::IndexMut<usize> for MutOpps {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut u32 {
        &mut self.0[index]
    }
}

impl fmt::Display for MutOpps {
    /// Display as n_mutation_channels by n_mutation_classes matrix
    /// matrix is stored in column-major order
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        const I: usize = n_mutation_channels;

        // print header
        try!(write!(f, "channel "));
        for &cl in MutImpact::iter() {
            // since the length of cl is not known at compile time, std::fmt does not
            // pad properly; therefore, pad manually
            try!(write!(f, " {} ", cl));
        }
        try!(writeln!(f, ""));

        let rownames = MutOpps::channels();
        for i in 0 .. I {
            try!(write!(f, "{} ", rownames[i]));
            for &cl in MutImpact::iter() {
                let j = cl as usize;
                try!(write!(f, "{:>4} ", self[j * I + i]));
            }
            try!(writeln!(f, ""));
        }

        // print class totals
        try!(write!(f, "# total "));
        let class_counts: Vec<u32> = MutImpact::iter().map(|&cl| self.sum_opp_in_class(cl)).collect();
        for c in class_counts {
            try!(write!(f, "{:>4} ", c));
        }
        try!(writeln!(f, ""));
        
        Ok(())
    }
}

impl MutOpps {
    #[inline]
    pub fn new() -> MutOpps {
        MutOpps([0; n_mutation_types])
    }
    
    // mutation channels are indexed by factors: impact class, stranded substitution, 5' context, and 3' context
    // each factor serve as a subscript into the mutation channel array
    // consider nucleotides in the order of A, C, G, T
    // stranded substitutions (12): A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T, T>A, T>C, T>G
    // 5' context (4): A, C, G, T
    // 3' context (4): A, C, G, T
    // number of channels = 12 * 4 * 4 = 192
    // later factors vary faster than early factors
    pub fn index(cl: MutImpact, nt_ref: u8, nt_alt: u8, nt_5p: u8, nt_3p: u8) -> usize {
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

    pub fn channels() -> Vec<String> {
        let cl = MutImpact::Synonymous;
        let mut names = vec![String::new(); n_mutation_channels as usize];
        for &nt_ref in nucleotides.iter() {
            for &nt_alt in nucleotides.iter() {
                if nt_ref != nt_alt {
                    for &nt_5p in nucleotides.iter() {
                        for &nt_3p in nucleotides.iter() {
                            let idx = MutOpps::index(cl, nt_ref, nt_alt, nt_5p, nt_3p);
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
    
    pub fn types() -> Vec<String> {
        let mut names = vec![String::new(); n_mutation_types as usize];
        for &cl in MutImpact::iter() {
            for &nt_ref in nucleotides.iter() {
                for &nt_alt in nucleotides.iter() {
                    if nt_ref != nt_alt {
                        for &nt_5p in nucleotides.iter() {
                            for &nt_3p in nucleotides.iter() {
                                let idx = MutOpps::index(cl, nt_ref, nt_alt, nt_5p, nt_3p);
                                names[idx] = format!(
                                    "{class}_{c5}{ref}{c3}_{c5}{alt}{c3}",
                                    class = cl,
                                    ref = (nt_ref as char).to_ascii_lowercase(),
                                    alt = (nt_alt as char).to_ascii_lowercase(),
                                    c5 = (nt_5p as char).to_ascii_lowercase(),
                                    c3 = (nt_3p as char).to_ascii_lowercase(),
                                );
                            }
                        }
                    }
                }
            }
        }
        names
    }
    
    fn sum_opp_in_class(&self, cl: MutImpact) -> u32 {
        const I: usize = n_mutation_channels;
        let j = cl as usize;
        let mut y = 0u32;
        for i in 0 .. I {
            y += self[j * I + i];
        }
        y
    }
}

