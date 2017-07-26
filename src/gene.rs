pub use std::collections::HashMap;
pub use bio::utils::Strand;

/// Gene.
/// All positions are 0-based.
#[derive(Debug)]
pub struct Gene {
    /// Gene name
    pub name: String,
    /// Chromosome or contig name
    pub chrom: String,
    /// Genomic start position
    pub start: u64,
    /// Genomic end position (exclusive)
    pub end: u64,
    /// Genomic strand
    pub strand: Strand,
    /// Collection of transcripts
    pub transcripts: HashMap<String, Transcript>,
}

/// Transcript.
/// All positions are 0-based.
#[derive(Debug)]
pub struct Transcript {
    /// Genomic start position
    pub start: u64,
    /// Genomic end position (exclusive)
    pub end: u64,
    /// Disjoint coding regions sorted by position
    pub coding_regions: Vec<Region>,
}

/// Region.
#[derive(Debug)]
pub struct Region {
    /// Number of nucleotides to remove before the start of the first complete codon
    pub phase: u8,
    /// Genomic start position
    pub start: u64,
    /// Genomic end position (exclusive)
    pub end: u64,
}