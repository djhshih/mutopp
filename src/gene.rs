pub use linked_hash_map::LinkedHashMap;
pub use bio::utils::Strand;

pub type Pos = u64;

/// Gene.
/// All positions are 0-based.
#[derive(Debug)]
pub struct Gene {
    /// Gene name
    pub name: String,
    /// Chromosome or contig name
    pub chrom: String,
    /// Genomic start position
    pub start: Pos,
    /// Genomic end position (exclusive)
    pub end: Pos,
    /// Genomic strand
    pub strand: Strand,
    /// Collection of transcripts
    pub transcripts: LinkedHashMap<String, Transcript>,
}

/// Transcript.
/// All positions are 0-based.
#[derive(Debug)]
pub struct Transcript {
    /// Genomic start position
    pub start: Pos,
    /// Genomic end position (exclusive)
    pub end: Pos,
    /// Disjoint coding regions in genomic coodinates but sorted 5' to 3' of the gene
    pub coding_regions: Vec<Region>,
}

/// Region.
/// All positions are 0-based.
#[derive(Debug)]
pub struct Region {
    /// Number of nucleotides to remove before the start of the first complete codon
    pub phase: u8,
    /// Genomic start position
    pub start: Pos,
    /// Genomic end position (exclusive)
    pub end: Pos,
}