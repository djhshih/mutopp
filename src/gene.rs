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

impl Gene {

    /// Return the transcript with the longest coding regions.
    pub fn canonical_transcript(&self) -> Option<&Transcript> {
        let mut canonical = None;
        let mut canonical_len = 0;
        for x in self.transcripts.iter() {
            let t = x.1;
            let current_len = t.coding_len();
            if current_len > canonical_len {
                canonical = Some(t);
                canonical_len = current_len;
            }
        }
        canonical
    }
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

impl Transcript {
    #[inline]
    pub fn len(&self) -> Pos {
        self.end - self.start
    }

    pub fn coding_len(&self) -> Pos {
        let mut n = 0;
        for r in &self.coding_regions {
            n += r.len(); 
        }
        n
    }
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

impl Region {
    #[inline]
    pub fn len(&self) -> Pos {
        self.end - self.start
    }
}