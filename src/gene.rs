pub use std::collections::HashMap;
pub use bio::utils::Strand;

#[derive(Debug)]
pub struct Gene {
    pub name: String,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub transcripts: HashMap<String, Transcript>,
}

#[derive(Debug)]
pub struct Transcript {
    pub start: u64,
    pub end: u64,
    /// Disjoint coding regions sorted by position
    pub coding_regions: Vec<Region>,
}

#[derive(Debug)]
pub struct Region {
    pub phase: i8,
    pub start: u64,
    pub end: u64,
}