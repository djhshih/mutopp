use std::fmt;
use std::ops;
use std::ascii::AsciiExt;
use std::slice;

use constants::*;

/// Genomic mutation opportunites
pub struct MutOpps([u32; N_NONSTRANDED_MUTATION_CHANNELS]);

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

impl MutOpps {
    #[inline]
    pub fn new() -> MutOpps {
        MutOpps([0; N_NONSTRANDED_MUTATION_CHANNELS])
    }
    
    // mutation channels are indexed by factors: impact class, nonstranded substitution, 5' context, and 3' context
    // each factor serve as a subscript into the mutation channel array
    // consider nucleotides in the order of A, C, G, T
    // nonstranded substitutions (6): C>A, C>G, C>T, T>A, T>C, T>G
    //                                G>T, G>C, G>A, A>T, A>G, A>C
    // 5' context (4): A, C, G, T
    // 3' context (4): A, C, G, T
    // number of channels = 6 * 4 * 4 = 96
    // later factors vary faster than early factors
    pub fn index(nt_ref: u8, nt_alt: u8, nt_5p: u8, nt_3p: u8) -> usize {
        assert!(nt_ref != nt_alt);

        let j = match nt_ref {
            b'A' => match nt_alt {
                b'C' => 5,
                b'G' => 4,
                b'T' => 3,
                _ => 0,
            },
            b'C' => match nt_alt {
                b'A' => 0,
                b'G' => 1,
                b'T' => 2,
                _ => 0,
            },
            b'G' => match nt_alt {
                b'A' => 2,
                b'C' => 1,
                b'T' => 0,
                _ => 0,
            },
            b'T' => match nt_alt {
                b'A' => 3,
                b'C' => 4,
                b'G' => 5,
                _ => 0,
            },
            _ => 0,
        };

        const K: usize = N_NUCLEOTIDES;
        let k = match nt_5p {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0,
        };

        const L: usize = N_NUCLEOTIDES;
        let l = match nt_3p {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0,
        };

        (j * K + k) * L + l
    }

    pub fn channels() -> Vec<String> {
        let mut names = vec![String::new(); N_NONSTRANDED_MUTATION_CHANNELS as usize];
        for &nt_ref in [b'C', b'T'].iter() {
            for &nt_alt in NUCLEOTIDES.iter() {
                if nt_ref != nt_alt {
                    for &nt_5p in NUCLEOTIDES.iter() {
                        for &nt_3p in NUCLEOTIDES.iter() {
                            let idx = MutOpps::index(nt_ref, nt_alt, nt_5p, nt_3p);
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
    
}
