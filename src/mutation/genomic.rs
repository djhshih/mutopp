use std::fmt;
use std::ops;
use std::ascii::AsciiExt;
use std::slice;

use std::path::Path;

use std::fs;
use std::io;
use std::io::BufRead;

use constants::*;
use seq;

/// Mutation signature spectrum
pub struct MutSpec([f64; N_NONSTRANDED_MUTATION_CHANNELS]);

impl MutSpec {
    #[inline]
    pub fn new() -> MutSpec {
        MutSpec([0.0; N_NONSTRANDED_MUTATION_CHANNELS])
    }

    pub fn from_file<P: AsRef<Path>>(path: P) -> MutSpec {
        let reader = match fs::File::open(path) {
            Err(why) => panic!("{:?}", why),
            Ok(f) => io::BufReader::new(f),
        };

        let mut mutspec = MutSpec::new();
        let mut i = 0;
        for line in reader.lines() {
            let l = line.unwrap();
            mutspec[i] = l.parse().unwrap(); 
            if i == MutSpec::len() {
                break;
            }
            i += 1;
        }

        mutspec
    }

    #[inline]
    pub fn iter(&self) -> slice::Iter<f64> {
        self.0.iter()
    }

    #[inline]
    pub fn index(nt_ref: u8, nt_alt: u8, nt_5p: u8, nt_3p: u8) -> usize {
        channel_index(nt_ref, nt_alt, nt_5p, nt_3p)
    }

    #[inline]
    pub fn channels() -> Vec<String> {
        channel_vector()
    }

    #[inline]
    pub fn len() -> usize {
        N_NONSTRANDED_MUTATION_CHANNELS
    }

}

impl fmt::Debug for MutSpec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(&self.0 as &[f64], f)
    }
}

impl ops::Index<usize> for MutSpec {
    type Output = f64;

    #[inline]
    fn index(&self, index: usize) -> &f64 {
        &self.0[index]
    }
}

impl ops::IndexMut<usize> for MutSpec {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        &mut self.0[index]
    }
}


/// Genomic mutation opportunites
pub struct MutOpps([u32; N_NONSTRANDED_MUTATION_CHANNELS]);

impl MutOpps {
    #[inline]
    pub fn new() -> MutOpps {
        MutOpps([0; N_NONSTRANDED_MUTATION_CHANNELS])
    }

    #[inline]    
    pub fn iter(&self) -> slice::Iter<u32> {
        self.0.iter()
    }
    
    #[inline]
    pub fn index(nt_ref: u8, nt_alt: u8, nt_5p: u8, nt_3p: u8) -> usize {
        channel_index(nt_ref, nt_alt, nt_5p, nt_3p)
    }

    #[inline]
    pub fn channels() -> Vec<String> {
        channel_vector()
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


fn channel_vector() -> Vec<String> {
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

// mutation channels are indexed by factors: impact class, nonstranded substitution, 5' context, and 3' context
// each factor serve as a subscript into the mutation channel array
// consider nucleotides in the order of A, C, G, T
// nonstranded substitutions (6): C>A, C>G, C>T, T>A, T>C, T>G
//                                G>T, G>C, G>A, A>T, A>G, A>C
// 5' context (4): A, C, G, T
// 3' context (4): A, C, G, T
// number of channels = 6 * 4 * 4 = 96
// later factors vary faster than early factors
fn channel_index(nt_ref: u8, nt_alt: u8, nt_5p: u8, nt_3p: u8) -> usize {
    assert!(nt_ref != nt_alt);

    // non-stranded mutation context
    let nt_ref_ns;
    let nt_alt_ns;
    let nt_5p_ns;
    let nt_3p_ns;

    match nt_ref {
        // mutation involves C or T: proceed as is
        b'C' | b'T' => {
            nt_ref_ns = nt_ref;
            nt_alt_ns = nt_alt;
            nt_5p_ns = nt_5p;
            nt_3p_ns = nt_3p;
        },
        // mutation does not involve C or T: need to reverse complement the mutation context
        _ => {
            nt_ref_ns = seq::complement(nt_ref);
            nt_alt_ns = seq::complement(nt_alt);
            // flip 5' and 3' as well as complement
            nt_5p_ns = seq::complement(nt_3p);
            nt_3p_ns = seq::complement(nt_5p);
        },
    }

    let j = match nt_ref_ns {
        b'A' => match nt_alt_ns {
            b'C' => 5,
            b'G' => 4,
            b'T' => 3,
            _ => 0,
        },
        b'C' => match nt_alt_ns {
            b'A' => 0,
            b'G' => 1,
            b'T' => 2,
            _ => 0,
        },
        b'G' => match nt_alt_ns {
            b'A' => 2,
            b'C' => 1,
            b'T' => 0,
            _ => 0,
        },
        b'T' => match nt_alt_ns {
            b'A' => 3,
            b'C' => 4,
            b'G' => 5,
            _ => 0,
        },
        _ => 0,
    };

    const K: usize = N_NUCLEOTIDES;
    let k = match nt_5p_ns {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    };

    const L: usize = N_NUCLEOTIDES;
    let l = match nt_3p_ns {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    };

    (j * K + k) * L + l
}
