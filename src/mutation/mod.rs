use std::fmt;
use std::ops;
use std::ascii::AsciiExt;
use std::slice;

use constants::*;

#[derive(Copy,Clone)]
pub enum MutImpact {
    Synonymous,
    Missense,
    StartOrStop,
    SpliceSite,
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
        };

        write!(f, "{}", c)
    }
}

/*
custom_derive! {
    #[derive(NewtypeDebug)]
    pub struct MutOpps([u32; n_mutation_types]);
}
*/

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

