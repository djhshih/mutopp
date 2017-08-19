use gene::Region;
use mutation::genomic;
use seq::{Nucleotide,DnaSeq};
use constants::*;

/// Contiguous genomic sequence.
pub struct Sequence {
    pub inner: DnaSeq,
}

impl Sequence {
    /// Assign channel to observed mutation.
    // TODO this method should be moved elsewhere: we do not need to entire sequence!
    pub fn assign_channels(&self, region: &Region, pos: u64, nt_ref: Nucleotide, nt_alt:Nucleotide) -> Option<usize> {
        let seq = &self.inner;
        let n = seq.len();

        if n < 3 {
            return None;
        }

        // check that region size matches sequence length
        if n != (region.end - region.start) as usize {
            return None;
        }

        // mutation cannot occur at first local index (0), because no 5' context is available
        // mutation cannot occur at last local index (n - 1), because no 3' context is available
        if pos > region.start && pos < region.end - 1 {
            // previous checks ensure that local_idx is value
            let local_idx = (pos - region.start) as usize;

            // check that mutation reference allele match refernce sequence
            if seq[local_idx] != nt_ref {
                return None;
            }

            // previous local index bound check ensures that nt_5p and nt_3p are well defined
            let nt_5p = seq[local_idx - 1];
            let nt_3p = seq[local_idx + 1];
            
            Some(genomic::MutOpps::index(nt_ref, nt_alt, nt_5p, nt_3p))
        } else {
            None
        }
    }

    /// Count mutation opportunities.
    pub fn count_opp(&self) -> Option<genomic::MutOpps> {
        let seq = &self.inner;
        let n = seq.len();

        // context of 3 nucleotide length
        if n < 3 {
            return None;
        }

        let mut x = genomic::MutOpps::new();
        self.accumulate_opp(&mut x);
        
        Some(x)
    }

    /// Accumulate mutation opportunities.
    pub fn accumulate_opp(&self, x: &mut genomic::MutOpps) {
        let seq = &self.inner;
        let n = seq.len();

        // context of 3 nucleotide length
        if n < 3 {
            return;
        }

        let mut nt_5p = seq[0];
        let mut nt_ref = seq[1];
        // nt_3p will be returned by the iterator below

        // iterate through each nucleotide, starting from the third
        let mut it = seq[2..].iter();
        loop {
            match it.next() {
                Some(&nt_3p) => {
                    // TODO avoid checking each nucleotide three times
                    if nt_5p != b'N' && nt_ref != b'N' && nt_3p != b'N' {
                        for &nt_alt in NUCLEOTIDES.iter() {
                            if nt_alt != nt_ref {
                                let idx = genomic::MutOpps::index(nt_ref, nt_alt, nt_5p, nt_3p);
                                x[idx] += 1;
                            }
                        }
                    }
                    // shift window for next iteration
                    nt_5p = nt_ref;
                    nt_ref = nt_3p;
                },
                None => break,
            }
        }
    }
}
