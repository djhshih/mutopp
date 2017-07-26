use std::io;
use std::fs;
use std::path::Path;
use std::convert::AsRef;
use std::result;
use std::num;

use csv;
use crc::crc32;

/// A SNV reader.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}

impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::new)
    }
}

impl<R: io::Read> Reader<R> {
    /// Read from a given reader.
    pub fn new(reader: R) -> Self {
        Reader {
            inner: csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .comment(Some(b'#'))
                .has_headers(true)
                .from_reader(reader),
        }
    }
    
    /// Iterate over records.
    pub fn records(&mut self) -> Records<R> {
        Records { inner: self.inner.records() }
    }
}

pub struct Records<'r, R: 'r + io::Read> {
    inner: csv::StringRecordsIter<'r, R>,
}

#[derive(Debug)]
pub enum Error {
    Csv(csv::Error),
    MissingField(String),
    ParseInt(num::ParseIntError),
}

pub type Result<T> = result::Result<T, Error>;

impl<'r, R: io::Read> Iterator for Records<'r, R> {
    type Item = Result<Record>;
    
    /// Get next record.
    /// Stop reading as soon as a problematic record is encountered.
    fn next(&mut self) -> Option<Result<Record>> {
        // StringRecordsIter next() produces doubly wrapped values:
        // Option<csv::Result<StringRecord>>
        // therefore, we need to map twice
        self.inner.next()
            .map(|res| {
                match res {
                    Err(err) => Err(Error::Csv(err)),
                    Ok(record) => Ok( Record {
                        chrom: try!(record.get(0)
                            .ok_or(Error::MissingField("chrom".to_owned()))
                            .map(|x| String::from(x))),
                        pos: try!(record.get(1)
                            .ok_or(Error::MissingField("pos".to_owned()))
                            .and_then(|x| x.parse::<u64>().map_err(Error::ParseInt))),
                        ref_nt: try!(record.get(2)
                            .ok_or(Error::MissingField("ref".to_owned()))
                            .map(|x| x.chars().nth(0).unwrap() as u8)),
                        alt_nt: try!(record.get(3)
                            .ok_or(Error::MissingField("alt".to_owned()))
                            .map(|x| x.chars().nth(0).unwrap() as u8)),
                        sample: record.get(4)
                            .and_then(|x| match x.parse::<u32>() {
                                // apply CRC32 hash to string
                                Err(_) => Some(crc32::checksum_ieee(x.as_bytes())),
                                Ok(i) => Some(i),
                            }).unwrap_or(0),
                    }),
                }
            })
    }
}

/// A SNV record.
pub struct Record {
    /// Chromosome or contig name
    chrom: String,
    /// Genomic position
    pos: u64,
    /// Reference nucleotide on the positive/reference strand
    ref_nt: Nucleotide,
    /// Observed alternate nucleotide on the positive/reference strand
    alt_nt: Nucleotide,
    /// Sample ID
    sample: u32,
}

#[cfg(test)]
mod tests {
    use super::*;

    const SNV_FILE: &'static [u8] = b"chrom\tpos\tref\talt\tsample
2\t1525\tA\tC\t1
3\t1602\tG\tC\tT002
";

    #[test]
    fn test_reader() {
        let chroms = ["2", "3"];
        let poss = [1525, 1602];
        let ref_nts = [b'A', b'G'];
        let alt_nts = [b'C', b'C'];
        let samples = [1, 0x56d1c50a];
        
        let mut reader = Reader::new(SNV_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.ok().expect("Error reading record");
            assert_eq!(record.chrom, chroms[i]);
            assert_eq!(record.pos, poss[i]);
            assert_eq!(record.ref_nt, ref_nts[i]);
            assert_eq!(record.alt_nt, alt_nts[i]);
            assert_eq!(record.sample, samples[i]);
        }
    }
}