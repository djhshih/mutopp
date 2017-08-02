use std::io;
use std::fs;
use std::path::Path;
use std::convert::AsRef;
use std::result;
use std::num;
use std::fmt;

use csv;
use crc::crc32;

use seq::{Nucleotide,Residue};

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
                            .and_then(|x| x.parse::<u64>().map_err(Error::ParseInt))
                            // convert from 1-based to 0-based coordinate
                            .map(|x| x - 1)),
                        nt_ref: try!(record.get(2)
                            .ok_or(Error::MissingField("ref".to_owned()))
                            .map(|x| x.chars().nth(0).unwrap() as u8)),
                        nt_alt: try!(record.get(3)
                            .ok_or(Error::MissingField("alt".to_owned()))
                            .map(|x| x.chars().nth(0).unwrap() as u8)),
                        sample_id: record.get(4)
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
    pub chrom: String,
    /// Genomic position
    pub pos: u64,
    /// Reference nucleotide on the positive/reference strand
    pub nt_ref: Nucleotide,
    /// Observed alternate nucleotide on the positive/reference strand
    pub nt_alt: Nucleotide,
    /// Sample ID
    pub sample_id: u32,
}


impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}@{}:g.{}{}>{}", self.sample_id, self.chrom, self.pos + 1, self.nt_ref as char, self.nt_alt as char)
    }
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
        let nts_ref = [b'A', b'G'];
        let nts_alt = [b'C', b'C'];
        let samples = [1, 0x56d1c50a];
        
        let mut reader = Reader::new(SNV_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.ok().expect("Error reading record");
            assert_eq!(record.chrom, chroms[i]);
            assert_eq!(record.pos + 1, poss[i]);
            assert_eq!(record.nt_ref, nts_ref[i]);
            assert_eq!(record.nt_alt, nts_alt[i]);
            assert_eq!(record.sample_id, samples[i]);
        }
    }
}