extern crate bio;
extern crate multimap;
#[macro_use] extern crate clap;
extern crate mutopp;

use std::io;
use std::str;
use std::fs;

use bio::io::fasta;
use bio::io::gff;

use multimap::MultiMap;

use clap::{Arg, App, SubCommand};

use mutopp::gene::{Gene,Transcript,Region,HashMap,Strand};
use mutopp::seq::{self,CodingSequence};
use mutopp::mutation::MutOpps;
use mutopp::utils;

type FastaIndexedReader = fasta::IndexedReader<fs::File>;
type GffReader = gff::Reader<fs::File>;

fn main() {
    let matches = App::new("mutopp")
        .version(crate_version!())
        .author("David J. H. Shih <djh.shih@gmail.com>")
        .about("Mutation opportunity")
        .subcommand(SubCommand::with_name("enumerate")
            .about("enumerate all mutation opportunities")
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("gff")
                .help("input gene annotation .gff3 file")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file")
                .required(true)))
        .subcommand(SubCommand::with_name("count")
            .about("count mutation types in observed mutations")
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("gff")
                .help("input gene annotation .gff3 file")
                .required(true))
            .arg(Arg::with_name("mut")
                .help("input aggregated mutations tab-separated file (fields: chrom pos ref alt [sample])")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file")
                .required(true)))
        .subcommand(SubCommand::with_name("extract")
            .about("extract sequence from .fasta file")
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("coord")
                .short("c")
                .long("coord")
                .help("target coordinate to extract")
                .takes_value(true))
            .arg(Arg::with_name("gff")
                .short("g")
                .long("gff")
                .help("genes to extract")
                .conflicts_with("coord")
                .takes_value(true)))
        .get_matches();

    match matches.subcommand_name() {
        Some("count") => {
            let m = matches.subcommand_matches("count").unwrap();
            let ref fasta_fn = m.value_of("fasta").unwrap();
            let ref gff_fn = m.value_of("gff").unwrap();
            let ref out_fn = m.value_of("output").unwrap();
            
            let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };

            let mut gff = match gff::Reader::from_file(gff_fn, gff::GffType::GFF3) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };

            let genes = Genes::from_gff(&mut gff);
            println!("number of valid protein-coding genes: {}", genes.len());
            println!("{:?}", genes);
            
            // Open a file in write-only mode, returns `io::Result<fs::File>`
            let mut out_file = match fs::File::create(&out_fn) {
                Err(why) => panic!("Could not create {}: {:?}",
                                out_fn,
                                why),
                Ok(file) => file,
            };

            if let Err(why) = genes.write_opps(&mut ifasta, &mut out_file) {
                panic!("Could not write mutation opportunities to file: {}", why);
            }
        },
        Some("extract") => {
            let m = matches.subcommand_matches("extract").unwrap();
            let ref fasta_fn = m.value_of("fasta").unwrap();
            
            let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };
            
            if let Some(ref gff_fn) = m.value_of("gff") {
                let mut gff = match gff::Reader::from_file(gff_fn, gff::GffType::GFF3) {
                    Err(why) => panic!("{:?}", why),
                    Ok(reader) => reader,
                };

                let genes = Genes::from_gff(&mut gff);
                println!("number of valid protein-coding genes: {}", genes.len());
                println!("{:?}", genes);
                
                // TODO extract and output sequence for genes
            }

            if let Some(coord) = m.value_of("coord") {
                // read only part of the fasta
                if let Some((seqname, start, end)) = utils::parse_coordinate(&coord) {
                    //println!("{} {} {}", seqname, start, end);
                    let mut seq: Vec<u8> = Vec::new();
                    match ifasta.read(&seqname, start, end, &mut seq) {
                        Err(why) => panic!("{:?}", why),
                        Ok(()) => seq::print_seq(&seq),
                    }
                } else {
                    panic!("Invalid coordinate: {}", coord);
                }
            }
            
        },
        None => println!("Type `mutopp help` for help."),
        _ => panic!("Invalid subcommand"),
    }
}

/// Assume that CDS are in order for features on the plus and the minus strand.
/// For features on the minus strand specifically, the positions of the CDS is in descending order
fn get_transcript_cds_from_fasta(reader: &mut FastaIndexedReader, transcript: &Transcript, chromosome: &str, strand: Strand, padding: usize) -> CodingSequence {
    let n = transcript.coding_regions.len();
    let mut seqs: seq::DnaSeqs = vec![Vec::new(); n];
    for (i, region) in transcript.coding_regions.iter().enumerate() {
        if let Err(why) = reader.read(chromosome, region.start - padding as u64, region.end + padding as u64, &mut seqs[i]) {
            panic!("{:?}", why);
        }
        if strand == Strand::Reverse {
            seq::reverse_complement(&mut seqs[i]);
        }
    }
    
    CodingSequence { seqs: seqs, padding: padding }
}

fn attribute_filter(attrs: &MultiMap<String, String>, key: &str, target: &str) -> bool {
    match attrs.get(key) {
        Some(value) => (value == target),
        None => false,
    }
}

fn parse_phase(phase: &str) -> i8 {
    match phase {
        "0" => 0,
        "1" => 1,
        "2" => 2,
        _ => -1,
    }
}

#[derive(Debug)]
struct Genes(HashMap<String, Gene>);

impl Genes {
    fn new() -> Genes {
        Genes(HashMap::new())
    }
    
    /// Construct genes from GFF reader.
    ///
    /// GFF coordinates are 1-based, closed. Convert the coordinates to 0-based, half-open.
    fn from_gff(reader: &mut GffReader) -> Genes {
        let mut genes = Genes::new();

        for r in reader.records() {
            let record = r.expect("Error reading GFF3 record");
            let attrs = &record.attributes();
            match record.feature_type() {
                "gene" => {
                    if attribute_filter(attrs, "gene_type", "protein_coding") {
                        genes.0.insert(
                            attrs.get("gene_id").unwrap().clone(),
                            Gene {
                                name: attrs.get("gene_name").unwrap().clone(),
                                chromosome: record.seqname().to_owned(),
                                start: *record.start() - 1,
                                end: *record.end(),
                                strand: record.strand().unwrap(),
                                transcripts: HashMap::new(),
                            }
                        );
                    }
                },
                "transcript" => {
                    if attribute_filter(attrs, "transcript_type", "protein_coding") {
                        // assume that gene, if valid, has already been inserted
                        if let Some(gene) = genes.0.get_mut(attrs.get("gene_id").unwrap()) {
                            gene.transcripts.insert(
                                attrs.get("transcript_id").unwrap().clone(),
                                Transcript {
                                    start: *record.start() - 1,
                                    end: *record.end(),
                                    coding_regions: Vec::new(),
                                }
                            );
                        }
                    }
                },
                "CDS" => {
                    // assume that gene and transcript, if valid, have already been inserted
                    // assume that frame of the first CDS is always 0
                    if let Some(gene) = genes.0.get_mut(attrs.get("gene_id").unwrap()) {
                        if let Some(transcript) = gene.transcripts.get_mut(attrs.get("transcript_id").unwrap()) {
                            transcript.coding_regions.push(
                                Region { phase: parse_phase(record.frame()), start: *record.start() - 1, end: *record.end() }
                            );
                        }
                    }
                },
                _ => {} // ignore all other fields
            }
        }

        // filter genes without valid transcripts
        Genes(
            genes.0.into_iter()
                .filter(|&(_, ref v)| !v.transcripts.is_empty())
                .collect()
        )
    }

    #[inline]
    fn len(&self) -> usize {
        self.0.len()
    }
    
    /// Write the mutation opporunity array of each gene to file.
    /// mut is needed for ifasta due to indexing
    fn write_opps(&self, mut ifasta: &mut FastaIndexedReader, out: &mut fs::File) -> io::Result<()> {
        use std::io::Write;

        let sep = "\t";

        let mut header = vec![String::from("gene"), String::from("symbol"), String::from("transcript")];
        header.extend(MutOpps::types());

        try!(writeln!(out, "{}", header.join(sep)));

        for (gid, gene) in self.0.iter() {
            for (tid, transcript) in gene.transcripts.iter() {
                let padding = 3;
                let cds = get_transcript_cds_from_fasta(&mut ifasta, transcript, &gene.chromosome, gene.strand, padding);
                let pseq = if cds.seqs.len() > 0 && cds.seqs[0].len() > padding + 6 {
                    str::from_utf8(&cds.seqs[0][0 .. padding]).unwrap().to_lowercase() +
                        str::from_utf8(&cds.seqs[0][padding .. padding + 6]).unwrap()
                } else {
                    String::new()
                };
                print!("{} {} {} {}... ", gid, gene.name, tid, pseq);
                // sequence may not have 9 nucleotides!
                //println!("{} {} {} {}...", gid, gene.name, tid, str::from_utf8(&cds.seqs[0][0..9]).unwrap());
                if let Some(opp) = cds.count_opp() {
                    let mut line = vec![gid.clone(), gene.name.clone(), tid.clone()].join(sep);
                    for o in opp.iter() {
                        line.push_str(sep);
                        line.push_str(&o.to_string());
                    }
                    try!(writeln!(out, "{}", line));
                    println!("succeeded");
                } else {
                    println!("failed");
                }
            }
        }
        
        Ok(())
    }
}