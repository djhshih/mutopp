extern crate bio;
extern crate multimap;
extern crate linked_hash_map;
#[macro_use] extern crate clap;
extern crate mutopp;

use std::io;
use std::str;
use std::fs;

use bio::io::fasta;
use bio::io::gff;

use bio::data_structures::interval_tree::{IntervalTree};

use multimap::MultiMap;
use linked_hash_map::LinkedHashMap;

use clap::{Arg, App, SubCommand};

use mutopp::gene::{Pos, Gene,Transcript,Region,Strand};
use mutopp::seq::{self,CodingSequence};
use mutopp::mutation::MutOpps;
use mutopp::utils;
use mutopp::io::snv;

type FastaIndexedReader = fasta::IndexedReader<fs::File>;
type GffReader = gff::Reader<fs::File>;
type SnvReader = snv::Reader<fs::File>;

const CDS_PADDING: usize = 3;

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
            .arg(Arg::with_name("snv")
                .help("input aggregated mutations tab-separated file (fields: chrom pos ref alt [sample])")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file")
                .required(true)))
        .subcommand(SubCommand::with_name("annotate")
            .about("annotate observed mutations")
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("gff")
                .help("input gene annotation .gff3 file")
                .required(true))
            .arg(Arg::with_name("snv")
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
        Some("enumerate") => {
            let m = matches.subcommand_matches("enumerate").unwrap();
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
        Some("count") => {
            let m = matches.subcommand_matches("count").unwrap();
            let ref fasta_fn = m.value_of("fasta").unwrap();
            let ref gff_fn = m.value_of("gff").unwrap();
            let ref snv_fn = m.value_of("snv").unwrap();
            let ref out_fn = m.value_of("output").unwrap();
            
            let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };
            
            let mut snv = match snv::Reader::from_file(snv_fn) {
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

            if let Err(why) = genes.write_counts(&mut snv, &mut ifasta, &mut out_file) {
                panic!("Could not write mutation type counts to file: {}", why);
            }
        },
        Some("annotate") => {
            let m = matches.subcommand_matches("annotate").unwrap();
            let ref fasta_fn = m.value_of("fasta").unwrap();
            let ref gff_fn = m.value_of("gff").unwrap();
            let ref snv_fn = m.value_of("snv").unwrap();
            let ref out_fn = m.value_of("output").unwrap();
            
            let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };
            
            let mut snv = match snv::Reader::from_file(snv_fn) {
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

            if let Err(why) = genes.write_annotations(&mut snv, &mut ifasta, &mut out_file) {
                panic!("Could not write mutation annotations to file: {}", why);
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

                for (_, gene) in genes.map.iter() {
                    for (tid, transcript) in gene.transcripts.iter() {
                        let padding = 3;
                        let cds = get_transcript_cds_from_fasta(&mut ifasta, transcript, &gene.chrom, gene.strand, padding);
                        println!(">{} gene_name={};padding={}", tid, gene.name, cds.padding);
                        for seq in cds.seqs {
                            seq::print_seq_padded(&seq, cds.padding);
                        }
                    }
                }
                
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
fn get_transcript_cds_from_fasta(reader: &mut FastaIndexedReader, transcript: &Transcript, chrom: &str, strand: Strand, padding: usize) -> CodingSequence {
    let n = transcript.coding_regions.len();
    let mut seqs: seq::DnaSeqs = vec![Vec::new(); n];
    for (i, region) in transcript.coding_regions.iter().enumerate() {
        if let Err(why) = reader.read(chrom, region.start - padding as u64, region.end + padding as u64, &mut seqs[i]) {
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

fn parse_phase(phase: &str) -> u8 {
    match phase {
        "0" => 0,
        "1" => 1,
        "2" => 2,
        _ => 0,
    }
}

#[derive(Debug)]
struct Genes {
    /// map maps gene id to gene struct
    map: LinkedHashMap<String, Gene>,
    /// index maps chromosome to interval tree
    /// interval tree stores interval and gene id
    index: LinkedHashMap<String, IntervalTree<Pos, String>>,
}

impl Genes {
    pub fn new() -> Genes {
        // upperbound estimates
        let ngenes: usize = 20000;
        let nchroms: usize = 25;
        Genes {
            map: LinkedHashMap::with_capacity(ngenes),
            index: LinkedHashMap::with_capacity(nchroms),
        }
    }
    
    /// Construct genes from GFF3 reader.
    ///
    /// GFF coordinates are 1-based, closed. Convert the coordinates to 0-based, half-open.
    /// Assume that genes are sorted by genomic start coordinate, the parent of each feature is
    /// always listed before the feature itself (e.g. gene before children transcripts, transcript
    /// before children CDS and other feature regions), and the CDS region are sorted in 5' to 3' order
    /// of the transcript.
    pub fn from_gff(reader: &mut GffReader) -> Genes {
        let mut genes = Genes::new();

        for r in reader.records() {
            let record = r.expect("Error reading GFF3 record");
            let attrs = &record.attributes();
            match record.feature_type() {
                "gene" => {
                    if attribute_filter(attrs, "gene_type", "protein_coding") {
                        genes.map.insert(
                            attrs.get("gene_id").unwrap().clone(),
                            Gene {
                                name: attrs.get("gene_name").unwrap().clone(),
                                chrom: record.seqname().to_owned(),
                                start: *record.start() - 1,
                                end: *record.end(),
                                strand: record.strand().unwrap(),
                                transcripts: LinkedHashMap::new(),
                            }
                        );
                    }
                },
                "transcript" => {
                    if attribute_filter(attrs, "transcript_type", "protein_coding") {
                        // assume that gene, if valid, has already been inserted
                        if let Some(gene) = genes.map.get_mut(attrs.get("gene_id").unwrap()) {
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
                    // CDS regions are in genomic coordinate but are listed in order of 5' to 3' of the gene
                    if let Some(gene) = genes.map.get_mut(attrs.get("gene_id").unwrap()) {
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

        // remove genes without valid transcripts
        let map_filtered = genes.map.into_iter()
            .filter(|&(_, ref v)| !v.transcripts.is_empty())
            .collect();
            
        let mut genes = Genes {
            map: map_filtered,
            index: genes.index,
        };
        genes.build_index();

        genes
    }
    
    #[inline]
    pub fn len(&self) -> usize {
        self.map.len()
    }
    
    /// Build interval tree index for genes.
    fn build_index(&mut self) {
        for (gid, gene) in self.map.iter() {
            let chrom_index = self.index.entry(gene.chrom.clone()).or_insert(IntervalTree::new());
            chrom_index.insert(gene.start .. gene.end, gid.clone());
        }
    }
    
    /// Find genes that overlap with genomic interval.
    pub fn find_overlap(&self, chrom: &str, start: Pos, end: Pos) -> Vec<String> {
        let mut overlap = Vec::new();

        if let Some(chrom_index) = self.index.get(chrom) {
            for r in chrom_index.find(start .. end) {
                let gid = r.data();
                overlap.push(gid.clone());
                //let gene = self.map.get(gid).expect("Gene id is not found.");
                //overlap.push(gene);
            }
        }
        
        overlap
    }
    
    /// Write the mutation opporunity array of each gene to file.
    /// mut is needed for ifasta due to indexing
    pub fn write_opps(&self, mut ifasta: &mut FastaIndexedReader, out: &mut fs::File) -> io::Result<()> {
        use std::io::Write;

        let sep = "\t";

        let mut header = vec![String::from("transcript_id"), String::from("gene_symbol")];
        header.extend(MutOpps::types());

        try!(writeln!(out, "{}", header.join(sep)));
        
        for (gid, gene) in self.map.iter() {
            for (tid, transcript) in gene.transcripts.iter() {
                let cds = get_transcript_cds_from_fasta(&mut ifasta, transcript, &gene.chrom, gene.strand, CDS_PADDING);
                let pseq = if cds.seqs.len() > 0 && cds.seqs[0].len() > CDS_PADDING + 6 {
                    str::from_utf8(&cds.seqs[0][0 .. CDS_PADDING]).unwrap().to_lowercase() +
                        str::from_utf8(&cds.seqs[0][CDS_PADDING .. CDS_PADDING + 6]).unwrap()
                } else {
                    String::new()
                };
                print!("{} {} {} {} ... ", gid, tid, gene.name, pseq);
                // sequence may not have 9 nucleotides!
                //println!("{} {} {} {}...", gid, gene.name, tid, str::from_utf8(&cds.seqs[0][0..9]).unwrap());
                if let Some(opp) = cds.count_opp() {
                    let mut line = vec![tid.clone(), gene.name.clone()].join(sep);
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
    
    pub fn write_counts(&self, mut snv: &mut SnvReader, mut ifasta: &mut FastaIndexedReader, out: &mut fs::File) -> io::Result<()> {
        use std::io::Write;
        
        let sep = "\t";
        
        // mutation type counts map -- key1: sample_id; key2: transcript_id
        let mut muts_map: LinkedHashMap<u32, LinkedHashMap<String, MutOpps>> = LinkedHashMap::new();
        
        // count mutation types
        for s in snv.records() {
            let record = s.ok().expect("Error reading SNV record");
            print!("{} ... ", record);
            
            let sample_muts_map = muts_map.entry(record.sample_id).or_insert(LinkedHashMap::new());
            let mut in_coding = false;
            let overlap = self.find_overlap(&record.chrom, record.pos, record.pos + 1);
            for gid in overlap.iter() {
                let gene = self.map.get(gid).unwrap();
                for (tid, transcript) in gene.transcripts.iter() {
                    let cds = get_transcript_cds_from_fasta(&mut ifasta, transcript, &gene.chrom, gene.strand, CDS_PADDING);
                    if let Some(effect) = cds.assign_mutation(gene, transcript, record.pos, record.nt_ref, record.nt_alt) {
                        let t_muts = sample_muts_map.entry(tid.clone()).or_insert(MutOpps::new());
                        t_muts[effect.type_index] += 1;
                        in_coding = true;
                    }
                }
            }
            
            if in_coding {
                println!("in");
            } else {
                println!("out");
            }
        }
        
        // write header
        let mut header = vec![String::from("sample_id"), String::from("transcript_id")];
        header.extend(MutOpps::types());
        try!(writeln!(out, "{}", header.join(sep)));
        
        // write opportunities to file
        for (sid, sample_muts_map) in muts_map.iter() {
            for (tid, t_muts) in sample_muts_map.iter() {
                let mut line = vec![sid.to_string(), tid.clone()].join(sep);
                for x in t_muts.iter() {
                    line.push_str(sep);
                    line.push_str(&x.to_string());
                }
                try!(writeln!(out, "{}", line));
            }
        }
        
        Ok(())
    }
    
    /// Write mutation annotations to file.
    pub fn write_annotations(&self, mut snv: &mut SnvReader, mut ifasta: &mut FastaIndexedReader, out: &mut fs::File) -> io::Result<()> {
        use std::io::Write;
        
        let sep = "\t";

        let header = vec![String::from("chrom"), String::from("pos"), String::from("ref"), 
        String::from("alt"), String::from("sample_id"), String::from("transcript_id"), String::from("gene_symbol"), 
        String::from("cdna_change"), String::from("protein_change"), String::from("impact"), String::from("type_index")];
        
        try!(writeln!(out, "{}", header.join(sep)));
        
        for s in snv.records() {
            let record = s.ok().expect("Error reading SNV record");
            println!("{}", record);

            let overlap = self.find_overlap(&record.chrom, record.pos, record.pos + 1);
            for gid in overlap.iter() {
                let gene = self.map.get(gid).unwrap();
                println!("{}: {{", gene.name);
                for (tid, transcript) in gene.transcripts.iter() {
                    print!("    {}: ", tid);
                    let cds = get_transcript_cds_from_fasta(&mut ifasta, transcript, &gene.chrom, gene.strand, CDS_PADDING);
                    if let Some(effect) = cds.assign_mutation(gene, transcript, record.pos, record.nt_ref, record.nt_alt) {
                        let line = vec![record.chrom.clone(), (record.pos + 1).to_string(),
                            str::from_utf8(&[record.nt_ref]).unwrap().to_owned(),
                            str::from_utf8(&[record.nt_alt]).unwrap().to_owned(),
                            record.sample_id.to_string(), tid.clone(), gene.name.clone(),
                            effect.cdna_change(), effect.protein_change(), effect.impact.to_string(), effect.type_index.to_string()].join(sep);
                        try!(writeln!(out, "{}", line));
                        println!("{},", effect.impact.to_string());
                    } else {
                        println!("none,");
                    }
                }
                println!("}}");
            }
        }
        
        Ok(())
    }
    
}