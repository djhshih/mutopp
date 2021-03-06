extern crate bio;
extern crate bio_types;
extern crate multimap;
extern crate linked_hash_map;
#[macro_use] extern crate clap;
extern crate mutopp;

extern crate rand;
use rand::distributions::IndependentSample;

use std::io;
use std::str;
use std::fs;

use bio::io::fasta;
use bio::io::gff;
use bio::io::bed;

use bio::data_structures::interval_tree::{IntervalTree};

use multimap::MultiMap;
use linked_hash_map::LinkedHashMap;

use clap::{Arg, App, SubCommand};

use mutopp::gene::{Pos, Gene,Transcript,Region,Strand};
use mutopp::seq::{self,Nucleotide};
use mutopp::mutation;
use mutopp::utils;
use mutopp::io::snv;
use mutopp::sample;
use mutopp::stats;

type FastaIndexedReader = fasta::IndexedReader<fs::File>;
type GffReader = gff::Reader<fs::File>;
type SnvReader = snv::Reader<fs::File>;
type BedReader = bed::Reader<fs::File>;

const CDS_PADDING: usize = 3;

fn main() {
    let matches = App::new("mutopp")
        .version(crate_version!())
        .author("David J. H. Shih <djh.shih@gmail.com>")
        .about("Mutation opportunity")
        .subcommand(SubCommand::with_name("enum-genes")
            .about("enumerate all mutation opportunities in genes")
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("gff")
                .help("input gene annotation .gff3 file")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file")
                .required(true)))
        .subcommand(SubCommand::with_name("count-gene-muts")
            .about("count mutation types of observed mutations")
            .arg(Arg::with_name("aggregate")
                .short("g")
                .long("aggregate")
                .help("aggregate mutations from different samples together"))
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("gff")
                .help("input gene annotation .gff3 file")
                .required(true))
            .arg(Arg::with_name("snv")
                .help("input mutations as tab-separated file (fields: chrom pos ref alt [sample])")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file")
                .required(true)))
        .subcommand(SubCommand::with_name("annot-gene-muts")
            .about("annotate observed mutations")
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("gff")
                .help("input gene annotation .gff3 file")
                .required(true))
            .arg(Arg::with_name("snv")
                .help("input mutations as tab-separated file (fields: chrom pos ref alt [sample])")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file")
                .required(true)))
        .subcommand(SubCommand::with_name("enum-regions")
            .about("enumerate all mutation opportunities")
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("bed")
                .help("input regions .bed file")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file")
                .required(true)))
        .subcommand(SubCommand::with_name("count-region-muts")
            .about("count mutation channels of observed mutations")
            .arg(Arg::with_name("aggregate")
                .short("g")
                .long("aggregate")
                .help("aggregate mutations from different samples together"))
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("bed")
                .help("input regions .bed file")
                .required(true))
            .arg(Arg::with_name("snv")
                .help("input mutations as tab-separated file (fields: chrom pos ref alt [sample])")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file")
                .required(true)))
        .subcommand(SubCommand::with_name("extract-fasta")
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
        .subcommand(SubCommand::with_name("gff-to-bed")
            .about("convert transcript information from GFF3 to BED format")
            .arg(Arg::with_name("gff")
                .help("input gene annotation .gff3 file")
                .required(true))
            .arg(Arg::with_name("bed")
                .help("output gene annotation .bed file")
                .required(true)))
        .subcommand(SubCommand::with_name("sample-snv")
            .about("sample coding SNV according to a mutation spectrum")
            .arg(Arg::with_name("fasta")
                .help("input indexed .fasta file (index .fasta.fai file must exist)")
                .required(true))
            .arg(Arg::with_name("gff")
                .help("input gene annotation .gff3 file")
                .required(true))
            .arg(Arg::with_name("spectrum")
                .help("input mutation spectrum .vtr file")
                .required(true))
            .arg(Arg::with_name("output")
                .help("output file for the sample")
                .required(true))
            .arg(Arg::with_name("sample-size")
                .help("number of samples to draw")
                .short("n")
                .long("size")
                .takes_value(true))
            .arg(Arg::with_name("output-imp")
                .help("output file for importance sample")
                .short("i")
                .long("output-imp")
                .takes_value(true)))
        .get_matches();

    match matches.subcommand_name() {
        Some("enum-genes") => {
            let m = matches.subcommand_matches("enum-genes").unwrap();
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
        Some("count-gene-muts") => {
            let m = matches.subcommand_matches("count-gene-muts").unwrap();
            let ref fasta_fn = m.value_of("fasta").unwrap();
            let ref gff_fn = m.value_of("gff").unwrap();
            let ref snv_fn = m.value_of("snv").unwrap();
            let ref out_fn = m.value_of("output").unwrap();
            let aggregate = m.is_present("aggregate");
            
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
            
            // Open a file in write-only mode, returns `io::Result<fs::File>`
            let mut out_file = match fs::File::create(&out_fn) {
                Err(why) => panic!("Could not create {}: {:?}",
                                out_fn,
                                why),
                Ok(file) => file,
            };

            if let Err(why) = genes.write_counts(&mut snv, &mut ifasta, &mut out_file, aggregate) {
                panic!("Could not write mutation type counts to file: {}", why);
            }
        },
        Some("annot-gene-muts") => {
            let m = matches.subcommand_matches("annot-gene-muts").unwrap();
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
        Some("enum-regions") => {
            let m = matches.subcommand_matches("enum-regions").unwrap();
            let ref fasta_fn = m.value_of("fasta").unwrap();
            let ref bed_fn = m.value_of("bed").unwrap();
            let ref out_fn = m.value_of("output").unwrap();
            
            let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };

            let mut bed = match bed::Reader::from_file(bed_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };

            let regions = Regions::from_bed(&mut bed);

            // Open a file in write-only mode, returns `io::Result<fs::File>`
            let mut out_file = match fs::File::create(&out_fn) {
                Err(why) => panic!("Could not create {}: {:?}",
                                out_fn,
                                why),
                Ok(file) => file,
            };

            if let Err(why) = regions.write_opps(&mut ifasta, &mut out_file) {
                panic!("Could not write mutation opportunities to file: {}", why);
            }
        },
        Some("count-region-muts") => {
            let m = matches.subcommand_matches("count-region-muts").unwrap();
            let ref fasta_fn = m.value_of("fasta").unwrap();
            let ref bed_fn = m.value_of("bed").unwrap();
            let ref snv_fn = m.value_of("snv").unwrap();
            let ref out_fn = m.value_of("output").unwrap();
            let aggregate = m.is_present("aggregate");
            
            let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };
            
            let mut snv = match snv::Reader::from_file(snv_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };

            let mut bed = match bed::Reader::from_file(bed_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };
            
            let regions = Regions::from_bed(&mut bed);
            
            // Open a file in write-only mode, returns `io::Result<fs::File>`
            let mut out_file = match fs::File::create(&out_fn) {
                Err(why) => panic!("Could not create {}: {:?}",
                                out_fn,
                                why),
                Ok(file) => file,
            };

            if let Err(why) = regions.write_counts(&mut snv, &mut ifasta, &mut out_file, aggregate) {
                panic!("Could not write mutation channel counts to file: {}", why);
            }
        },
        Some("extract-fasta") => {
            let m = matches.subcommand_matches("extract-fasta").unwrap();
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
                println!("; number of valid protein-coding genes: {}", genes.len());
                println!("; {:?}", genes.map);

                for (_, gene) in genes.map.iter() {
                    for (tid, transcript) in gene.transcripts.iter() {
                        let padding = 3;
                        let cds = get_transcript_cds_from_fasta(&mut ifasta, transcript, &gene.chrom, gene.strand, padding);
                        let lens: Vec<usize> = cds.seqs.iter().map(|x| x.len() - 2 * padding).collect();
                        let len_total: usize = lens.iter().sum();
                        let len_str: Vec<String> = lens.iter().map(|x| x.to_string()).collect();
                        println!(">{} gene_name={};lengths={};total={};padding={};", tid, gene.name, len_str.join(","), len_total, cds.padding);
                        for seq in cds.seqs {
                            seq::print_seq_padded(&seq, cds.padding);
                        }
                    }
                }
            }

            if let Some(coord) = m.value_of("coord") {
                // read only part of the fasta
                if let Some((seqname, start, end)) = utils::parse_coordinate(&coord) {
                    //println!("{} {} {}", seqname, start, end);
                    match ifasta.fetch(&seqname, start, end) {
                        Err(why) => panic!("{:?}", why),
                        Ok(()) => {
                            let mut seq: Vec<u8> = Vec::new();
                            ifasta.read(&mut seq);
                            seq::print_seq(&seq)
                        }
                    }
                } else {
                    panic!("Invalid coordinate: {}", coord);
                }
            }
            
        },
        Some("gff-to-bed") => {
            let m = matches.subcommand_matches("gff-to-bed").unwrap();
            let ref gff_fn = m.value_of("gff").unwrap();
            let ref out_fn = m.value_of("bed").unwrap();

            let mut gff = match gff::Reader::from_file(gff_fn, gff::GffType::GFF3) {
                    Err(why) => panic!("{:?}", why),
                    Ok(reader) => reader,
            };

            let genes = Genes::from_gff(&mut gff);
            println!("number of valid protein-coding genes: {}", genes.len());

            // Open a file in write-only mode, returns `io::Result<fs::File>`
            let mut out_file = match fs::File::create(&out_fn) {
                    Err(why) => panic!("Could not create {}: {:?}",
                                                    out_fn,
                                                    why),
                    Ok(file) => file,
            };

            if let Err(why) = genes.write_bed(&mut out_file) {
                panic!("Could not write gene annotations to BED file: {}", why);
            }
        },
        Some("sample-snv") => {
            let m = matches.subcommand_matches("sample-snv").unwrap();
            let ref fasta_fn = m.value_of("fasta").unwrap();
            let ref gff_fn = m.value_of("gff").unwrap();
            let ref spec_fn = m.value_of("spectrum").unwrap();
            let ref out_fn = m.value_of("output").unwrap();
            let nsamples = value_t!(m, "sample-size", usize).unwrap_or(10);
            
            let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };

            let mut gff = match gff::Reader::from_file(gff_fn, gff::GffType::GFF3) {
                Err(why) => panic!("{:?}", why),
                Ok(reader) => reader,
            };

            let mut out = match fs::File::create(&out_fn) {
                Err(why) => panic!("Could not create {}: {:?}", out_fn, why),
                Ok(file) => file,
            };

            let mut out_imp: fs::File;
            let out_imp_opt = match m.value_of("output-imp") {
                Some(out_imp_fn) => {
                    match fs::File::create(&out_imp_fn) {
                        Err(why) => panic!("Could not create {}: {:?}", out_imp_fn, why),
                        Ok(file) => {
                            out_imp = file;
                            Some(&mut out_imp)
                        },
                    }
                },
                None => None,
            };


            let genes = Genes::from_gff(&mut gff);
            println!("number of valid protein-coding genes: {}", genes.len());

            let mutspec = mutation::genomic::MutSpec::from_file(spec_fn);
            //println!("{:?}", mutspec);

            if let Err(why) = write_snv_sample(&mut ifasta, &genes, &mutspec, &mut out, nsamples, out_imp_opt) {
                panic!("Could not write sampled SNV to out file: {}", why);
            }
        },
        None => println!("Type `mutopp help` for help."),
        _ => panic!("Invalid subcommand"),
    }
}

/// Assume that CDS are in order for features on the plus and the minus strand.
/// For features on the minus strand specifically, the positions of the CDS is in descending order
fn get_transcript_cds_from_fasta(reader: &mut FastaIndexedReader, transcript: &Transcript, chrom: &str, strand: Strand, padding: usize) -> seq::coding::Sequence {
    let n = transcript.coding_regions.len();
    let mut seqs: Vec<seq::DnaSeq> = vec![Vec::new(); n];
    for (i, region) in transcript.coding_regions.iter().enumerate() {
        if let Err(why) = reader.fetch(chrom, region.start - padding as u64, region.end + padding as u64) {
            panic!("{:?}", why);
        }
        reader.read(&mut seqs[i]);
        if strand == Strand::Reverse {
            seq::reverse_complement(&mut seqs[i]);
        }
    }
    
    seq::coding::Sequence { seqs: seqs, padding: padding }
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
    /// vector of gene names
    names: Vec<String>,
}

impl Genes {
    pub fn new() -> Genes {
        // upperbound estimates
        let ngenes: usize = 20000;
        let nchroms: usize = 25;
        Genes {
            map: LinkedHashMap::with_capacity(ngenes),
            index: LinkedHashMap::with_capacity(nchroms),
            names: Vec::with_capacity(ngenes),
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
            names: genes.names,
        };
        genes.build_index();
        genes.build_keys();

        genes
    }
    
    #[inline]
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Draw a gene sample.
    pub fn sample<R: rand::Rng>(&self, rng: &mut R) -> &Gene {
        let range = rand::distributions::Range::new(0, self.len());
        &self.map[&self.names[range.ind_sample(rng)]]
    }

    /// Build vector of gene names.
    fn build_keys(&mut self) {
        for k in self.map.keys() {
            self.names.push(k.clone());
        }
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
        header.extend(mutation::coding::MutOpps::types());

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
    
    pub fn write_counts(&self, snv: &mut SnvReader, mut ifasta: &mut FastaIndexedReader, out: &mut fs::File, aggregate: bool) -> io::Result<()> {
        use std::io::Write;
        
        let sep = "\t";
        
        // mutation type counts map -- key1: sample_id; key2: transcript_id
        let mut muts_map: LinkedHashMap<u32, LinkedHashMap<String, mutation::coding::MutOpps>> = LinkedHashMap::new();
        
        // count mutation types
        for s in snv.records() {
            let record = s.ok().expect("Error reading SNV record");
            print!("{} ... ", record);
            
            let sample_id = if aggregate { 0 } else { record.sample_id };
            let sample_muts_map = muts_map.entry(sample_id).or_insert(LinkedHashMap::new());
            let mut in_coding = false;
            let overlap = self.find_overlap(&record.chrom, record.pos, record.pos + 1);
            for gid in overlap.iter() {
                let gene = self.map.get(gid).unwrap();
                for (tid, transcript) in gene.transcripts.iter() {
                    let cds = get_transcript_cds_from_fasta(&mut ifasta, transcript, &gene.chrom, gene.strand, CDS_PADDING);
                    if let Some(effect) = cds.assign_mutation(gene, transcript, record.pos, record.nt_ref, record.nt_alt) {
                        let t_muts = sample_muts_map.entry(tid.clone()).or_insert(mutation::coding::MutOpps::new());
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
        header.extend(mutation::coding::MutOpps::types());
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
    pub fn write_annotations(&self, snv: &mut SnvReader, mut ifasta: &mut FastaIndexedReader, out: &mut fs::File) -> io::Result<()> {
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

    /// Write annotation to bed file.
    pub fn write_bed(&self, out: &mut fs::File) -> io::Result<()> {
        use std::io::Write;
        
        let sep = "\t";

        // chrom chromStart chromEnd name score strand
        // NB blocks are *not* used
				// coordinates are 0-based, half-open

		let score = "0";
        for (_, gene) in self.map.iter() {
						let name = format!("gene:{}", gene.name);
						let strand = match gene.strand {
								Strand::Forward => "+",
								Strand::Reverse => "-",
								_ => ".",
						};
            let line = vec![gene.chrom.to_owned(), gene.start.to_string(), gene.end.to_string(), name, score.to_owned(), strand.to_owned()].join(sep);
            try!(writeln!(out, "{}", line));
            for (tid, transcript) in gene.transcripts.iter() { 
								let name = format!("transcript:{}", tid);
                let line = vec![gene.chrom.to_owned(), transcript.start.to_string(), transcript.end.to_string(), name, score.to_owned(), strand.to_owned()].join(sep);
                try!(writeln!(out, "{}", line));
                for cds in transcript.coding_regions.iter() {
										let cds_name = format!("CDS:{}", tid);
										let line = vec![gene.chrom.to_owned(), cds.start.to_string(), cds.end.to_string(), cds_name, score.to_owned(), strand.to_owned()].join(sep);
										try!(writeln!(out, "{}", line));
                }
            }
        }

      Ok(())
    }
    
}

#[derive(Debug,Clone)]
struct Snv {
    /// genomic position
    pos: u64,
    /// reference allele
    nt_ref: Nucleotide,
    /// alternative allele
    nt_alt: Nucleotide,
    /// sample ID
    sample_id: u32,
}

/// Collection of single nucleotide variants grouped by chromosome
#[derive(Debug)]
struct GroupedSnvs {
    /// map from chromosome name to sorted vector of snvs
    inner: LinkedHashMap<String, Snvs>,
}

impl GroupedSnvs {
    fn new() -> Self {
        let nchroms: usize = 25;
        GroupedSnvs {
            inner: LinkedHashMap::with_capacity(nchroms),
        }
    }

    fn from_snv(reader: &mut SnvReader) -> Self {
        let mut gsnvs = GroupedSnvs::new();
            
        for r in reader.records() {
            let record = r.ok().expect("Error reading SNV record");
            let snvs = gsnvs.inner.entry(record.chrom.to_owned()).or_insert(Snvs::new());
            snvs.push(Snv {
                pos: record.pos,
                nt_ref: record.nt_ref,
                nt_alt: record.nt_alt,
                sample_id: record.sample_id,
            });
        }
        
        // sort regions within each chromosome by position
        let sorted = gsnvs.inner.into_iter()
            .map(|(chrom, mut snvs)| {
                snvs.sort();
                (chrom, snvs)
            }).collect();

        GroupedSnvs { inner: sorted }
    }
}

#[derive(Debug)]
struct Snvs {
    inner: Vec<Snv>,
}

impl Snvs {
    /// Construct a new, empty SNV collection.
    fn new() -> Self {
        Snvs { inner: Vec::new() }
    }

    /// Append an SNV to the back of the collection.
    fn push(&mut self, snv: Snv) {
        self.inner.push(snv);
    }

    /// Sort SNVs by position.
    fn sort(&mut self) {
        self.inner.sort_by_key(|x| x.pos);
    }

    /// Group SNVs by region.
    /// Return a map from region to Snvs
    fn group_by_regions(&self, regions: &Vec<Region>) -> LinkedHashMap<std::ops::Range<Pos>, Snvs> {
        let mut map = LinkedHashMap::new();

        let opt_indices = self.find_containing_regions(regions);
        let mut sit = self.inner.iter();
        for &x in opt_indices.iter() {
            if let Some(snv) = sit.next() {
                if let Some(idx) = x {
                    let snvs = map.entry(regions[idx].start .. regions[idx].end).or_insert(Snvs::new());
                    snvs.push(snv.clone()); 
                }
            } else {
                // ran out of SNVs
                break;
            }
        }

        map
    }

    /// For each snv, find the region index that contains it, if any.
    /// snvs are sorted by position.
    /// regions are sorted by start position.
    fn find_containing_regions(&self, regions: &Vec<Region>) -> Vec<Option<usize>> {
        let mut sit = self.inner.iter();
        // since SNVs and regions are both sorted, we can shrink the region slice during the search
        let mut slice = &regions[0..];
        let mut result = Vec::with_capacity(self.inner.len());
        loop {
            match sit.next() {
                Some(snv) => {
                    // MAYBE check if current SNV matches previous SNV => skip work
                    // check if first region contains the SNV
                    let region = &slice[0];
                    if snv.pos >= region.start && snv.pos < region.end {
                        result.push(Some(0));
                    } else {
                        // find region containing SNV among remainig regions, if any
                        // binary search will find last of region s.t. query position <= region start
                        let option: Option<usize> = match slice[1..].binary_search_by_key(&snv.pos, |r| r.start) {
                            Err(idx) => {
                                if idx > 0 {
                                    // regions[idx-1].start < snv.pos < regions[idx].start    
                                    let prev_idx = idx - 1;
                                    let region = &slice[prev_idx];
                                    if idx > 1 {
                                        slice = &slice[prev_idx..];
                                    }
                                    if snv.pos < region.end {
                                        // region contains pos
                                        Some(prev_idx)
                                    } else {
                                        // since regions are not overlapping
                                        // no region contains pos
                                        None
                                    }
                                } else {
                                    // else pos < regions[0].start
                                    None
                                }
                            },
                            Ok(idx) => {
                                // pos == regions[idx].start    
                                // region contains pos
                                if idx > 0 {
                                    slice = &slice[idx..];
                                }
                                Some(idx)
                            },
                        };
                        result.push(option);
                    }
                },
                None => break,
            }
        }
        result
    }
}

/// Collection of non-overlapping genomic regions grouped by chromosome
#[derive(Debug)]
struct Regions {
    /// map from chromosome name to sorted vector of regions
    inner: LinkedHashMap<String, Vec<Region>>,
}

impl Regions {
    pub fn new() -> Self {
        let nchroms: usize = 25;
        Regions {
            inner: LinkedHashMap::with_capacity(nchroms),
        }
    }

    pub fn from_bed(reader: &mut BedReader) -> Self {
        let mut regions = Regions::new();

        for r in reader.records() {
            let record = r.expect("Error reading BED record");
            let vec = regions.inner.entry(record.chrom().to_owned()).or_insert(Vec::new());
            // BED file uses 0-based, half-open coordindates: no adjustment required
            vec.push(Region { phase: 0, start: record.start(), end: record.end() } );
        }

        // sort regions within each chromosome by start position
        let sorted = regions.inner.into_iter()
            .map(|(chrom, mut rs)| {
                rs.sort_by_key(|r| r.start);
                (chrom, rs)
            })
            .collect();

        Regions { inner: sorted }
    }

    pub fn write_counts(&self, snv: &mut SnvReader, ifasta: &mut FastaIndexedReader, out: &mut fs::File, aggregate: bool) -> io::Result<()> {
        use std::io::Write;
        
        let sep = "\t";
        
        // mutation channel counts map -- key: sample_id
        let mut muts_map: LinkedHashMap<u32, mutation::genomic::MutOpps> = LinkedHashMap::new();
        
        // both snvs and regions are sorted within each chromosome
        let snvs_map = GroupedSnvs::from_snv(snv);
        let regions_map = &self.inner;

        // get all regions on each chromosome
        for (chrom, snvs) in snvs_map.inner.iter() {
            // get all regions on the same chromosome
            if let Some(regions) = regions_map.get(chrom) {
                // group snvs by region and iterate through regions
                let regioned_snvs = snvs.group_by_regions(regions);
                for (region, snvs) in regioned_snvs.iter() {
                    // get region sequence once for all snvs within region
                    if let Err(why) = ifasta.fetch(chrom, region.start, region.end) {
                        panic!("{:?}", why);
                    }
                    let mut seq: seq::DnaSeq = Vec::new();
                    ifasta.read(&mut seq);
                    let seq = seq::genomic::Sequence{ inner: seq };

                    // count all mutations in region
                    for snv in snvs.inner.iter() {
                        let sample_id = if aggregate { 0 } else { snv.sample_id };
                        let sample_muts = muts_map.entry(sample_id).or_insert(mutation::genomic::MutOpps::new());
                        let region = Region { phase: 0, start: region.start, end: region.end };
                        match seq.assign_channels(&region, snv.pos, snv.nt_ref, snv.nt_alt) {
                            None => {
                                println!("{}:g.{}{}>{} is not in any region or has insufficient context available", chrom, snv.pos + 1, snv.nt_ref as char, snv.nt_alt as char);
                            },
                            Some(idx) => {
                                sample_muts[idx] += 1;
                            },
                        }
                    }
                }
            }
            println!("{} ... done", chrom);
        }
        
        // write header
        let mut header = vec![String::from("sample_id")];
        header.extend(mutation::genomic::MutOpps::channels());
        try!(writeln!(out, "{}", header.join(sep)));

        // write opportunities to file
        for (sid, muts) in muts_map.iter() {
            let mut line = sid.to_string();
            for x in muts.iter() {
                line.push_str(sep);
                line.push_str(&x.to_string());
            }
            try!(writeln!(out, "{}", line));
        }
        
        Ok(())
    }

    pub fn write_opps(&self, ifasta: &mut FastaIndexedReader, out: &mut fs::File) -> io::Result<()> {
        use std::io::Write;

        let sep = "\t";

        let header = mutation::genomic::MutOpps::channels();
        try!(writeln!(out, "{}", header.join(sep)));

        // accmulate mutation opportunities over all regions on all chromosomes
        let mut opp = mutation::genomic::MutOpps::new();
        for (chrom, regions) in self.inner.iter() {
            print!("{} ... ", chrom);
            for region in regions.iter() {
                if let Err(why) = ifasta.fetch(chrom, region.start, region.end) {
                    panic!("{:?}", why);
                }
                let mut seq: seq::DnaSeq = Vec::new();
                ifasta.read(&mut seq);
                let seq = seq::genomic::Sequence{ inner: seq };
                seq.accumulate_opp(&mut opp);
            }
            println!("done");
        }

        // write results to file
        let mut line = String::new();
        for x in opp.iter() {
            if !line.is_empty() {
                line.push_str(sep);
            }
            line.push_str(&x.to_string());
        }
        try!(writeln!(out, "{}", line));

        Ok(())
    }
}

struct MutSpecSnv {
    /// chromosome
    chrom: String,
    /// genomic position
    pos: u64,
    /// reference nucleotide
    nt_ref: Nucleotide,
    /// alternative nucleotide
    nt_alt: Nucleotide,
    /// 5' nucleotide
    nt_5p: Nucleotide,
    /// 3'nucleotide
    nt_3p: Nucleotide,
    /// channel index
    channel: usize,
}

fn write_snv_sample(mut ifasta: &mut FastaIndexedReader, genes: &Genes, mutspec: &mutation::genomic::MutSpec, out: &mut fs::File, nsamples: usize, mut out_imp_opt: Option<&mut fs::File>) -> io::Result<()> {
    use std::io::Write;
    let sep = "\t";

    if let Some(ref mut out_imp) = out_imp_opt {
        let header = vec![
            String::from("chrom"),
            String::from("pos"),
            String::from("ref"),
            String::from("alt"),
            String::from("context"),
            String::from("lweight")
        ];

        try!(writeln!(out_imp, "{}", header.join(sep)));
    }

    let ngenes = genes.len();

    //let batch_size = 1e6;

    let mut rng = rand::thread_rng();

    let size_factor = 10;

    // number of importance samples
    let nsamples_imp = (nsamples * size_factor) as usize;

    let mut snvs: Vec<MutSpecSnv> = Vec::with_capacity(nsamples_imp);
    let mut weights: Vec<f64> = Vec::with_capacity(nsamples_imp);

    // draw the importance samples
    for _ in 0..nsamples_imp {
        // sample a gene uniformly
        let gene = genes.sample(&mut rng);
        //println!("Sampled gene: {}", gene.name);

        // select the canonical transcript
        let transcript;
        match gene.canonical_transcript() {
            Some(x) => {
                transcript = x;
            },
            None => {
                break;
            },
        }

        // sample a region from the transcript uniformly
        let nregions = transcript.coding_regions.len();
        let r = rand::distributions::Range::new(
            0, nregions
        ).ind_sample(&mut rng);
        let region = &transcript.coding_regions[r];
        //println!("Sampled region: {:?}", region);

        // TODO allow sampling of splice sites

        // sample a position from the region uniformly
        let nbases = region.end - region.start;
        let n = rand::distributions::Range::new(
            0, nbases
        ).ind_sample(&mut rng);
        let position = region.start + n;
        //println!("Sampled position: {}", position);

        // read reference nucleotide at position
        if let Err(why) = ifasta.fetch(&gene.chrom, position - 1, position + 2) {
            panic!("{:?}", why);
        }
        let mut seq: seq::DnaSeq = vec![b'N'; 1];
        ifasta.read(&mut seq);

        let nt_5p = seq[0];
        let nt_ref = seq[1];
        let nt_3p = seq[2];
        //println!("Sequence: {}", str::from_utf8(&seq).unwrap());

        // sample an alternative nucleotide uniformly
        let a = rand::distributions::Range::new(
            0, 3
        ).ind_sample(&mut rng);
        let nt_alt = match nt_ref {
            b'A' => match a {
                0 => b'C',
                1 => b'G',
                2 => b'T',
                _ => b'N',
            },
            b'C' => match a {
                0 => b'A',
                1 => b'G',
                2 => b'T',
                _ => b'N',
            },
            b'G' => match a {
                0 => b'A',
                1 => b'C',
                2 => b'T',
                _ => b'N',
            },
            b'T' =>  match a {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'N',
            },
            _ => b'N',
        };
        //println!("nt_alt: {}", str::from_utf8(&[nt_alt]).unwrap());

        // w(x) = p(x) / q(x)
        // where p(x) is the target distribution
        //       q(x) is the sampling distribution
        let q = (1.0 / ngenes as f64) * (1.0 / nregions as f64) * 
            (1.0 / nbases as f64) * (1.0 / 3.0);
        // p(x) is defined by the input mutation spectrum
        let channel = mutation::genomic::MutSpec::index(nt_ref, nt_alt, nt_5p, nt_3p);
        let p = mutspec[channel];
        let lweight = (p / q).ln();
        //println!("p: {}", p);
        //println!("q: {}", q);
        //println!("weight: {}", weight);
        
        snvs.push(
            MutSpecSnv {
                chrom: gene.chrom.clone(),
                pos: position,
                nt_ref: nt_ref,
                nt_alt: nt_alt,
                nt_5p: nt_5p,
                nt_3p: nt_3p,
                channel: channel,
            }
        );

        weights.push(lweight);

        if let Some(ref mut out_imp) = out_imp_opt {
            let mut context = seq.clone();
            if nt_ref != b'C' && nt_ref  != b'T' {
                seq::reverse_complement(&mut context);
            }

            let line = vec![
                gene.chrom.clone(),
                // convert from 0-based to 1-based
                (position + 1).to_string(),
                str::from_utf8(&[nt_ref]).unwrap().to_owned(),
                str::from_utf8(&[nt_alt]).unwrap().to_owned(),
                str::from_utf8(&context).unwrap().to_owned(),
                lweight.to_string(),
            ].join(sep);

            try!(writeln!(out_imp, "{}", line));
        }
    }

    
    // normalize the log importance weights and anti-log transform
    let lweight_sum = stats::log_sum_exp(&weights);
    for w in weights.iter_mut() {
        *w = (*w - lweight_sum).exp();
    }

    // weights now contain sampling probabilities,
    // which we use to resample the SNVs
    let sample_idx = sample::sample(&weights, nsamples, false, &mut rng);

    println!("Jensen-Shannon Divergence (before rejection): {}", evaluate_sample(&sample_idx, &snvs, &mutspec.0));

    // refine the sample using rejection sampling
    let channel_values: Vec<usize> = sample_idx.iter().map(|&i| snvs[i].channel).collect();
    let accept_idx = sample::resample_factor_distrib(&channel_values, &mutspec.0, 0.01, &mut rng);

    println!("Accepted {} samples after rejection sampling", accept_idx.len());

    let sample_idx: Vec<usize> = accept_idx.iter().map(|&a| sample_idx[a]).collect();

    println!("Jensen-Shannon Divergence (after rejection): {}", evaluate_sample(&sample_idx, &snvs, &mutspec.0));

    // write resampled SNVs to file 

    let header = vec![
        String::from("chrom"),
        String::from("pos"),
        String::from("ref"),
        String::from("alt"),
        String::from("context"),
    ];

    try!(writeln!(out, "{}", header.join(sep)));

    for &i in sample_idx.iter() {
        let MutSpecSnv { ref chrom, pos: position, nt_ref, nt_alt, nt_5p, nt_3p, channel: _ } = snvs[i];

        let mut context = vec![nt_5p, nt_ref, nt_3p];
        if nt_ref != b'C' && nt_ref  != b'T' {
            seq::reverse_complement(&mut context);
        }

        let line = vec![
            chrom.clone(),
            // convert from 0-based to 1-based
            (position + 1).to_string(),
            str::from_utf8(&[nt_ref]).unwrap().to_owned(),
            str::from_utf8(&[nt_alt]).unwrap().to_owned(),
            str::from_utf8(&context).unwrap().to_owned(),
        ].join(sep);


        try!(writeln!(out, "{}", line));
    }

    Ok(())
}

fn evaluate_sample(sample_idx: &[usize], snvs: &[MutSpecSnv], target: &[f64]) -> f64 {
    // count mutations in each channels for the resampled data
    let mut counts = vec![0 as u64; mutation::genomic::MutSpec::len()];
    for &i in sample_idx.iter() {
        counts[ snvs[i].channel ] += 1;
    }
    let mut counts_sum = 0.0;
    for &c in counts.iter() {
        counts_sum += c as f64;
    }

    // observed mutation spectrum in the resampled data
    let mutspec_sample: Vec<f64> = counts.iter().map(|&x| (x as f64) / counts_sum).collect();

    stats::js_divergence(target, &mutspec_sample)
}

