extern crate bio;
extern crate multimap;
extern crate mutopp;

use std::str;
use std::env;
use std::process;
use std::fs::File;

use multimap::MultiMap;

use bio::io::fasta;
use bio::io::gff;

use mutopp::gene::{Gene,Transcript,Region,HashMap,Strand};
use mutopp::seq::{self,CodingSequence};
use mutopp::mutation::MutOpps;

type FastaIndexedReader = fasta::IndexedReader<std::fs::File>;
type GffReader = gff::Reader<std::fs::File>;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 4 {
        println!("usage: mutopp <fasta> <gff3> <out>");
        process::exit(1);
    }

    let fasta_fn = &args[1];
    let gff_fn = &args[2];
    let out_fn = &args[3];

    // read the entire fasta
    /*
    let fasta = match fasta::Reader::from_file(fasta_fn) {
        Err(why) => panic!("{:?}", why),
        Ok(reader) => reader,
    };

    for (_, r) in fasta.records().enumerate() {
        let record = r.expect("Error reading FASTA record");
        let desc = match record.desc() {
            Some(x) => x,
            None  => "",
        };

        println!("{}", record.id());
        println!("{}", desc);
        print_cds(record.seq());

        /*
        (record.seq());

        let cds = CodingSequence { seqs: vec![record.seq()], padding: 0 };
        let opp = cds.count_opp().unwrap();
        print_opp_matrix(&opp);
        */

        println!("");
    }
    */

    // extract gene annotation from gff into struct Gene

    let mut gff = match gff::Reader::from_file(gff_fn, gff::GffType::GFF3) {
        Err(why) => panic!("{:?}", why),
        Ok(reader) => reader,
    };

    let genes = genes_from_gff(&mut gff);
    println!("number of valid protein-coding genes: {}", genes.len());
    //println!("{:?}", genes);

    let mut ifasta = match fasta::IndexedReader::from_file(fasta_fn) {
        Err(why) => panic!("{:?}", why),
        Ok(reader) => reader,
    };

    /*
    let coord = &args[2];
    // read only part of the fasta
    if let Some((seqname, start, end)) = parse_coordinate(coord) {
        println!("{} {} {}", seqname, start, end);
        let mut seq: Vec<u8> = Vec::new();
        match ifasta.read(&seqname, start, end, &mut seq) {
            Err(why) => panic!("{:?}", why),
            Ok(()) => print_seq(&seq),
        }
    } else {
        panic!("Invalid coordinate: {}", coord);
    }
    */

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut out_file = match File::create(&out_fn) {
        Err(why) => panic!("Could not create {}: {:?}",
                           out_fn,
                           why),
        Ok(file) => file,
    };

    write_opp(genes, &mut ifasta, &mut out_file);

    /*
    let t = Transcript {
       start: 0,
       end: 1403,
       coding_regions: vec![Region{start: 4, end: 10}, Region{start: 15, end: 18}],
    };
    */


    /*
    print_mutation_types();

    print_mutation_channels();

    let xs = mutation_channels();
    for (i, x) in xs.iter().enumerate() {
        println!("{} {}", i, x);
    }
    */
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

/// Construct genes from GFF reader.
///
/// GFF coordinates are 1-based, closed. Convert the coordinates to 0-based, half-open.
fn genes_from_gff(reader: &mut GffReader) -> HashMap<String, Gene> {
    let mut genes: HashMap<String, Gene> = HashMap::new();

    for r in reader.records() {
        let record = r.expect("Error reading GFF3 record");
        let attrs = &record.attributes();
        match record.feature_type() {
            "gene" => {
                if attribute_filter(attrs, "gene_type", "protein_coding") {
                    genes.insert(
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
                if attribute_filter(attrs, "transcript_support_level", "1") && attribute_filter(attrs, "transcript_type", "protein_coding") {
                    // assume that gene, if valid, has already been inserted
                    if let Some(gene) = genes.get_mut(attrs.get("gene_id").unwrap()) {
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
                if let Some(gene) = genes.get_mut(attrs.get("gene_id").unwrap()) {
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
    genes.into_iter()
        .filter(|&(_, ref v)| !v.transcripts.is_empty())
        .collect()
}

/// Parse genomic coordinate.
///
/// Assume coordinate is specified in 1-based, closed format.
/// Internally, coordinate is converted to 0-based, half-open format.
///
fn parse_coordinate(x: &str) -> Option<(String, u64, u64)> {
    match x.find(':') {
        Some(i) => {
            let seqname = &x[..i];
            let range = &x[i+1..];
            match range.find('-') {
                Some(j) => {
                    let start = &range[..j];
                    let end = &range[j+1..];
                    if let Ok(start) = start.parse() {
                        if let Ok(end) = end.parse() {
                            if start > 0 && end > start {
                                Some((String::from(seqname), start-1, end))
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                },
                None => None,
            }
        },
        None => None,
    }
}

fn write_opp(genes: HashMap<String, Gene>, mut ifasta: &mut FastaIndexedReader, out: &mut File) {
    use std::io::Write;
    use std::error::Error;

    let sep = "\t";

    let mut header = vec![String::from("gene"), String::from("symbol"), String::from("transcript")];
    header.extend(MutOpps::types());

    if let Err(why) = writeln!(out, "{}", header.join(sep)) {
        panic!("couldn't write to output file: {}", why.description())
    }

    for (gid, gene) in genes.iter() {
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
                writeln!(out, "{}", line);
                println!("succeeded");
            } else {
                println!("failed");
            }
        }
    }

}
