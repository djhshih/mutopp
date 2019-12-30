# mutopp

Tool for categorizing mutations and enumerating mutation opportunities.

## Dependencies

- [rust](https://www.rust-lang.org/en-US/install.html)

## Installation

```
$ cargo install
```

## Usage

See a list of available commands using

```
$ mutopp help
```

```
USAGE:
    mutopp [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    annot-gene-muts      annotate observed mutations
    count-gene-muts      count mutation types of observed mutations
    count-region-muts    count mutation channels of observed mutations
    enum-genes           enumerate all mutation opportunities in genes
    enum-regions         enumerate all mutation opportunities
    extract-fasta        extract sequence from .fasta file
    gff-to-bed           convert transcript information from GFF3 to BED
                         format
    help                 Prints this message or the help of the given
                         subcommand(s)
```

