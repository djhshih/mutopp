// Number of non-degenerate DNA nucleotides
pub const N_NUCLEOTIDES: usize = 4;

// Collection of non-degenerate DNA nucleotides
pub const NUCLEOTIDES: [u8; N_NUCLEOTIDES] = [b'A', b'C', b'G', b'T'];

// Number of valid mutation impact classes
pub const N_MUTATION_CLASSES: usize = 4;

// Number of substitution types, considering strand
pub const N_STRANDED_SUBSTITUTION_TYPES: usize = N_NUCLEOTIDES * (N_NUCLEOTIDES - 1);

// Number of mutation channels, considering strand, 5' one-nucleotide context, and 3' one-nucleotide context
pub const N_STRANDED_MUTATION_CHANNELS: usize = N_STRANDED_SUBSTITUTION_TYPES * N_NUCLEOTIDES * N_NUCLEOTIDES;

// Number of mutation types, considering mutation impact and mutation channel
pub const N_STRANDED_MUTATION_TYPES: usize = N_MUTATION_CLASSES * N_STRANDED_MUTATION_CHANNELS;

// Number of substitution types, ignoring strand
pub const N_NONSTRANDED_SUBSTITUTION_TYPES: usize = N_STRANDED_SUBSTITUTION_TYPES / 2;

// Number of mutation channels, ignoring strand, considering 5' one-nucleotide context, and 3' one-nucleotide context
pub const N_NONSTRANDED_MUTATION_CHANNELS: usize = N_NONSTRANDED_SUBSTITUTION_TYPES * N_NUCLEOTIDES * N_NUCLEOTIDES;