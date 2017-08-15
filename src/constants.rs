pub const nucleotides: [u8; 4] = [b'A', b'C', b'G', b'T'];
pub const n_nucleotides: usize = 4;

// TODO rename as N_STRANDED_SUBSTITUTION_TYPES
pub const n_substitution_types: usize = n_nucleotides * (n_nucleotides - 1);

// TODO rename as N_STRANDED_MUTATION_CHANNELS
pub const n_mutation_channels: usize = n_substitution_types * n_nucleotides * n_nucleotides;

pub const n_mutation_classes: usize = 4;
pub const n_mutation_types: usize = (n_mutation_classes * n_mutation_channels) as usize;

pub const N_NONSTRANDED_SUBSTITUTION_TYPES: usize = n_substitution_types / 2;

pub const N_NONSTRANDED_MUTATION_CHANNELS: usize = N_NONSTRANDED_SUBSTITUTION_TYPES * n_nucleotides * n_nucleotides;