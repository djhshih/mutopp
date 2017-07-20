pub const nucleotides: [u8; 4] = [b'A', b'C', b'G', b'T'];
pub const n_nucleotides: usize = 4;

pub const n_substitution_types: usize = n_nucleotides * (n_nucleotides - 1);

pub const n_mutation_channels: usize = n_substitution_types * n_nucleotides * n_nucleotides;

pub const n_mutation_classes: usize = 4;
pub const n_mutation_types: usize = (n_mutation_classes * n_mutation_channels) as usize;