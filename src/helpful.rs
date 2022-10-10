pub const CODON_MAPPING: [&[u8; 3]; 64] = [
    // Make it sorted so we can use binary search
    b"AAA", b"AAC", b"AAG", b"AAT", b"ACA", b"ACC", b"ACG", b"ACT", b"AGA", b"AGC", b"AGG", b"AGT",
    b"ATA", b"ATC", b"ATG", b"ATT", b"CAA", b"CAC", b"CAG", b"CAT", b"CCA", b"CCC", b"CCG", b"CCT",
    b"CGA", b"CGC", b"CGG", b"CGT", b"CTA", b"CTC", b"CTG", b"CTT", b"GAA", b"GAC", b"GAG", b"GAT",
    b"GCA", b"GCC", b"GCG", b"GCT", b"GGA", b"GGC", b"GGG", b"GGT", b"GTA", b"GTC", b"GTG", b"GTT",
    b"TAA", b"TAC", b"TAG", b"TAT", b"TCA", b"TCC", b"TCG", b"TCT", b"TGA", b"TGC", b"TGG", b"TGT",
    b"TTA", b"TTC", b"TTG", b"TTT",
];

// Map from codon to amino acid
pub const AMINO_MAPPING: [u8; 64] = [
    // Mapping for a sorted codon array
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S', b'I', b'I', b'M', b'I',
    b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P', b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L',
    b'E', b'D', b'E', b'D', b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'*', b'C', b'W', b'C', b'L', b'F', b'L', b'F',
];

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Amino {
    A,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    K,
    L,
    M,
    N,
    P,
    Q,
    R,
    S,
    T,
    V,
    W,
    Y,
    Stop,
}
