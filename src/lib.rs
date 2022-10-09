#![feature(is_sorted)]
///! ORF finding
/// 
/// Algorithm:
///     Given a sequence, find stop codons (use jetscii)
///     Calculate distance between stop codons
///     If greater than MIN_ORF_LENGTH, translate intervening sequence

use lazy_static::lazy_static;

use jetscii::{SubstringConst, Substring};
use pulp::Arch;


lazy_static! {
    static ref STOP_CODONS: [SubstringConst; 3] = [Substring::new("TAG"), Substring::new("TAA"), Substring::new("TGA")];
}

const CODON_MAPPING: [&[u8; 3]; 64] =  [
    // Make it sorted so we can use binary search
    b"AAA", b"AAC", b"AAG", b"AAT", b"ACA", b"ACC", b"ACG", b"ACT", b"AGA", b"AGC", b"AGG", b"AGT", b"ATA", b"ATC", b"ATG", b"ATT",
    b"CAA", b"CAC", b"CAG", b"CAT", b"CCA", b"CCC", b"CCG", b"CCT", b"CGA", b"CGC", b"CGG", b"CGT", b"CTA", b"CTC", b"CTG", b"CTT",
    b"GAA", b"GAC", b"GAG", b"GAT", b"GCA", b"GCC", b"GCG", b"GCT", b"GGA", b"GGC", b"GGG", b"GGT", b"GTA", b"GTC", b"GTG", b"GTT",
    b"TAA", b"TAC", b"TAG", b"TAT", b"TCA", b"TCC", b"TCG", b"TCT", b"TGA", b"TGC", b"TGG", b"TGT", b"TTA", b"TTC", b"TTG", b"TTT",
];

// Map from codon to amino acid
const AMINO_MAPPING: [u8; 64] = [
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

pub fn find_stop_codon_intervals(sequence: &[u8]) -> Vec<usize> {
    let mut stop_codon_intervals = Vec::new();
    for reading_frame in 0..3 {
        for stop_codon in STOP_CODONS.iter() {

        }
    }

    stop_codon_intervals
}

pub fn translate_sequence(sequence: &[u8]) -> Vec<Amino> {
    let mut result = Vec::new();
    for codon in sequence.chunks_exact(3) {
        //let codon: &[u8] = &sequence[i..i + 3];
        let amino = AMINO_MAPPING[CODON_MAPPING.binary_search_by(|&c| c.as_slice().cmp(codon)).unwrap()];
        let amino = match amino {
            b'A' => Amino::A,
            b'C' => Amino::C,
            b'D' => Amino::D,
            b'E' => Amino::E,
            b'F' => Amino::F,
            b'G' => Amino::G,
            b'H' => Amino::H,
            b'I' => Amino::I,
            b'K' => Amino::K,
            b'L' => Amino::L,
            b'M' => Amino::M,
            b'N' => Amino::N,
            b'P' => Amino::P,
            b'Q' => Amino::Q,
            b'R' => Amino::R,
            b'S' => Amino::S,
            b'T' => Amino::T,
            b'V' => Amino::V,
            b'W' => Amino::W,
            b'Y' => Amino::Y,
            b'*' => Amino::Stop,
            _ => unreachable!(),
        };
        result.push(amino);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_confirm_is_sorted() {
        assert!(CODON_MAPPING.is_sorted());
    }

    #[test]
    fn sizeof_aminoacid_enum() {
        assert_eq!(std::mem::size_of::<Amino>(), 1);
    }

    #[test]
    fn test_translate_sequence() {
        let sequence = b"ATGCGA";
        let result = translate_sequence(sequence);
        assert_eq!(result, vec![Amino::M, Amino::R]);

        let sequence = b"ATGCGATACGCTT";
        let result = translate_sequence(sequence);
        assert_eq!(result, vec![Amino::M, Amino::R, Amino::Y, Amino::A]);

        // Test a sequence that is not an ORF
        let sequence = b"ATGCGATACGCTTATGCGATACGCTT";
        let result = translate_sequence(sequence);
        assert_eq!(result, vec![Amino::M, Amino::R, Amino::Y, Amino::A, Amino::Y, Amino::A, Amino::I, Amino::R]);

    }
}
