#![feature(is_sorted)]
///! ORF finding
///
/// Algorithm:
///     Given a sequence, find stop codons (use jetscii)
///     Calculate distance between stop codons
///     If greater than MIN_ORF_LENGTH, translate intervening sequence
use lazy_static::lazy_static;

lazy_static! {
    static ref STOP_CODONS: [ByteSubstringConst; 3] = [
        ByteSubstring::new(b"TAG"),
        ByteSubstring::new(b"TAA"),
        ByteSubstring::new(b"TGA")
    ];
}

pub mod helpful;
pub use helpful::*;

use jetscii::{ByteSubstring, ByteSubstringConst};
use pulp::Arch;

pub fn find_stop_codons(sequence: &[u8]) -> Vec<usize> {
    let mut sequence = sequence.to_vec();
    sequence.make_ascii_uppercase();

    let mut found_stop_codons = Vec::new();

    for stop_codon in STOP_CODONS.as_slice().iter() {
        let mut offset = 0;
        while let Some(pos) = stop_codon.find(&sequence[offset..]) {
            println!("Pos: {}", pos + offset);
            found_stop_codons.push(pos + offset);
            offset += pos + 3;
        }

        if offset >= sequence.len() {
            break;
        }
    }

    found_stop_codons.sort_unstable();
    found_stop_codons
}

/// Returns reading_frame, start, end
pub fn stop_codons_to_intervals(
    stop_codons: &[usize],
    min_orf_length: usize,
    sequence_length: usize,
) -> Vec<(usize, usize, usize)> {
    let mut intervals = Vec::new();
    let mut reading_frames_and_stops = Vec::new();

    // Start of the sequence and end of the sequence...
    for i in 0..3 {
        reading_frames_and_stops.push((i, 0));
        if i == 0 {
            let end = sequence_length - ((sequence_length - i) % 3);
            reading_frames_and_stops.push((i, end));
        } else {
            let end = sequence_length - ((sequence_length - i) % 3);
            reading_frames_and_stops.push((i, end));
        }
    }

    for stop in stop_codons.iter() {
        let reading_frame = stop % 3;
        reading_frames_and_stops.push((reading_frame, stop + 3));
    }

    reading_frames_and_stops.sort_by_key(|x| (x.0, x.1));

    for x in reading_frames_and_stops.windows(2) {
        let (reading_frame, start) = x[0];
        let (reading_frame2, end) = x[1];
        if reading_frame == reading_frame2 {
            let length = end - start;
            if length >= min_orf_length {
                intervals.push((reading_frame, start, end));
            }
        }
    }

    intervals
}

pub fn translate_sequence(sequence: &[u8]) -> Vec<Amino> {
    let mut result = Vec::new();

    for codon in sequence.chunks_exact(3) {
        //let codon: &[u8] = &sequence[i..i + 3];
        println!("Codon: {:?}", std::str::from_utf8(codon).unwrap());
        let amino = AMINO_MAPPING[CODON_MAPPING
            .binary_search_by(|&c| c.as_slice().cmp(codon))
            .unwrap()];
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

pub fn complement(c: &mut u8) {
    if *c != b'N' {
        if *c & 2 != 0 {
            *c ^= 4;
        } else {
            *c ^= 21;
        }
    }
}

pub fn revcomp(sequence: &mut [u8]) {
    let arch = Arch::new();
    arch.dispatch(|| {
        sequence.reverse();
        sequence.make_ascii_uppercase();
        sequence.iter_mut().for_each(complement);
    });
}

pub fn translate_interval(
    sequence: &[u8],
    (_, start, stop): &(usize, usize, usize)
) -> Vec<Amino> {
    let seq_to_translate = &sequence[*start..*stop];
    let mut translated = translate_sequence(seq_to_translate);
    translated
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
        assert_eq!(
            result,
            vec![
                Amino::M,
                Amino::R,
                Amino::Y,
                Amino::A,
                Amino::Y,
                Amino::A,
                Amino::I,
                Amino::R
            ]
        );
    }

    #[test]
    fn test_find_stops() {
        let sequence = b"ccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgataaaagcgtaccgtccccggcTCCTGTGAACTCCCTCCCCCCTTTTACCGCTTGATTGTACTAGAGTCGGTTTGTAGTCATCTTCCAGTTCTAAGCGTGTCACTCATGTGAAGTTAGTAGTTTCGTTTGTAGATCTGGCTAACTAAGGtgattaagttttatataattttttttaagttttgctaaaaatcattttagaagatattttttaaaaatttatgttcttttatgtggtcctttctcaaaatatattgtactgtatatttattataataaagtaccgtatttagttattaaaaatcagCTTGATCGTGTAATAAaaacacaggaaaaaaataaaatttatacaacaaattgtgaaatattaatattacaatgataaaaaataaagttatgaattaaaaatagccTAGGGCCTAGGCTATTAGAAATTGACTTACACATTTGATAACGTGACTCACTATAATGatagattttaatgtattataaaataatttggcaaCCAAAGTACGATTTTATCCATGTTTATGCATACCGCGTCTGACTATAAATCGTTGCATGTGTGGAGTACGTTTTCTTTGATCAGGTCTGTCAAACTTCATCGACAATTTAGCTTTAGACCGGTCTTCAGATAAATTGTTGCATTATGTGGGGGCGTAAGAGAACAGTAAGAACTCAATAATctaatagcctaatttttatttttggacatTTCAAAGCACTCTAAGAGAAAAGTGACACTCAATTAGCTGTAGGATATATGAGTACTGATATAGTTTTTTATGACCCTGttttagactttgaaattttctttcattgtcatTATTTACGGGATTTGTACAAACTGTGCTGGTGCAAGGCTATTTCAGTAACGAGTGTGTCTGATTTATGGCCCAGGAAATTTAAGACTATCTCAATTCTGAATTGAATGGGTTAACAAAGGCAGCAGTTATAAcctgaaataaacaattttaaaataacacaaatgctgtccaaaatattttttaattttataaaatgttatttagccTAATACCAGTAGTAGGTGACCCTAGTACAGCCATCAGTAGCCTATGATTCAGCAATATGATAAGATACAACAAATGTTTAATTGTATTCAAACTAAATATGTAAAaccttgaacattttttttgtcatatacaataaatactattAACAAAATACATTCCTCCGTGCCCCCTCCTCTCCCCTGGAATTAATAATTGctgattttgatttaaattaattttctgatttaagtcatatctttaattataattttgggtCAAATTATCTCACAAACCATgcataatcatattaaaaaaaatgcgctgtatttgtacttttaaaaaaaaatcactactcGTGAGAATCATGAGCAAAATATTCTAAAGTGGAAACGGCACTAAGGTGAACTAAGCAACTTAGTGCAAAActaaatgttagaaaaaatatCCTACACTGCATAAACTATTTTGcaccataaaaaaaagttatgtgtgGGTCTAAAATAATTTGCTGAGCAATTAATGATTTCTAAATGATGCTAAAGTGAACCATTGTAatgttatatgaaaaataaatacacaattaagATCAACACAGTGAAATAACATTGATTGGGTGATTTCAAATGGGGTCTATctgaataatgttttatttaacagtaatttttatttctatcaatttttagtaatatctacaaatattttgttttaggcTGCCAGAAGATCGGCGGTGCAAGGTCAGAGGTGAGATGTTAGGTGGTTCCACCAACTGCACGGAAGAGCTGCCCTCTGTCATTCAAAATTTGACAGGTACAAACAGactatattaaataagaaaaacaaactttttaaaggCTTGACCATTAGTGAATAGGTTATATGCTTATTATTTCCATTTAGCTTTTTGAGACTAGTATGATTAGACAAATCTGCTTAGttcattttcatataatattgaGGAACAAAATTTGTGAGATTTTGCTAAAATAACTTGCTTTGCTTGTTTATAGAGGCacagtaaatcttttttattattattataattttagattttttaatttttaaataagtgataGCATAtgctgtattattaaaaatttaagaactttaaagtatttacagtagcctatattatgtaaTAGGCTAGCCTACTTTATTGTTCGGccaattctttttcttattcatCTTATCATTATTagcagattattatttattactagtttaaaagcacgtcaaaaatgacggacattaatttcttcctccttaggttgctatacttaaatgccggtcaaaaatttaaaagcccgtcaaaaataacagacagcgacatctatggacagaaaatagagagaaacaTCTTCGGGCaacatcggctcgccgaagatggccgagcagtatgacagacatcatgaccaaactttctctttattatagttagattagatagaagattagaagactagtttaaaagcccatca";
        let result = find_stop_codons(sequence);

        let results = stop_codons_to_intervals(&result, 1, sequence.len());
        assert!(results.contains(&(0, 0, 9)));
        assert!(results.contains(&(0, 78, 186)));
        assert!(results.contains(&(0, 2550, 2559)));
        assert!(results.contains(&(1, 547, 619)));
        assert!(results.contains(&(1, 2545, 2560)));
        assert!(results.contains(&(2, 2537, 2558)));
    }

    #[test]
    fn test_revcomp() {
        let mut sequence = b"ACGTCACTACGATCGATCAGACTNTNNATCGATCATCA".to_vec();
        let rev = b"TGATGATCGATNNANAGTCTGATCGATCGTAGTGACGT";
        revcomp(sequence.as_mut_slice());
        assert_eq!(sequence, rev);
    }

    #[test]
    fn test_translate_interval() {
        let sequence = b"ccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgataaaagcgtaccgtccccggcTCCTGTGAACTCCCTCCCCCCTTTTACCGCTTGATTGTACTAGAGTCGGTTTGTAGTCATCTTCCAGTTCTAAGCGTGTCACTCATGTGAAGTTAGTAGTTTCGTTTGTAGATCTGGCTAACTAAGGtgattaagttttatataattttttttaagttttgctaaaaatcattttagaagatattttttaaaaatttatgttcttttatgtggtcctttctcaaaatatattgtactgtatatttattataataaagtaccgtatttagttattaaaaatcagCTTGATCGTGTAATAAaaacacaggaaaaaaataaaatttatacaacaaattgtgaaatattaatattacaatgataaaaaataaagttatgaattaaaaatagccTAGGGCCTAGGCTATTAGAAATTGACTTACACATTTGATAACGTGACTCACTATAATGatagattttaatgtattataaaataatttggcaaCCAAAGTACGATTTTATCCATGTTTATGCATACCGCGTCTGACTATAAATCGTTGCATGTGTGGAGTACGTTTTCTTTGATCAGGTCTGTCAAACTTCATCGACAATTTAGCTTTAGACCGGTCTTCAGATAAATTGTTGCATTATGTGGGGGCGTAAGAGAACAGTAAGAACTCAATAATctaatagcctaatttttatttttggacatTTCAAAGCACTCTAAGAGAAAAGTGACACTCAATTAGCTGTAGGATATATGAGTACTGATATAGTTTTTTATGACCCTGttttagactttgaaattttctttcattgtcatTATTTACGGGATTTGTACAAACTGTGCTGGTGCAAGGCTATTTCAGTAACGAGTGTGTCTGATTTATGGCCCAGGAAATTTAAGACTATCTCAATTCTGAATTGAATGGGTTAACAAAGGCAGCAGTTATAAcctgaaataaacaattttaaaataacacaaatgctgtccaaaatattttttaattttataaaatgttatttagccTAATACCAGTAGTAGGTGACCCTAGTACAGCCATCAGTAGCCTATGATTCAGCAATATGATAAGATACAACAAATGTTTAATTGTATTCAAACTAAATATGTAAAaccttgaacattttttttgtcatatacaataaatactattAACAAAATACATTCCTCCGTGCCCCCTCCTCTCCCCTGGAATTAATAATTGctgattttgatttaaattaattttctgatttaagtcatatctttaattataattttgggtCAAATTATCTCACAAACCATgcataatcatattaaaaaaaatgcgctgtatttgtacttttaaaaaaaaatcactactcGTGAGAATCATGAGCAAAATATTCTAAAGTGGAAACGGCACTAAGGTGAACTAAGCAACTTAGTGCAAAActaaatgttagaaaaaatatCCTACACTGCATAAACTATTTTGcaccataaaaaaaagttatgtgtgGGTCTAAAATAATTTGCTGAGCAATTAATGATTTCTAAATGATGCTAAAGTGAACCATTGTAatgttatatgaaaaataaatacacaattaagATCAACACAGTGAAATAACATTGATTGGGTGATTTCAAATGGGGTCTATctgaataatgttttatttaacagtaatttttatttctatcaatttttagtaatatctacaaatattttgttttaggcTGCCAGAAGATCGGCGGTGCAAGGTCAGAGGTGAGATGTTAGGTGGTTCCACCAACTGCACGGAAGAGCTGCCCTCTGTCATTCAAAATTTGACAGGTACAAACAGactatattaaataagaaaaacaaactttttaaaggCTTGACCATTAGTGAATAGGTTATATGCTTATTATTTCCATTTAGCTTTTTGAGACTAGTATGATTAGACAAATCTGCTTAGttcattttcatataatattgaGGAACAAAATTTGTGAGATTTTGCTAAAATAACTTGCTTTGCTTGTTTATAGAGGCacagtaaatcttttttattattattataattttagattttttaatttttaaataagtgataGCATAtgctgtattattaaaaatttaagaactttaaagtatttacagtagcctatattatgtaaTAGGCTAGCCTACTTTATTGTTCGGccaattctttttcttattcatCTTATCATTATTagcagattattatttattactagtttaaaagcacgtcaaaaatgacggacattaatttcttcctccttaggttgctatacttaaatgccggtcaaaaatttaaaagcccgtcaaaaataacagacagcgacatctatggacagaaaatagagagaaacaTCTTCGGGCaacatcggctcgccgaagatggccgagcagtatgacagacatcatgaccaaactttctctttattatagttagattagatagaagattagaagactagtttaaaagcccatca";
        let mut sequence = sequence.to_vec();
        sequence.make_ascii_uppercase();
        let translated = translate_interval(&sequence, &(2, 2537, 2558));
        assert_eq!(translated, vec![
            Amino::K,
            Amino::T,
            Amino::S,
            Amino::L,
            Amino::K,
            Amino::A,
            Amino::H,
        ])



    }
}
