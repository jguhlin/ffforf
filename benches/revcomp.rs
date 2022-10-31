// TODO

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pulp::Arch;

use ffforf::*;

/// Modified input in place
pub fn revcomp_pulp(sequence: &mut [u8]) {
    let arch = Arch::new();
    arch.dispatch(|| {
        sequence.reverse();
        sequence.make_ascii_uppercase();
        for i in 0..sequence.len() {
            sequence[i] = match sequence[i] {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            };
        }
    });
}

/// Modified input in place
pub fn revcomp(sequence: &mut [u8]) {
    sequence.reverse();
    sequence.make_ascii_uppercase();
    for i in 0..sequence.len() {
        sequence[i] = match sequence[i] {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        };
    }
}

pub fn complement(c: &mut u8) {
    let val = *c;
    let new_val = if val != b'N' {
        if val & 2 != 0 {
            val ^ 4
        } else {
            val ^ 21
        }
    } else {
        val
    };
    *c = new_val;
}

pub fn complement_alt(c: &mut u8) {
    let val = *c;
    if val != b'N' {
        if val & 2 != 0 {
            *c = val ^ 4
        } else {
            *c = val ^ 21
        };
    };
}

// This is the one used
pub fn revcomp_xor_pulp_alt_complement(sequence: &mut [u8]) {
    let arch = Arch::new();
    arch.dispatch(|| {
        sequence.reverse();
        sequence.make_ascii_uppercase();
        sequence.iter_mut().for_each(complement_alt);
    });
}

// Reverse complement through bit twiddling
pub fn revcomp_xor(sequence: &mut [u8]) {
    sequence.reverse();
    sequence.make_ascii_uppercase();
    for i in 0..sequence.len() {
        complement(&mut sequence[i]);
    }
}

// Reverse complement through bit twiddling
pub fn revcomp_xor_alt_complement(sequence: &mut [u8]) {
    sequence.reverse();
    sequence.make_ascii_uppercase();
    for i in 0..sequence.len() {
        complement_alt(&mut sequence[i]);
    }
}

// This is the one used
pub fn revcomp_xor_pulp(sequence: &mut [u8]) {
    let arch = Arch::new();
    arch.dispatch(|| {
        sequence.reverse();
        sequence.make_ascii_uppercase();
        sequence.iter_mut().for_each(complement);
    });
}

fn criterion_benchmark(c: &mut Criterion) {
    let sequence = b"ccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgataaaagcgtaccgtccccggcTCCTGTGAACTCCCTCCCCCCTTTTACCGCTTGATTGTACTAGAGTCGGTTTGTAGTCATCTTCCAGTTCTAAGCGTGTCACTCATGTGAAGTTAGTAGTTTCGTTTGTAGATCTGGCTAACTAAGGtgattaagttttatataattttttttaagttttgctaaaaatcattttagaagatattttttaaaaatttatgttcttttatgtggtcctttctcaaaatatattgtactgtatatttattataataaagtaccgtatttagttattaaaaatcagCTTGATCGTGTAATAAaaacacaggaaaaaaataaaatttatacaacaaattgtgaaatattaatattacaatgataaaaaataaagttatgaattaaaaatagccTAGGGCCTAGGCTATTAGAAATTGACTTACACATTTGATAACGTGACTCACTATAATGatagattttaatgtattataaaataatttggcaaCCAAAGTACGATTTTATCCATGTTTATGCATACCGCGTCTGACTATAAATCGTTGCATGTGTGGAGTACGTTTTCTTTGATCAGGTCTGTCAAACTTCATCGACAATTTAGCTTTAGACCGGTCTTCAGATAAATTGTTGCATTATGTGGGGGCGTAAGAGAACAGTAAGAACTCAATAATctaatagcctaatttttatttttggacatTTCAAAGCACTCTAAGAGAAAAGTGACACTCAATTAGCTGTAGGATATATGAGTACTGATATAGTTTTTTATGACCCTGttttagactttgaaattttctttcattgtcatTATTTACGGGATTTGTACAAACTGTGCTGGTGCAAGGCTATTTCAGTAACGAGTGTGTCTGATTTATGGCCCAGGAAATTTAAGACTATCTCAATTCTGAATTGAATGGGTTAACAAAGGCAGCAGTTATAAcctgaaataaacaattttaaaataacacaaatgctgtccaaaatattttttaattttataaaatgttatttagccTAATACCAGTAGTAGGTGACCCTAGTACAGCCATCAGTAGCCTATGATTCAGCAATATGATAAGATACAACAAATGTTTAATTGTATTCAAACTAAATATGTAAAaccttgaacattttttttgtcatatacaataaatactattAACAAAATACATTCCTCCGTGCCCCCTCCTCTCCCCTGGAATTAATAATTGctgattttgatttaaattaattttctgatttaagtcatatctttaattataattttgggtCAAATTATCTCACAAACCATgcataatcatattaaaaaaaatgcgctgtatttgtacttttaaaaaaaaatcactactcGTGAGAATCATGAGCAAAATATTCTAAAGTGGAAACGGCACTAAGGTGAACTAAGCAACTTAGTGCAAAActaaatgttagaaaaaatatCCTACACTGCATAAACTATTTTGcaccataaaaaaaagttatgtgtgGGTCTAAAATAATTTGCTGAGCAATTAATGATTTCTAAATGATGCTAAAGTGAACCATTGTAatgttatatgaaaaataaatacacaattaagATCAACACAGTGAAATAACATTGATTGGGTGATTTCAAATGGGGTCTATctgaataatgttttatttaacagtaatttttatttctatcaatttttagtaatatctacaaatattttgttttaggcTGCCAGAAGATCGGCGGTGCAAGGTCAGAGGTGAGATGTTAGGTGGTTCCACCAACTGCACGGAAGAGCTGCCCTCTGTCATTCAAAATTTGACAGGTACAAACAGactatattaaataagaaaaacaaactttttaaaggCTTGACCATTAGTGAATAGGTTATATGCTTATTATTTCCATTTAGCTTTTTGAGACTAGTATGATTAGACAAATCTGCTTAGttcattttcatataatattgaGGAACAAAATTTGTGAGATTTTGCTAAAATAACTTGCTTTGCTTGTTTATAGAGGCacagtaaatcttttttattattattataattttagattttttaatttttaaataagtgataGCATAtgctgtattattaaaaatttaagaactttaaagtatttacagtagcctatattatgtaaTAGGCTAGCCTACTTTATTGTTCGGccaattctttttcttattcatCTTATCATTATTagcagattattatttattactagtttaaaagcacgtcaaaaatgacggacattaatttcttcctccttaggttgctatacttaaatgccggtcaaaaatttaaaagcccgtcaaaaataacagacagcgacatctatggacagaaaatagagagaaacaTCTTCGGGCaacatcggctcgccgaagatggccgagcagtatgacagacatcatgaccaaactttctctttattatagttagattagatagaagattagaagactagtttaaaagcccatca";
    let mut sequence = sequence.to_vec();
    sequence.make_ascii_uppercase();

    let mut seq1 = sequence.clone();
    let mut seq2 = sequence.clone();
    let mut seq3 = sequence.clone();
    let mut seq4 = sequence.clone();

    revcomp(&mut seq1);
    revcomp_pulp(&mut seq2);
    revcomp_xor(&mut seq3);
    revcomp_xor_pulp(&mut seq4);

    assert!(seq1 == seq2);
    assert!(seq2 == seq3);
    assert!(seq3 == seq4);

    let mut group = c.benchmark_group("RevCmp");
    group.throughput(criterion::Throughput::Bytes(sequence.len() as u64));

    group.bench_function("revcomp_xor_pulp", |b| {
        b.iter(|| {
            revcomp_xor_pulp(black_box(&mut sequence));
        })
    });

    group.bench_function("revcomp_xor_pulp_alt_complement", |b| {
        b.iter(|| {
            revcomp_xor_pulp_alt_complement(black_box(&mut sequence));
        })
    });

    group.bench_function("revcomp_xor", |b| {
        b.iter(|| {
            revcomp_xor(black_box(&mut sequence));
        })
    });

    group.bench_function("revcomp_xor_alt_complement", |b| {
        b.iter(|| {
            revcomp_xor(black_box(&mut sequence));
        })
    });

    group.bench_function("revcomp", |b| {
        b.iter(|| {
            revcomp(black_box(&mut sequence));
        })
    });

    group.bench_function("revcomp_pulp", |b| {
        b.iter(|| {
            revcomp_pulp(black_box(&mut sequence));
        })
    });

    group.finish();
}

criterion_group! {
    name=parse_fasta;
    config = Criterion::default().significance_level(0.05).measurement_time(std::time::Duration::from_secs(10));
    targets=criterion_benchmark
}

criterion_main!(parse_fasta);
