// TODO: Command line options
// TODO: Does this automatically handle reverse complement

use ffforf::*;
use needletail::{parse_fastx_file, Sequence};

pub fn main() {
    // Get filename from arguments
    let args: Vec<String> = std::env::args().collect();
    let filename = &args[1];

    // Read file
    let mut reader =
        parse_fastx_file(&filename).expect("Invalid fasta/x filename or file not found");

    // Parse file
    while let Some(record) = reader.next() {
        let record = record.expect("invalid record");
        let id = record.id();
        let sequence = record.seq();

        let results = find_stop_codons(&sequence);
        let results = stop_codons_to_intervals(&results, 10, sequence.len());

        for result in results {
            let translated = translate_interval(&sequence, &result);
            // Remove last entity (stop codon)
            let translated = &translated[0..translated.len() - 1];
            if translated.len() < 50 {
                continue;
            }
            let id = format!(
                "{}_{}_{}_{}",
                unsafe { std::str::from_utf8_unchecked(&id) },
                result.0,
                result.1,
                result.2
            );
            println!(">{}", &id);
            let as_str: String = translated.iter().map(|x| Into::<char>::into(*x)).collect();
            println!("{}", as_str);
        }

        // Reverse complement

        let rc = sequence.reverse_complement();

        let results = find_stop_codons(&rc);
        let results = stop_codons_to_intervals(&results, 10, rc.len());

        for result in results {
            let translated = translate_interval(&rc, &result);
            let translated = &translated[0..translated.len() - 1];
            if translated.len() < 50 {
                continue;
            }

            let id = format!(
                "{}_rc_{}_{}_{}",
                unsafe { std::str::from_utf8_unchecked(&id) },
                result.0,
                result.1,
                result.2
            );
            println!(">{}", &id);
            let as_str: String = translated.iter().map(|x| Into::<char>::into(*x)).collect();
            println!("{}", as_str);
        }
    }
}
