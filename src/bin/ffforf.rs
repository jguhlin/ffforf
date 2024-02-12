// TODO: It's fast enough right now, but needletail is better.
// TODO: Command line options
// TODO: Does this automatically handle reverse complement

use ffforf::*;

pub fn main() {
    // Get filename from arguments
    let args: Vec<String> = std::env::args().collect();
    let filename = &args[1];

    // Read file
    let file = std::fs::File::open(filename).unwrap();
    let mut reader = std::io::BufReader::new(file);
    let reader = fffx::fasta::Fasta::from_buffer(&mut reader);

    // Parse file
    for sequence in reader {
        let id = sequence.as_ref().unwrap().id.as_ref().unwrap().clone();
        let mut seqobj = sequence.unwrap();
        let sequence = seqobj.sequence.as_mut().unwrap();

        let results = find_stop_codons(sequence);
        let results = stop_codons_to_intervals(&results, 10, sequence.len());

        for result in results {
            let translated = translate_interval(sequence, &result);
            // Remove last entity (stop codon)
            let translated = &translated[0..translated.len()-1];
            if translated.len() < 50 {
                continue;                
            }
            let id = format!("{}_{}_{}_{}", id, result.0, result.1, result.2);
            println!(">{}", &id);
            let as_str: String = translated.iter().map(|x| Into::<char>::into(*x)).collect();
            println!("{}", as_str);
        }

        // Reverse complement

        revcomp(sequence);

        let results = find_stop_codons(sequence);
        let results = stop_codons_to_intervals(&results, 10, sequence.len());

        for mut result in results {
            let translated = translate_interval(sequence, &result);
            let translated = &translated[0..translated.len()-1];
            if translated.len() < 50 {
                continue;                
            }

            let id = format!("{}_rc_{}_{}_{}", id, result.0, result.1, result.2);
            println!(">{}", &id);
            let as_str: String = translated.iter().map(|x| Into::<char>::into(*x)).collect();
            println!("{}", as_str);
        }
    }
}
