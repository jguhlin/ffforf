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
        let mut seqobj = sequence.unwrap();
        let sequence = seqobj.sequence.as_mut().unwrap();
        let results = find_stop_codons(sequence);
        let results = stop_codons_to_intervals(&results, 10, sequence.len());
        for result in results {
            let translated = translate_interval(sequence, &result);
            todo!();
        }
    }
}
