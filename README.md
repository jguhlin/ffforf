# ffforf

`ffforf` is a Rust library and command-line tool for identifying open reading frames. It uses the jetscii crate for efficient searching of stop codons, needletail for fast FASTA parsing, and translates the ORFs into amino acid sequences. 

## Installation

To use the ORF Finder library in your Rust project, add it to your `Cargo.toml` file:

```toml
[dependencies]
ffforf = "0.3.0"
```

To install the `ffforf` binary, just `cargo install ffforf`

## Running the binary

To run 
```bash
ffforf genome.fna > translated_sequences.faa
```

Output looks like:
```
>Chr18_rc_2_3557_3872
ISTNLCTFLCSDTEFTPRVTNAKDSDTFDGILTLNNRQKHAERIAYNRGAGSGIGGGRGPGRPPITEIPLEELLACEEPEAKAARTRRRGATLALTALGRYIFN
```
Which is the landmark (contig/chromosome), whether reverse complement or not, reading frame, and start and end genomic coordinates.

### Note
- It does not fail gracefully, but if more people use it I will add in more command line arguments, help messages, etc...

- Min ORF size is 50, can be changed by altering the source code, for now. Please open an issue and I'll fix it right away.

- Unknown sequence gets read through, such as ```TTFLYLNYIITXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXVIGSYKEHFSVSTRDKPHVTKGRKERCNGNRITYYVIQNNFPALPTVSILYSLFQTQMIGRKNFA```

## Using the Library

To use the ORF Finder library in your Rust code, import the crate and use the `find_all_orfs` function:

```rust
use orf_finder::{find_all_orfs, Orf, Strand};

fn main() {
    let sequence = b"ATGCTAGTAACTAGCGTAA";
    let min_orf_length = 5;
    let orfs = find_all_orfs(sequence, min_orf_length);

    for orf in orfs {
        println!(
            "ORF: Start: {}, End: {}, Strand: {}, Reading Frame: {}",
            orf.start, orf.end, orf.strand, orf.reading_frame
        );
    }
}
```

You can also check out src/bin/ffforf.rs for more potential usage.

## License
This project is licensed under the MIT License. A copy can typically be found easily...