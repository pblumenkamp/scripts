import gzip
import sys

from typing import Callable, Dict, TextIO


def main(kmer_length: int, file_path: str):
    count_kmers_at_start: Callable
    #Detect if FASTA or FASTQ
    with gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path) as file:
        for line in file:
            if not line.strip() or line.startswith("#"):
                continue
            #FASTA file
            if line.startswith(">"):
                count_kmers_at_start = kmers_at_start_in_fasta
                break
            #FASTQ file
            elif line.startswith("@"):
                next(file) # sequence
                line_with_plus: str = next(file)
                #Third line must start with plus else invalid FASTQ file
                if line_with_plus.startswith("+"):
                    count_kmers_at_start = kmers_at_start_in_fastq
                else:
                    print("Unknown file format: {}".format(file_path), file=sys.stderr)
                break
            else:
                print("Unknown file format: {}".format(file_path), file=sys.stderr)
                break
    statistics: Dict[str, int] = count_kmers_at_start(file_path, kmer_length)  # filename: {sequence: count} 
    
    print("{}\t{}".format("kmer", "count"))
    for kmer in sorted(statistics.items(), key=lambda x: x[1]):
        print("{}\t{}".format(kmer[0], kmer[1]))


# List all variations of the first k bases (kmer) in each sequence and count the number of each kmer
# Fasta variant
def kmers_at_start_in_fasta(file_path: str, kmer_length: int) -> Dict[str, int]:
    statistics: Dict[str, int] = {}
    with gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path) as fasta:
        sequence: str = ""
        for line in fasta:
            if not line.strip() or line.startswith("#"):
                continue
            if line.startswith(">"):
                if sequence:
                    kmer: str = sequence[0:kmer_length]
                    statistics[kmer] = statistics.get(kmer, 0) + 1
                    sequence = ""
            else:
                sequence += line.strip()
        # flush last sequence
        if sequence:
            kmer = sequence[0:kmer_length]
            statistics[kmer] = statistics.get(kmer, 0) + 1

    return(statistics)

# List all variations of the first k bases (kmer) in each sequence and count the number of each kmer
# Fastq variant
def kmers_at_start_in_fastq(file_path: str, kmer_length: int) -> Dict[str, int]:
    statistics: Dict[str, int] = {}
    with gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path) as fastq:
        # Ignore all lines but the second of four
        for line in fastq:
            if not line.strip() or line.startswith("#"):
                continue
            sequence = next(fastq).strip()
            statistics[sequence[0:kmer_length]] = statistics.get(sequence[0:kmer_length], 0) + 1
            next(fastq) # + line
            next(fastq) # quality line
    return(statistics)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("", file=sys.stderr)
        print("Summarize kmers at start of reads", file=sys.stderr)
        print("Usage: python3 {} <kmer length> <fasta/fastq>".format(__file__), file=sys.stderr)
        print("", file=sys.stderr)
        sys.exit(1)
    kmer_length: int = int(sys.argv[1])
    if kmer_length < 1:
        print("ERROR: kmer length must be >=1", file=sys.stderr)
        print("", file=sys.stderr)
        sys.exit(2)
    main(int(sys.argv[1]), sys.argv[2])