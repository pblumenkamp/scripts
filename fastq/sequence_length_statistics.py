import gzip
import sys

from typing import Callable, Dict, List, TextIO


def main(files: List[str]):
    statistics: Dict[str, Dict[int, int]] = {}  # filename: {length: count}
    for file in files:
        f: TextIO
        get_sequence_lengths: Callable
        #Detect if FASTA or FASTQ
        with gzip.open(file, "rt") if file.endswith(".gz") else open(file) as f:
            for line in f:
                if not line.split() or line.startswith("#"):
                    continue
                #FASTA file
                if line.startswith(">"):
                    get_sequence_lengths = sequence_lengths_in_fasta
                    break
                #FASTQ file
                elif line.startswith("@"):
                    next(f) # sequence
                    line_with_plus: str = next(f)
                    #Third line must start with plus else invalid FASTQ file
                    if line_with_plus.startswith("+"):
                        get_sequence_lengths = sequence_lengths_in_fastq
                    else:
                        print("Unknown file format: {}".format(file), file=sys.stderr)
                    break
                else:
                    print("Unknown file format: {}".format(file), file=sys.stderr)
                    break
        statistics[file] = get_sequence_lengths(file)
    
    # Find length of longest read between all files (calculate length of TSV)
    length_of_longest_read: int = 0
    for file_stats in statistics.values():
        length_of_longest_read = max(length_of_longest_read, max(*file_stats.keys()))

    print("sequence_length\t{}".format("\t".join(files)))
    for read_length in range(length_of_longest_read):
        counts: List[int] = [str(statistics[file].get(read_length+1, 0)) for file in files]
        print("{}\t{}".format(read_length+1, "\t".join(counts)))


def sequence_lengths_in_fasta(file: str) -> Dict[int, int]:
    statistics: Dict[int, int] = {}
    with gzip.open(file, "rt") if file.endswith(".gz") else open(file) as fasta:
        sequence: str = ""
        for line in fasta:
            if not line.split() or line.startswith("#"):
                continue
            if line.startswith(">"):
                statistics[len(sequence)] = statistics.get(len(sequence), 0) + 1
                sequence = ""
            else:
                sequence += line.strip()
    return(statistics)


def sequence_lengths_in_fastq(file: str) -> Dict[int, int]:
    statistics: Dict[int, int] = {}
    with gzip.open(file, "rt") if file.endswith(".gz") else open(file) as fastq:
        for line in fastq:
            if not line.split() or line.startswith("#"):
                continue
            sequence = next(fastq).strip()
            statistics[len(sequence)] = statistics.get(len(sequence), 0) + 1
            next(fastq) # + line
            next(fastq) # quality line
    return(statistics)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("", file=sys.stderr)
        print("Summarize read lengths as TSV", file=sys.stderr)
        print("Usage: python3 {} <fasta/fastq>...".format(__file__), file=sys.stderr)
        print("", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1:])