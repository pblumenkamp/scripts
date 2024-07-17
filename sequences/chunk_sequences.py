import argparse
import gzip
import itertools
import sys

from typing import Callable, Dict, Iterator, List, TextIO


def main(file: str, chunkSize: int, slidingWindow: bool):
    f: TextIO
    get_sequence_lengths: Callable
    #Detect if FASTA or FASTQ
    with gzip.open(file, "rt") if file.endswith(".gz") else open(file) as f:
        for line in f:
            if not line.split() or line.startswith("#"):
                continue
            #FASTA file
            if line.startswith(">"):
                f.seek(0)
                chunk_fasta(f, chunkSize, slidingWindow)
                break
            #FASTQ file
            elif line.startswith("@"):
                next(f) # sequence
                line_with_plus: str = next(f)
                #Third line must start with plus else invalid FASTQ file
                if line_with_plus.startswith("+"):
                    f.seek(0)
                    chunk_fastq(f, chunkSize, slidingWindow)
                else:
                    print("Unknown file format: {}".format(file), file=sys.stderr)
                break
            else:
                print("Unknown file format: {}".format(file), file=sys.stderr)
                break


def chunk_fasta(file: TextIO, chunkSize: int, slidingWindow: bool) -> None:
    header: str = ""
    sequence: str = ""
    for line in file:
        if not line.strip() or line.startswith("#"):
            continue
        if line.startswith(">"):
            if header:
                if slidingWindow:
                    for i in range(len(sequence) - chunkSize + 1):
                        print("{} - {}".format(header, i))
                        print(sequence[i:i+chunkSize])
                else:
                    for i,seq in enumerate(itertools.batched(sequence, chunkSize)):
                        print("{} - {}".format(header, i))
                        print("".join(seq))
            header = line.strip()
            sequence = ""
        else:
            sequence += line.strip()
    
    if header:
        if slidingWindow:
            for i in range(len(sequence) - chunkSize + 1):
                print("{} - {}".format(header, i))
                print(sequence[i:i+chunkSize])
        else:
            for i,seq in enumerate(itertools.batched(sequence, chunkSize)):
                print("{} - {}".format(header, i))
                print("".join(seq))


def chunk_fastq(file: TextIO, chunkSize: int, slidingWindow: bool) -> None:
    for line in file:
        if not line.split() or line.startswith("#"):
            continue
        header: str = line.strip()
        sequence: str = next(file).strip()
        next(file) # + line
        quality: str = next(file) # quality line
        if slidingWindow:
            for i in range(len(sequence) - chunkSize + 1):
                print("{} - {}".format(header, i))
                print(sequence[i:i+chunkSize])
                print("+")
                print(quality[i:i+chunkSize])
        else:
            for i in range(0, len(sequence), chunkSize):
                print("{} - {}".format(header, int(i/chunkSize)))
                print(sequence[i:i+chunkSize])
                print("+")
                print(quality[i:i+chunkSize])


if __name__ == "__main__":
    if sys.version_info[0] != 3 or sys.version_info[1] < 12:
        print("This script requires Python version 3.12")
        sys.exit(1)
    parser = argparse.ArgumentParser(
                    prog=__file__,
                    description='Split FASTA/Q sequences into chunks')
    parser.add_argument('fasta', metavar="FASTA/Q")
    parser.add_argument('-c', '--chunk-size', dest="chunk_size", required=True, type=int)
    parser.add_argument('-s', '--sliding-window', dest="sliding_window", action='store_true')
    args = parser.parse_args()
    main(args.fasta, args.chunk_size, args.sliding_window)