"""
Script for converting WIG files with fixed step entries to a WIG file with variable steps.
https://m.ensembl.org/info/website/upload/wig.html
"""

from types import SimpleNamespace
from typing import Dict
import sys
import argparse


def main():
    args = parse_arguments()
    out_file = sys.stdout if args.output is None else open(args.output, "w")
    with open(args.wigfile) as in_file:
        #possible states: fixedStep, fixedStepSpan, variableStep
        STATE = SimpleNamespace(start="start", fixedStep="fixedStep", fixedStepSpan="fixedStepSpan", variableStep="variableStep")
        state: STATE = STATE.start
        current_wig_parameters: Dict[str, str] = None
        current_pos: int = None
        for line in in_file:
            # Detect if next block is fixedStep (with or without span) or variableStep,
            # save all parameters of this block (chromosome, start, step, span),
            # and print new variableStep header
            if line.startswith("fixedStep"):
                current_wig_parameters = {kv.split("=")[0]: kv.split("=")[1] for kv in line.strip().split(" ")[1:]}
                current_pos = int(current_wig_parameters["start"])
                if "span" in current_wig_parameters:
                    state = STATE.fixedStepSpan
                    print("variableStep chrom={} span={}".format(current_wig_parameters["chrom"], current_wig_parameters["span"]), file=out_file)
                else:
                    state = STATE.fixedStep
                    print("variableStep chrom={}".format(current_wig_parameters["chrom"]), file=out_file)
            elif line.startswith("variableStep"):
                state = STATE.variableStep
                print(line.strip(), file=out_file)
            else:
                # start data block

                # ignore all line before first declaration line (could be extended to support track definition lines)
                if state == STATE.start:
                    continue
                # print data values
                elif state == STATE.fixedStepSpan:
                    print("{} {}".format(current_pos, line.strip()), file=out_file)
                    current_pos += int(current_wig_parameters["step"])
                elif state == STATE.fixedStep:
                    print("{} {}".format(current_pos, line.strip()), file=out_file)
                    current_pos += int(current_wig_parameters["step"])
                elif state == STATE.variableStep:
                    print(line.strip(), file=out_file)
    
    if args.output is not None:
        out_file.close()


def parse_arguments():
    arg = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    arg.add_argument("wigfile", help="Input WIG file")
    arg.add_argument("-o", "--output", type=str, required=False, dest="output",
                   help="Write output to file [Default: STDOUT]")

    return(arg.parse_args())


if __name__ == '__main__':
    
    if sys.version_info<(3,5,0):
        sys.stderr.write("You need python 3.5 or later to run this script\n")
        sys.exit(1)
    main()
