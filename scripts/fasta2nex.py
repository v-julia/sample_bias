import argparse
import os
from Bio import AlignIO
from Bio import Alphabet

def fasta2nex(file_name):

    out_name = os.path.splitext(file_name)[0]+'.nex'
    AlignIO.convert(file_name, "fasta", out_name, "nexus",alphabet=Alphabet.generic_dna)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file in fasta format", required=True)
    args = parser.parse_args()
    fasta2nex(args.input_file)