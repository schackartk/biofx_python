#!/usr/bin/env python3
"""
Author : Kenneth Schackart <schackartk1@gmail.com>
Date   : 2021-05-13
Purpose: Print the reverse complement of DNA
"""

import argparse
from Bio.Seq import Seq
import os
from typing import NamedTuple, TextIO


class Args(NamedTuple):
    """ Command-line arguments """
    dna: str


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Print the reverse complement of DNA',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dna',
                        metavar='DNA',
                        help='Input sequence or file')


    args = parser.parse_args()

    if os.path.isfile(args.dna):
        args.dna = open(args.dna).read().rstrip()

    return Args(dna=args.dna)


# --------------------------------------------------
def main() -> None:
    """ Do the stuff """

    args = get_args()
    dna = Seq(args.dna)

    revc = dna.reverse_complement()

    print(revc)

# --------------------------------------------------
if __name__ == '__main__':
    main()
