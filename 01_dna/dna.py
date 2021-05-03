#!/usr/bin/env python3
"""
Author : ken <ken@localhost>
Date   : 2021-05-03
Purpose: Tetranucleotide frequency
"""

import argparse
from typing import NamedTuple, TextIO, Tuple
import os


class Args(NamedTuple):
    """ Command-line arguments """
    dna: str

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Tetranucleotide frequency',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dna',
                        metavar='DNA',
                        help='Input DNA sequence')

    args = parser.parse_args()

    if os.path.isfile(args.dna):
        args.dna = open(args.dna).read().rstrip()

    return Args(args.dna)

# --------------------------------------------------
def count_bases(dna: str) -> Tuple[int, int, int, int]:
    """ Count each of the bases in DNA """

    cnt_a, cnt_c, cnt_g, cnt_t = 0, 0, 0, 0
    for base in dna:
        if base == 'A':
            cnt_a += 1
        elif base == 'C':
            cnt_c += 1
        elif base == 'G':
            cnt_g += 1
        elif base == 'T':
            cnt_t += 1

    return (cnt_a, cnt_c, cnt_g, cnt_t)

# --------------------------------------------------
def test_count_bases() -> None:
    """ Test count """

    assert count_bases('') == (0, 0, 0, 0)
    assert count_bases('123XYZ') == (0, 0, 0, 0)
    assert count_bases('A') == (1, 0, 0, 0)
    assert count_bases('C') == (0, 1, 0, 0)
    assert count_bases('G') == (0, 0, 1, 0)
    assert count_bases('T') == (0, 0, 0, 1)
    assert count_bases('AAAACCCGGT') == (4, 3, 2, 1)

# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()  
    print('{} {} {} {}'.format(*count_bases(args.dna)))

# --------------------------------------------------
if __name__ == '__main__':
    main()
