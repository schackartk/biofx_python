#!/usr/bin/env python3
"""
Author : Kenneth Schackart <schackartk1@gmail.com>
Date   : 2021-09-29
Purpose: Calculate hamming distance
"""

import argparse
from typing import NamedTuple


class Args(NamedTuple):
    """ Command-line arguments """
    seq1: str
    seq2: str


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Calculate Hamming distance',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('seq1', metavar='str', help='Sequence 1')

    parser.add_argument('seq2', metavar='str', help='Sequence 2')

    args = parser.parse_args()

    return Args(args.seq1, args.seq2)


# --------------------------------------------------
def get_dist(args: Args) -> int:
    """ Calculate Hamming distance """

    seq1 = args.seq1
    seq2 = args.seq2
    dist = abs(len(seq1) - len(seq2))

    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] != seq2[i]:
            dist += 1

    return dist


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()

    print(get_dist(args))


# --------------------------------------------------
if __name__ == '__main__':
    main()
