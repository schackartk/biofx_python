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
def get_dist(seq1: str, seq2: str) -> int:
    """ Calculate Hamming distance """

    dist = abs(len(seq1) - len(seq2))

    dist += sum(map(lambda x: x[0] != x[1], zip(seq1, seq2)))

    return dist


# --------------------------------------------------
def test_get_dist() -> None:
    """ test get_dist() """

    assert get_dist('ATGC', 'ATGC') == 0
    assert get_dist('ATGC', 'ATGG') == 1
    assert get_dist('ATGC', 'ATG') == 1
    assert get_dist('ATGC', 'AGG') == 2


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()

    print(get_dist(args.seq1, args.seq2))


# --------------------------------------------------
if __name__ == '__main__':
    main()
