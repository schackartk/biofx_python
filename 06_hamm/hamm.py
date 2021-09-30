#!/usr/bin/env python3
"""
Author : Kenneth Schackart <schackartk1@gmail.com>
Date   : 2021-09-29
Purpose: Calculate hamming distance
"""

import argparse
from itertools import zip_longest
from typing import NamedTuple, Tuple


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
def not_same(t: Tuple) -> bool:
    """ Compare tuple elements """

    return t[0] != t[1]


# --------------------------------------------------
def test_not_same() -> None:
    """ Test not_same()"""

    assert not_same(('A', 'A')) is False
    assert not_same(('A', 'T')) is True
    assert not_same(('A', None)) is True


# --------------------------------------------------
def get_dist(seq1: str, seq2: str) -> int:
    """ Calculate Hamming distance """

    return sum(map(not_same, zip_longest(seq1, seq2)))


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
