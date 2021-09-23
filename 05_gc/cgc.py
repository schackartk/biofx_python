#!/usr/bin/env python3
"""
Author : Kenneth Schackart <schackartk1@gmail.com>
Date   : 2021-07-21
Purpose: Compute GC content
"""

import argparse
import sys
from typing import NamedTuple, TextIO
from Bio import SeqIO


class Args(NamedTuple):
    """ Command-line arguments """
    file: TextIO


# --------------------------------------------------
class MySeq(NamedTuple):
    """ Sequence """
    gc: float
    name: str


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Compute GC content',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file',
                        help='Input sequence file',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        default=sys.stdin,
                        nargs='?')

    args = parser.parse_args()

    return Args(args.file)


# --------------------------------------------------
def get_gc(seq: str) -> float:
    """ Calculate GC content"""

    return (seq.count('G') +
    seq.count('C')) * 100 / len(seq) if seq else 0


# --------------------------------------------------
def test_get_gc() -> None:
    """ Test get_gc()"""

    assert get_gc('C') == 100.
    assert get_gc('G') == 100.
    assert get_gc('CGCGCG') == 100.
    assert get_gc('ATATAT') == 0.
    assert get_gc('ATGC') == 50.


# --------------------------------------------------
def main() -> None:
    """ Do the stuff """

    args = get_args()

    high = MySeq(0., '')

    for rec in SeqIO.parse(args.file, 'fasta'):
        pct = get_gc(rec.seq)

        if pct > high.gc:
            high = MySeq(pct, rec.id)

    print(f'{high.name} {high.gc:.6f}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
