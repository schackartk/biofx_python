#!/usr/bin/env python3
"""
Author : Kenneth Schackart <schackartk1@gmail.com>
Date   : 2021-07-21
Purpose: Compute GC content
"""

import argparse
import re
import sys
from typing import NamedTuple, TextIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
def get_gc(rec: SeqRecord) -> MySeq:
    """ Calculate GC content"""

    pct = 0.

    if seq := str(rec.seq):
        gc = len(re.findall('[GC]', seq.upper()))
        pct = (gc * 100) / len(seq)

    return MySeq(pct, rec.id)


# --------------------------------------------------
def test_get_gc() -> None:
    """ Test get_gc()"""

    assert get_gc(SeqRecord(Seq(''), id='123')) == (0.0, '123')
    assert get_gc(SeqRecord(Seq('C'), id='ABC')) == (100.0, 'ABC')
    assert get_gc(SeqRecord(Seq('G'), id='XYZ')) == (100.0, 'XYZ')
    assert get_gc(SeqRecord(Seq('ACTG'), id='JKL')) == (50.0, 'JKL')
    assert get_gc(SeqRecord(Seq('GGCC'), id='DEF')) == (100.0, 'DEF')


# --------------------------------------------------
def main() -> None:
    """ Do the stuff """

    args = get_args()

    high = MySeq(0., '')

    for seq in map(get_gc, SeqIO.parse(args.file, 'fasta')):
        if seq.gc > high.gc:
            high = seq

    print(f'{high.name} {high.gc:.6f}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
