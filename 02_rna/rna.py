#!/usr/bin/env python3
"""
Author : Kenneth Schackart <schackartk1@gmail.com>
Date   : 2021-05-13
Purpose: Translate DNA to RNA
"""

import argparse
import os
from typing import NamedTuple, List, TextIO


class Args(NamedTuple):
    """ Command-line arguments """
    files: List[TextIO]
    out_dir: str

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Convert DNA string into RNA',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file',
                        metavar='FILE',
                        help='Input DNA file(s)',
                        type=argparse.FileType('rt'),
                        nargs='+')

    parser.add_argument('-o',
                        '--out_dir',
                        help='Output directory',
                        metavar='DIR',
                        type=str,
                        default='out')

    args = parser.parse_args()

    return Args(files=args.file, out_dir=args.out_dir)


# --------------------------------------------------
def main() -> None:
    """ Do the stuff """

    args = get_args()
    files = args.files
    out_dir = args.out_dir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    num_files, num_seqs = 0, 0
    for fh in files:
        num_files += 1
        out_file = os.path.join(out_dir, os.path.basename(fh.name))
        out_fh = open(out_file, 'wt')
        
        for dna in fh:
            num_seqs += 1
            print(dna.replace('T', 'U'), file = out_fh, end = '')

        out_fh.close()
        

    print(f'Done, wrote {num_seqs} sequence{"" if num_seqs == 1 else "s"}',
    f'in {num_files} file{"" if num_files == 1 else "s"}',
    f'to directory "{out_dir}".')



# --------------------------------------------------
if __name__ == '__main__':
    main()
