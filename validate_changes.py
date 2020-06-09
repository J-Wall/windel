#!/usr/bin/env python

import sys

from pyfaidx import Fasta

before_genome = Fasta(sys.argv[1])
after_genome = Fasta(sys.argv[2])
with open(sys.argv[3]) as f:
    changes = [l.split() for l in f.readlines()]

checksize = int(sys.argv[4])

prev_line = None

for before_loc, after_loc in changes:
    before_chr, before_idx = before_loc.split(":")
    after_chr, after_idx = after_loc.split(":")
    before_idx = int(before_idx.split("-")[0])
    after_idx = int(after_idx.split("-")[0])

    before_seq = before_genome[before_chr][before_idx - checksize : before_idx].seq
    after_seq = after_genome[after_chr][after_idx - checksize : after_idx].seq
    line = f"{before_loc} {after_loc} {before_seq} {after_seq}"

    if before_seq != after_seq:
        print(prev_line)
        print(line)
        break

    prev_line = line

