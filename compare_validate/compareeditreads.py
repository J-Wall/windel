#!/usr/bin/env python

import fire
import numpy as np
import pysam


def get_read_indel_dict(alignmentfile, contig, location):
    """
    Get the number of insertions and deletions for reads which overlap a specified
    location.

    Parameters
    ----------
    alignmentfile : pysam.AlignmentFile
    contig : str
    location : int
        0-based index within contig

    Returns
    -------
    dict of query_name : (n_insertions, n_deletions) pairs
    """
    d = {}
    for read in alignmentfile.fetch(contig=contig, start=location, stop=location + 1):
        d[read.query_name] = tuple(read.get_cigar_stats()[1][1:3])

    return d


def locus_venn(d0, d1):
    """
    Make a venn diagram of reads which align to a particular locus and its equivalent
    locus before and after correction (eg. by Windel)

    Parameters
    ----------
    d0 : dict
        dict of query_name : (n_insertions, n_deletions) pairs
    d1 : dict
        dict of query_name : (n_insertions, n_deletions) pairs

    Returns
    -------
    tuple of number of reads in d0 only, both d0 and d1, and d1 only
    """
    return (
        len(set(d0).difference(d1)),
        len(set(d0).intersection(d1)),
        len(set(d1).difference(d0)),
    )


def loci_venn(alignmentfile_before, alignmentfile_after, loci_before, loci_after):
    """
    Aggregates locus_venn over a list of loci

    Parameters
    ----------
    alignmentfile_before : pysam.AlignmentFile
    alignmentfile_after : pysam.AlignmentFile
    loci_before : iterable
        iterates over (str, int) tuples of (contig, location) pairs
    loci_after : iterable
        iterates over (str, int) tuples of (contig, location) pairs

    Returns
    -------
    3-tuple
        Venn diagram made from summing the venn diagrams at each loci
    """
    before_only, both, after_only = 0, 0, 0
    for locus_before, locus_after in zip(loci_before, loci_after):
        d_before = get_read_indel_dict(alignmentfile_before, *locus_before)
        d_after = get_read_indel_dict(alignmentfile_after, *locus_after)

        locus_before_only, locus_both, locus_after_only = locus_venn(d_before, d_after)

        before_only += locus_before_only
        both += locus_both
        after_only += locus_after_only

    return before_only, both, after_only


def parse_changes(changesfile):
    """
    Parses a changes file (output of generatewindelchanges.py)

    Parameters
    ----------
    changesfile : str
        Path to changes file (output of generatewindelchanges.py)

    Returns
    -------
    loci_before : list
        Contains (str, int) tuples of (contig, location) pairs
    loci_after : list
        Contains (str, int) tuples of (contig, location) pairs
    """
    with open(changesfile, "r") as f:
        lines = f.readlines()

    loci_before = []
    loci_after = []
    for line in lines:
        locus_before, locus_after = line.split()

        contig_before, location_before = locus_before.split(":")
        location_before = int(location_before.split("-")[0])
        loci_before.append((contig_before, location_before))

        contig_after, location_after = locus_after.split(":")
        location_after = int(location_after.split("-")[0])
        loci_after.append((contig_after, location_after))

    return loci_before, loci_after


def compare_edit_reads(before_bam, after_bam, changesfile):
    """
    Parameters
    ----------
    before_bam
        Path to bam file of reads mapped to before genome
    after_bam
        Path to bam file of reads mapped to after genome
    changesfile
        Path to changes file (output of generatewindelchanges.py)
    """
    alignmentfile_before = pysam.AlignmentFile(before_bam, "rb")
    alignmentfile_after = pysam.AlignmentFile(after_bam, "rb")
    loci_before, loci_after = parse_changes(changesfile)

    venn = loci_venn(alignmentfile_before, alignmentfile_after, loci_before, loci_after)
    print(f"Before:\t{venn[0]}\nBoth:\t{venn[1]}\nAfter:\t{venn[2]}")


def main():
    fire.Fire(compare_edit_reads)


if __name__ == "__main__":
    main()
