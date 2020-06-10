#!/usr/bin/env python

import fire
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


def iterate_indel_dicts(alignmentfile, loci):
    """
    Returns a generator which calls get_read_indel_dict for multiple loci

    Parameters
    ----------
    alignmentfile : pysam.AlignmentFile
    loci : iterator
        With 2-tuple elements of (contig, location) pairs

    Yields
    ------
    output of get_read_indel_dict
    """
    for locus in loci:
        yield get_read_indel_dict(alignmentfile, *locus)


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


def loci_venn(d_before_iter, d_after_iter):
    """
    Aggregates locus_venn over a list of loci

    Parameters
    ----------
    d_before_iter : iterable of dicts
        See iterate_indel_dicts
    d_after_iter : iterable of dicts
        See iterate_indel_dicts

    Returns
    -------
    3-tuple
        Venn diagram made from summing the venn diagrams at each loci
    """
    before_only, both, after_only = 0, 0, 0
    for d_before, d_after in zip(d_before_iter, d_after_iter):
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


def count_indels(indel_dict_iter):
    """
    Parameters
    ----------
    indel_dict_iter : iterator of dicts
        each element is a dict, with read_id : (n_insertions, n_deletions) key : values

    Returns
    -------
    tot_insertions : int
        total number of insertions
    tot_deletions : int
        total number of deletions
    """
    tot_insertions, tot_deletions = 0, 0
    for d in indel_dict_iter:
        try:
            n_insertions, n_deletions = zip(*d.values())
        except ValueError:
            continue
        tot_insertions += sum(n_insertions)
        tot_deletions += sum(n_deletions)

    return tot_insertions, tot_deletions


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

    d_before_list = list(iterate_indel_dicts(alignmentfile_before, loci_before))
    d_after_list = list(iterate_indel_dicts(alignmentfile_after, loci_after))

    venn = loci_venn(d_before_list, d_after_list)
    print("Number of reads mapping to edit loci")
    print("Note: Reads may be counted twice if edit-loci are close together")
    print(f"Before:\t{venn[0]}\nBoth:\t{venn[1]}\nAfter:\t{venn[2]}")
    print()

    tot_insertions_before, tot_deletions_before = count_indels(d_before_list)
    tot_insertions_after, tot_deletions_after = count_indels(d_after_list)
    print("Indel counts for edit-mapped reads\nBefore:")
    print(f"\tInsertions:\t{tot_insertions_before}")
    print(f"\tDeletions:\t{tot_deletions_before}")
    print("After:")
    print(f"\tInsertions:\t{tot_insertions_after}")
    print(f"\tDeletions:\t{tot_deletions_after}")


def main():
    fire.Fire(compare_edit_reads)


if __name__ == "__main__":
    main()
