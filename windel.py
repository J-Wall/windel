import collections
import csv
import re
import time
from multiprocessing import Pool

import fire
import numpy as np
import pysam
from scipy import stats


class WindowedIndels:
    """
    Iterates over a pysam.AlignmentFile.pileup and accumulates indels over a
    sliding window of odd width. Optionally supresses indels at positions which
    have a non-maximal number of indels within the surrounding window (default)

    Parameters
    ----------
    pileup : iterator
        Output of pysam.AlignmentFile.pileup
    width : int
        Width of window. Must be odd.
    supress_nonmax : bool
        If True (default), supress non-maximal indels.

    Yields
    ------
    ref_name : str
        Name of reference scaffold/contig.
    position : int
        Position in assembly.
    coverage : int
        Coverage at position.
    w_coverage : int
        Windowed coverage at position.
    best_nucs : str
        Most common insertion.
    best_count : int
        Count of best_nucs.
    indels : dict
        Dictionary with fields `queryname` : `insertion_length` (-ve for deletion).
    w_indels : dict
        Window-accumulated version of indels.
    """

    def __init__(self, pileup, width, supress_nonmax=True):
        """
        Parameters
        ----------
        pileup : iterator
            Output of pysam.AlignmentFile.pileup
        width : int
            Width of window. Must be odd.
        supress_nonmax : bool
            If True (default), supress non-maximal indels.
        """
        assert width % 2 == 1, "width must be odd."

        self.pileup = pileup
        self.width = width
        self.supress_nonmax = supress_nonmax

        self._temp_column = None
        self._temp_pos = 0
        self._pos = np.inf
        self._ref_name = ""
        self._ringbuf = [self._next_pos() for i in range(width - 1)]
        self._ringbuf.append(None)

        self._i = width - 1
        self._centre = width // 2 - 1

    def __iter__(self):
        return self

    def __next__(self):
        """
        Returns
        -------
        ref_name : str
            Name of reference scaffold/contig.
        position : int
            Position in assembly.
        coverage : int
            Coverage at position.
        w_coverage : int
            Windowed coverage at position.
        best_nucs : str
            Most common insertion.
        best_count : int
            Count of best_nucs.
        indels : dict
            Dictionary with fields `queryname` : `insertion_length` (-ve for deletion).
        w_indels : dict
            Window-accumulated version of indels.
        """
        self._ringbuf[self._i] = self._next_pos()
        self._i = (self._i + 1) % self.width
        self._centre = (self._centre + 1) % self.width

        if self.supress_nonmax and len(self._ringbuf[self._centre][6]) != max(
            [len(entry[6]) for entry in self._ringbuf]
        ):
            return self._ringbuf[self._centre] + ({},)
        else:
            counter = collections.Counter()
            for entry in self._ringbuf:
                counter.update(entry[6])
            return self._ringbuf[self._centre] + (dict(counter),)

    def _next_pos(self):
        """
        Returns
        -------
        ref_name : str
            Name of reference scaffold/contig.
        position : int
            Position in assembly.
        coverage : int
            Coverage at position.
        w_coverage : int
            Windowed coverage at position.
        best_nucs : str
            Most common insertion.
        best_count : int
            Count of best_nucs.
        indels : dict
            Dictionary with fields `queryname` : `insertion_length` (-ve for deletion).
        """
        if self._pos < self._temp_pos - 1:
            self._pos += 1
            return self._ref_name, self._pos, 0, 0, "", 0, {}

        self._temp_column = next(self.pileup)
        self._temp_pos = self._temp_column.pos
        self._ref_name = self._temp_column.reference_name

        if self._pos < self._temp_pos - 1:
            self._pos += 1
            return self._ref_name, self._pos, 0, 0, "", 0, {}

        w_coverage = 0
        indels = {}
        column = self._temp_column
        for read in column.pileups:
            if (
                read.alignment.reference_start > column.reference_pos - self.width // 2
                or read.alignment.reference_end < column.reference_pos + self.width // 2
            ):
                continue
            else:
                w_coverage += 1
                if read.indel != 0:
                    indels[read.alignment.query_name] = read.indel

        query_sequences = column.get_query_sequences(add_indels=True)
        indel_seqs = [
            re.sub("[0-9]", "", q[1:].upper()) for q in query_sequences if len(q) > 1
        ]

        mode, mode_count = stats.mode(indel_seqs)
        mode = mode[0] if len(mode) > 0 else ""
        mode_count = mode_count[0] if len(mode_count) > 0 else ""

        self._pos = column.reference_pos
        return (
            column.reference_name,
            column.pos,
            column.nsegments,
            w_coverage,
            mode,
            mode_count,
            indels,
        )


def summarise_column(
    ref_name, position, coverage, w_coverage, best_nucs, best_count, indels, w_indels
):
    """
    Takes the ouput of one iteration of WindowedIndels and gives a summary statistics
    relating to the genomic position.

    Parameters
    ----------
    ref_name : str
        Name of reference scaffold/contig.
    position : int
        Position in assembly.
    coverage : int
        Coverage at position.
    w_coverage : int
        Windowed coverage at position.
    best_nucs : str
        Most common insertion.
    best_count : int
        Count of best_nucs.
    indels : dict
        Dictionary with fields `queryname` : `insertion_length` (-ve for deletion).
    w_indels : dict
        Window-accumulated version of indels.

    Returns
    -------
    summary : dict
        Summary statistics
    """
    length_mode, length_mode_count = stats.mode(list(indels.values()))
    length_mode, length_mode_count = (
        (int(length_mode[0]), int(length_mode_count[0]))
        if len(length_mode > 0)
        else (0, 0)
    )

    w_length_mode, w_length_mode_count = stats.mode(list(w_indels.values()))
    w_length_mode, w_length_mode_count = (
        (int(w_length_mode[0]), int(w_length_mode_count[0]))
        if len(w_length_mode > 0)
        else (0, 0)
    )
    insert_seq = best_nucs.strip("+-")

    summary = dict(
        ref_name=ref_name,
        position=position,
        prop_indel=len(indels) / w_coverage if w_coverage > 0 else 0,
        prop_mode=float(length_mode_count / w_coverage) if w_coverage > 0 else 0,
        length_mode=length_mode,
        prop_w_indel=len(w_indels) / w_coverage if w_coverage > 0 else 0,
        prop_w_mode=float(w_length_mode_count / w_coverage) if w_coverage > 0 else 0,
        w_length_mode=w_length_mode,
        insertion=True if best_nucs.startswith("+") else False,
        deletion=True if length_mode < 0 else False,
        indel_length=len(insert_seq),
        insert_seq=insert_seq.strip("N"),
        lengths_match=bool(length_mode == len(insert_seq)),
    )

    return summary


def summarise(windowed_indels):
    """
    Iterate over windowed_indels and call summarise_column on each column.

    Parameters
    ----------
    windowed_indels : WindowedIndels

    Returns
    -------
    Generative iterator of summary dicts.
    """
    for column in windowed_indels:
        summary = summarise_column(*column)
        yield summary


def get_summary_arrays(
    windowed_indels,
    labels=[
        "ref_name",
        "position",
        "prop_indel",
        "prop_mode",
        "length_mode",
        "prop_w_indel",
        "prop_w_mode",
        "w_length_mode",
        "insertion",
        "deletion",
        "indel_length",
        "insert_seq",
        "lengths_match",
    ],
):
    """
    Call summarise(windowed_indels) and output the desired fields to arrays

    Parameters
    ----------
    windowed_indels : WindowedIndels
    labels : list
        Names of arrays to return.
    """
    summary_lists = [[] for i in range(len(labels))]
    for column_summary in summarise(windowed_indels):
        for i, label in enumerate(labels):
            try:
                summary_lists[i].append(column_summary[label])
            except KeyError:
                q = "'"
                raise ValueError(
                    f"""{label} is not in summary dict. Possible values are:
                        {', '.join([q + k + q for k in column_summary.keys()])}"""
                )
    return tuple(np.array(l) for l in summary_lists)


def threshold(alignmentfile, region=None, prop_thresh=0.0, window=11):
    """
    Performs windowed indel analysis.

    Parameters
    ----------
    alignmentfile : str
        path to alignment file (bam/sam)
    region : str
        regions can be specified as: RNAME[:STARTPOS[-ENDPOS]] (same as samtools)
    prop_thresh : float
        threshold for the proportion of windowed reads having an indel at a position for
        an edit to be made (default 0.9)
    window : int
        width of window to use, must be odd (default 11)

    Returns
    -------
    summaries : list
        A list of indel summary dicts.
    """
    print(
        f"{time.strftime('%Y-%m-%d %H:%M:%S')}:",
        f"processing {alignmentfile}, region={region}",
    )
    t0 = time.perf_counter()
    samfile = pysam.AlignmentFile(alignmentfile, "rb")
    pileup = samfile.pileup(region=region)
    windowed_indels = WindowedIndels(pileup, window)

    summaries = []
    for summary in summarise(windowed_indels):
        if summary["prop_w_indel"] > prop_thresh:
            summaries.append(summary)

    samfile.close()

    t1 = time.perf_counter()
    t = t1 - t0
    m, s = t // 60, t % 60
    h, m = m // 60, m % 60
    print(
        f"{time.strftime('%Y-%m-%d %H:%M:%S')}:",
        f"finished processing {region}, {len(summaries)} possible indel sites found.",
        f"[{h:02.0f}:{m:02.0f}:{s:05.2f}]",
    )
    return summaries


def _threshold_1arg(args_dict):
    return threshold(**args_dict)


def process_alignments(
    referencefile,
    alignmentfile,
    outfile,
    fasta_out=None,
    prop_thresh=0.0,
    window=11,
    processes=None,
):
    """
    Performs windowed indel analysis in a multiprocessing fashion.

    Parameters
    ----------
    referencefile : str
        path to reference fasta file
    alignmentfile : str
        path to alignment file (bam/sam)
    outfile : str
        path to output file (csv)
    fasta_out : str
        Optional path to corrected fasta output (default None).
    prop_thresh : float
        threshold for the proportion of windowed reads having an indel at a position for
        an edit to be made (default 0.9)
    window : int
        width of window to use, must be odd (default 11)
    processes : int
        number of processes to spawn (default is the number of cpus)
    """
    pool = Pool(processes=processes)

    fasta = pysam.FastaFile(referencefile)
    if fasta_out is not None:
        fasta_out_file = open(fasta_out, "w")
    with open(outfile, "w") as f:
        writer = csv.DictWriter(
            f,
            [
                "ref_name",
                "position",
                "prop_indel",
                "prop_mode",
                "length_mode",
                "prop_w_indel",
                "prop_w_mode",
                "w_length_mode",
                "insertion",
                "deletion",
                "indel_length",
                "insert_seq",
                "lengths_match",
            ],
            dialect="excel-tab",
        )
        writer.writeheader()

        for summaries, seqname in zip(
            pool.imap(
                _threshold_1arg,
                (
                    {
                        "alignmentfile": alignmentfile,
                        "region": region,
                        "prop_thresh": prop_thresh,
                        "window": window,
                    }
                    for region in fasta.references
                ),
            ),
            fasta.references,
        ):
            writer.writerows(summaries)
            if fasta_out is not None:
                positions = np.array([s["position"] for s in summaries], dtype="i8")
                insert_nucs = np.array([s["insert_seq"] for s in summaries], dtype=str)
                insertion_mask = np.array(
                    [s["insertion"] for s in summaries], dtype=bool
                )
                deletion_mask = np.array([s["deletion"] for s in summaries], dtype=bool)
                deletion_lengths = np.array(
                    [s["indel_length"] for s in summaries], dtype="i8"
                )[deletion_mask]
                masked_seq = mask_deletions(
                    fasta.fetch(seqname), positions[deletion_mask], deletion_lengths
                )
                corrected_seq = correct_insertions(
                    masked_seq, positions[insertion_mask], insert_nucs[insertion_mask]
                )
                corrected_seq = corrected_seq.replace("-", "")
                print(f">{seqname}", file=fasta_out_file)
                print(
                    "\n".join(
                        (
                            corrected_seq[i : i + 80]
                            for i in range(0, len(corrected_seq), 80)
                        )
                    ),
                    file=fasta_out_file,
                )

    fasta.close()
    if fasta_out is not None:
        fasta_out_file.close()
    pool.close()
    pool.join()


def correct_insertions(sequence, positions, insert_nucs):
    """
    Inserts nucleotides into sequence after the specified 0-based positions

    Parameters
    ----------
    sequence : str
        Nucleotide sequence to be modified.
    positions : array_like
        Positions in sequence after which insertions are to be made (0-based).
    insert_nucs : array_like
        Array of insertion strings

    Returns
    -------
    corrected_sequence : str
        Corrected sequence.
    """
    insert_array = np.array(insert_nucs)
    s = np.array([c for c in sequence], dtype=insert_array.dtype)
    return "".join(np.insert(s, np.array(positions) + 1, insert_array))


def mask_deletions(sequence, positions, deletion_lengths):
    """
    Replaces regions in `sequence` specified by `positions` and `deletion_lengths` with
    '-'s. Should call this before `correct_insertions`.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence to be modified.
    positions : array_like
        Positions in sequence after which deletions are to be made (0-based).
    deletion_lengths : array_like
        Lengths of deletions. Should have the same number of elements as `positions`.

    Returns
    -------
    corrected_sequence : str
        Corrected sequence
    """
    s = np.array([c for c in sequence], dtype=str)
    for position, length in zip(positions, deletion_lengths):
        s[position + 1 : position + length + 1] = "-"

    return "".join(s)


def main():
    fire.Fire(process_alignments)


if __name__ == "__main__":
    main()
