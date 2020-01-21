Windel
------

Sliding-window INDEL correction

```
NAME
    windel.py - Performs windowed indel analysis in a multiprocessing fashion.

SYNOPSIS
    windel.py REFERENCEFILE ALIGNMENTFILE OUTFILE <flags>

DESCRIPTION
    Performs windowed indel analysis in a multiprocessing fashion.

POSITIONAL ARGUMENTS
    REFERENCEFILE
        path to reference fasta file
    ALIGNMENTFILE
        path to alignment file (bam/sam)
    OUTFILE
        path to output file (csv)

FLAGS
    --fasta_out=FASTA_OUT
        Optional path to corrected fasta output (default None).
    --prop_thresh=PROP_THRESH
        threshold for the proportion of windowed reads having an indel at a position for an edit to be made (default 0.9)
    --window=WINDOW
        width of window to use, must be odd (default 11)
    --processes=PROCESSES
        number of processes to spawn (default is the number of cpus)

NOTES
    You can also use flags syntax for POSITIONAL ARGUMENTS
```
