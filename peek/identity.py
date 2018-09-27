import contextlib
import re
import sys

import click
import pysam
from tqdm import tqdm


def get_discont_pos(line, min_gap=6):
    '''Yields splice site and position on the reference sequence.

    Yields pairs of (x, y) where

    x is the first base of the poly(N) stretch
    y is first base after the poly(N) stretch

    http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/
    Operator    Description
    D   Deletion; the nucleotide is present in the reference but not in the read
    H   Hard Clipping; the clipped nucleotides are not present in the read.
    I   Insertion; the nucleotide is present in the read but not in the reference.
    M   Match; can be either an alignment match or mismatch. 
    N   Skipped region; a region of nucleotides is not present in the read
    P   Padding; padded area in the read and not in the reference
    S   Soft Clipping;  the clipped nucleotides are present in the read
    X   Read Mismatch; the nucleotide is present in the reference
    =   Read Match; the nucleotide is present in the reference
    '''
    import re

    c = re.compile('([0-9]+)([SIDMN]+)')

    pos = line.reference_start
    # position on reference in the alignment 
    # the (sam, bam) coordinate system starts w/ (0, 1)
    # pysam uses Python's 0-based indexing
    # http://pysam.readthedocs.io/en/latest/faq.html#pysam-coordinates-are-wrong
    for match in c.finditer(line.cigarstring):
        n, op = int(match.group(1)), match.group(2)   
        # op .. operator, n .. number of bases

        if op == 'N' and n >= min_gap:
            yield pos , pos + n  
            # previous versions were
            # pos + 1, pos + n
            # this creates an off-by-1 error, verified w/
            # line.get_aligned_pairs(with_seq=True)
            pos += n

        elif op in ['M', 'D', 'N']:  # not I, mistake in journal 2017-12-28T1805
            pos += n

        else:
            pass


@contextlib.contextmanager
def smart_open(filename=None):
    '''
    stackoverflow.com/questions/17602878
    '''
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


@click.command()
@click.option(
    '--bam', required=True, help='Reads aligned to reference.')
@click.option(
    '-o', default=None, help='Where to write output.')
@click.option(
    '--header', is_flag=True, help='Include header.')
@click.option(
    '--maxn', default=6, type=int, help='More than this number of consecutive N bases is considered a splice site and will not contribute to error calculation.')
def identity(bam, o, header, maxn):
    '''Calculate per-base identity (1 - error rate).
    \b
    > The definition of 'identity' is the same as how BLAST defines it: the 
    number of matches in the alignment divided by the total bases in the 
    alignment (including gaps). -- Ryan Wick, 
    https://github.com/rrwick/Basecalling-comparison

    Usage:

    \b
    peek identity --bam test.bam --header -o eval.tsv
    head -n 4 eval.tsv
    # name length_query    length_reference    identity
    # 7b599a8b-69fc-4c87-9aea-1575d01223be    821 878 0.911
    # aa51cfbc-07de-4375-b11e-b97f8eaa4d67    856 919 0.944
    # f8514967-2988-45e5-9b39-33ad039f17d7    795 845 0.931
    '''
    with smart_open(o) as fh:  # fh .. file handle
        if header:
            print('{}\t{}\t{}\t{}'.format(
                'name', 'length_query', 'length_reference', 'identity'), file=fh)

        for line in tqdm(pysam.AlignmentFile(bam, 'rb')):
            cigar = line.cigarstring
            
            if not line.is_unmapped:
                matches = sum(
                    [int(i.strip('M')) for i in re.findall('\d*M', line.cigarstring)])
                
                # we explude the soft clipped (label: 4) subsequences from the 
                # calculation of alignment length -- seems reasonable
                aln_len = sum(
                    [count for label, count in line.cigartuples if label != 4])

                # in spliced mappings, we exclude large gaps (i.e. splice
                # junctions)
                for p in get_discont_pos(line, min_gap=maxn):
                    donor, acceptor = p
                    aln_len -= (acceptor - donor)

                identity = round(matches / aln_len, 3)
            
            print('{}\t{}\t{}\t{:.3f}'.format(
                line.qname, line.qlen, line.rlen, identity), file=fh)


if __name__ == '__main__':
    identity()

