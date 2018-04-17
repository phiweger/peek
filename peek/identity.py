import contextlib
import re
import sys

import click
import pysam


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
def identity(bam, o, header):
    '''
    \b
    > The definition of 'identity' is the same as how BLAST defines it: the 
    number of matches in the alignment divided by the total bases in the 
    alignment (including gaps). -- Ryan Wick, 
    https://github.com/rrwick/Basecalling-comparison

    Usage:

    nanopeek identity --bam test.bam --header -o eval.tsv 
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

        for line in pysam.AlignmentFile(bam, 'rb'):
            cigar = line.cigarstring
            
            if not line.is_unmapped:
                matches = sum(
                    [int(i.strip('M')) for i in re.findall('\d*M', line.cigarstring)])
                
                # we explude the soft clipped (4) subsequences from the 
                # calculation of alignment length -- seems reasonable
                aln_len = sum(
                    [count for label, count in line.cigartuples if label != 4])

                identity = round(matches / aln_len, 3)
            
            print('{}\t{}\t{}\t{:.3f}'.format(
                line.qname, line.qlen, line.rlen, identity), file=fh)


if __name__ == '__main__':
    identity()

