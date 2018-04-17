#!/usr/bin/env python3


import contextlib
import itertools
import sys
import warnings

from Bio import SeqIO
import click
import numpy as np
from tqdm import tqdm


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


@click.option(
    '-f', required=True, help='Raw fastq file from e.g. Albacore basecaller')
@click.option(
    '-o', default=None, help='Where to write output?')
@click.option(
    '-n', default=None, type=int, help='How many reads?')
@click.option(
    '-l', '--label', required=False, help='Label')
@click.option(
    '--header', is_flag=True, help='Include header')
@click.option(
    '--summary', is_flag=True, help='Write summary to stderr')
@click.command()
def stats(f, n, o, label, header, summary):
    '''Return basic stats from basecalled fastq.

    stats: read name, length, median quality, label

    Usage:

    \b
    peek stats -f bar.fq -o bar.txt -l bad
    # 3964e59a-f07d-46ff-bff5-1e2c03e4254a    981 10.0    bad
    # c0dd59b3-bd29-488f-a78f-fa75d6fc6a3d    987 12.0    bad
    # 3f6a3cd0-c788-40c2-900a-4ea3a6eb424c    869 13.0    bad

    \b
    # only summary
    peek stats --summary -f bar.fq > /dev/null

    \b
    # or read from stdin
    seqtk sample -s42 bar.fq 10000 | peek stats -f - -o bar.txt -l bad
    cat workspace/pass/barcode01/* | peek stats --summary -f - > /dev/null

    \b
    # concatenate multiple fastq stats and plot using R
    peek stats --label bc01 -f BC01.fastq --summary > read_len_dist.tsv
    peek stats --label bc02 -f BC02.fastq --summary >> read_len_dist.tsv
    R

    \b
    library(ggplot2)
    library(readr)
    df <- read_tsv('read_len_dist.tsv', col_names=c('name', 'length', 'quality', 'label'))
    ggplot(df, aes(x=length, colour=label)) + 
        geom_histogram(bins=100) + 
        scale_x_log10(labels=scales::comma)
    '''
    counter = 0

    if f == '-':
        fastq = SeqIO.parse(sys.stdin, 'fastq')
        click.echo('Streaming reads ...', file=sys.stderr)
    else:
        fastq = SeqIO.parse(f, 'fastq')
        click.echo('Evaluating reads ...', file=sys.stderr)


    longest = (0, '')
    collect = []
    discard = 0
    nbases = 0

    with smart_open(o) as fh:  # fh .. file handle

        if header:
            print('{}\t{}\t{}\t{}'.format(
                'name', 'length', 'quality', 'label'), file=fh)
        
        for r in tqdm(itertools.islice(fastq, n)):  # if n == None, defaults to all
            if len(r) == 0:
                # This can happen when a bam file is converted to fastq,
                # e.g. by bedtools bamtofastq. The use case here is to
                # look at the read length distribution of a subset of reads
                # that map some target.
                discard += 1
                continue
            phred = r.letter_annotations['phred_quality']
            qual = np.median(phred)
            print('{}\t{}\t{}\t{}'.format(r.name, len(r), qual, label), file=fh)

            nbases += len(r)

            if longest[0] < len(r):
                longest = (len(r), r.name)
            collect.append(len(r))
            
            counter += 1


    # click.echo('Found {} reads.'.format(counter), file=sys.stderr)
    # tqdm shows how many reads were taken into account
    if summary:
        click.echo('gigabases\t{:>5.3f}'.format(round(nbases/1e9, 3)), file=sys.stderr)
        click.echo('longest\t\t{} ({})'.format(*longest), file=sys.stderr)
        click.echo('mean\t\t{:.0f}'.format(round(np.mean(collect)), 0), file=sys.stderr)
        click.echo('std\t\t{:.0f}'.format(np.std(collect)), file=sys.stderr)
        click.echo('25%\t\t{:.0f}'.format(np.percentile(collect, q=25)), file=sys.stderr)
        click.echo('50%\t\t{:.0f}'.format(np.percentile(collect, q=50)), file=sys.stderr)
        click.echo('75%\t\t{:.0f}'.format(np.percentile(collect, q=75)), file=sys.stderr)
    
    if discard:
        click.echo('Discarded {} reads w/ length 0.'.format(discard), file=sys.stderr)
    click.echo('Done.', file=sys.stderr)


if __name__ == '__main__':
    stats()
