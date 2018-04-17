import contextlib
import sys

from Bio import SeqIO
import click
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


@click.command()
@click.option(
    '-f', required=True, help='Some fastq.')
@click.option(
    '-o', default=None, help='Where to write output.')
@click.option(
    '-l', '--low', required=True, type=int, help='Minimum length to filter.')
@click.option(
    '-h', '--high', required=True, type=int, help='Maximum length to filter.')
def filtlen(f, o, low, high):
    '''
    Filter reads based on length.
    '''
    if f == '-':
        fastq = SeqIO.parse(sys.stdin, 'fastq')
        click.echo('Streaming reads ...', file=sys.stderr)
    else:
        fastq = SeqIO.parse(f, 'fastq')
        click.echo('Evaluating reads ...', file=sys.stderr)

    with smart_open(o) as fh:  # fh .. file handle
        for record in tqdm(fastq):
            if low <= len(record) <= high:
                SeqIO.write(record, fh, 'fastq')


if __name__ == '__main__':
    filtlen()
