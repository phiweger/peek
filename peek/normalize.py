'''
Aim: Draw random reads from fastq and stop when uniform coverage over a reference
genome achieved.

use:

- mappy
- brentp code
- zibra code

maybe give and argument like: have 90% of thing covered by 100x with a maximum of 
y bases

or:

- take a fq of debarcoded samples (i.e. this script will be able to run)
over the multiplexed samples in parallel

- map primers
- then map randomly drawn reads until we reach nbases target (don't choose best)
- choose 2 closest primers -- read assigned to those
- crop softclip from reads
- count primer pair +1 
- discard softclip (likely some remaining adapter)
'''

from collections import defaultdict
import csv
import numpy as np
import re
import sys

from Bio import SeqIO
import click
import mappy as mp
# https://github.com/lh3/minimap2/blob/master/python/minimap2.py
from skbio.sequence import DNA
# conda install matplotlib  # stackoverflow.com/questions/31373163
# otherwise catch RuntimeError and import again, odd, but solves problem
from skbio.alignment import local_pairwise_align_ssw
from tqdm import tqdm


def load_primer(fp):
    '''
    in: primal_scheme.csv from http://primal.zibraproject.org/

    Primer Name,Sequence,Pool,Length,Tm,GC%,Start,End
    ppv_1_LEFT,ACTCAACACAACATACAAAATTTTATGCG,1,29,60.53,31.03,11,40
    '''
    d = {}
    with open(fp, 'r') as file:
        f = csv.reader(file)
        header = next(f)
        for row in f:
            name, seq, pool, length, tm, gc, start, end = row
            d[name] = seq
    return d


def load_reference(fp):
    '''
    only support for one contig for now

    out: ('EU117116.1', 'AAAATATAAAAACT...')
    '''
    rname, rseq, _ = next(mp.fastx_read(fp))
    return rname, rseq


def get_primer_positions(primer_seqs, reference_seq):
    # hash map to hold start, stop positions for primers
    d = {}
    
    for p in primer_seqs.items():
        qname, qseq = p
        if 'RIGHT' in qname:  # mind the reverse complement
            qseq = str(DNA(qseq).reverse_complement())

        # align primer to reference using (striped) Smith-Waterman
        msa, aln_score, pos = local_pairwise_align_ssw(
            DNA(qseq), DNA(reference_seq))
        
        _, rpos = pos
        pstart, pend = rpos
        pspan = range(pstart, pend + 1)  # pspan .. primer span
        # + 1 bc/ the alignment is inclusive of last position while the fn 
        # range (Python in general) is not
    
        # contains start, end position of primer on ref
        d[pstart] = qname
        d[pend] = qname
    return d


def closest(l, pos):
    '''
    stackoverflow.com/questions/12141150
    '''
    return min(l, key=lambda x: abs(x - pos))


def get_number(p):
    return re.findall('\d+', p)[0]


def match(p1, p2):
    return re.findall('\d+', p1) == re.findall('\d+', p2)


def collect_valid_amplicons(
    fp_query, fp_reference, fp_out, primer_pos, coverage):
    '''
    in: 
    
    - dict w/ primer coordinates on reference (bed format?)
    - reads.fq
    - reference.fa

    out:

    - reads.valid.fq

    KM273015.1  11  40  ppv2000_1_LEFT  1
    KM273015.1  2160    2182    ppv2000_1_RIGHT 1
    '''
    # nreads = 100  # aimed for coverage
    counter = defaultdict(int)
    a = mp.Aligner(fp_reference, preset='map-ont')  # load or build index
    d = primer_pos

    # q = SeqIO.parse(query, 'fastq')
    # TODO: sampling logic has to go here, i.e. pbrent's code
    # https://github.com/brentp/bio-playground/blob/master/fileindex/examples/bench.py
    # https://github.com/zibraproject/zika-pipeline/blob/master/scripts/align_trim.py

    with open(fp_out, 'w+') as out:
        for name, seq, qual in tqdm(mp.fastx_read(fp_query)): 

            for hit in a.map(seq):
                if hit.is_primary:  
                # only non-softclipped parts of seq 
               
                    rstart = hit.r_st
                    rend = hit.r_en
                    pstart = closest(d.keys(), rstart)
                    pend = closest(d.keys(), rend)
                    p1 = d[pstart]
                    p2 = d[pend]
                    # print(hit.r_st, hit.r_en, closest(d.keys(), hit.r_st), closest(d.keys(), hit.r_en))

                    if match(p1, p2):  
                        # check plausibility, 137048 / 153418 = 0.89
                        # cleaning routine
                        # logic from:
                        # https://github.com/zibraproject/zika-pipeline/blob/master/scripts/align_trim.py
                        # if the alignment starts before the end of the primer, trim to that position
                        
                        if counter[get_number(p1)] < coverage:
                            counter[get_number(p1)] += 1  # no need to init, defaults to 0
                            # TODO: pop off list primer pair that reached 100
                            # will save time

                            # write seq w/o softclipped parts; remember to
                            # trim quality string as well
                           
                            # have a condition that checks for soft clip
                            out.write('@{}\n{}\n+\n{}\n'.format(
                                name, 
                                seq[hit.q_st:hit.q_en], 
                                qual[hit.q_st:hit.q_en]))

                            # stop going through reads once target coverage
                            # achieved for all primer pairs (amplicons)
                            if all([v == coverage for v in counter.values()]):
                                # print(counter)
                                # print(coverage)
                                sys.exit(0)

                            break  # next read in fastq



@click.command()
@click.option(
    '--query', required=True, help='Path to reads in fastq (adapters removed).')
@click.option(
    '--reference', required=True, help='Path to reference in fasta.')
@click.option(
    '--primers', default=None, help='Path to primers in csv.')
@click.option(
    '--out', default=None, help='Where to write output.')
@click.option(
    '--coverage', default=100, help='Desired subsampled coverage.')
def normalize(query, reference, primers, out, coverage):
    '''Sample a read set to normalize coverage accross amplicons.

    This script performs multiple tasks:

    \b
    1. Align primers to reference.
    2. Align reads to reference.
    3. Filter amplicons that can be assigned to one specific primer pair.
    4. Trim softclipped start, end.
    5. Return reads such that coverage is even accross reference.

    Usage

    PRIMERS='primer_amplicon_length1000_overlap200.primal.csv'
    TARGET='EU117116.1.fasta'
    QUERY='test.fq'  # porechopped version of '2017-11-30_amplicon_1000.fq' (default params)
    OUT='sample.norm.fq'

    \b
    nanopeek normalize \\
        --query $QUERY --reference $REF --out $OUT \\
        --coverage 50 --primers $PRIMERS
    '''
    
    p = load_primer(primers)
    rname, rseq = load_reference(reference)
    pos = get_primer_positions(p, rseq)
    collect_valid_amplicons(query, reference, out, pos, coverage)


if __name__ == '__main__':
    normalize()



