from collections import defaultdict
from datetime import datetime
from itertools import islice
from operator import itemgetter
import os
import subprocess
import sys

from Bio import SeqIO
import click
import numpy as np
import maya
import pandas as pd
import peek
from tqdm import tqdm


def create_folder(directory):
    '''
    https://gist.github.com/keithweaver/562d3caa8650eefe7f84fa074e9ca949
    '''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
        sys.exit(-1)


@click.option(
    '-i', '--infile', required=True, 
    help='Raw fastq file from e.g. Albacore basecaller')
@click.option(
    '-o', '--outdir', required=True, default=None, 
    type=click.Path(exists=False, resolve_path=True), 
    help='Output directory')
@click.option(
    '-f', '--force', required=False, is_flag=True, 
    help='Override output directory')
@click.option(
    '-p', '--plot', required=False, is_flag=True, 
    help='Plot results?')
@click.command()
def pores(infile, outdir, force, plot):
    '''Check run yield over time.

    peek pores -f -p -i zebra.fastq -o results
    '''

    print(outdir)

    if os.path.exists(outdir) and not force:
        print('Cannot clobber output directory w/o --force.')
        sys.exit(-1)


    create_folder(outdir)

    fq = SeqIO.parse(infile, format='fastq')
    d = defaultdict(list)
    
    mn = maya.now()
    length = []
    records = []


    print('Calculating read stats ...')

    for i in tqdm(fq):
        uid, runid, r, ch, start = i.description.split(' ')  # ch .. channel
        ch = ch.split('=')[-1]  # ch=170
        t = maya.parse(start.split('=')[-1])  # start_time=2017-10-01T02:05:51Z
        l = len(i.seq)
    
        if t < mn:  # find the time the run started at
            mn = t
    
        length.append(l)
        # for length summary, could be replaced by probabilistic data structure
    
        d[ch].append(t)
        records.append([t, ch, l])
    
    
    df = pd.DataFrame.from_records(records)
    df.columns = 'time channel length'.split(' ')
    
    hours = []
    for i in df['time']:
        delta = i.datetime() - mn.datetime()
        hours.append(delta.total_seconds()/3600)
    
    df['time'] = hours
    df.to_csv(outdir + '/pores.csv', index=None)
    

    print('\nSummary:')
    print('{:<10}{} gb'.format('yield:', round(sum(length)/1e9, 3)))
    print('{:<10}{} nt'.format('longest:', max(length)))
    print('{:<10}{:.0f} nt'.format('median:', np.median(length), 0))
    # print('\nDone.')
    
    
    # Kaplam-Meier curve of active channels.
    # Get the maximum time for each channel and sort channels accordingly.
    dd = {}
    for k, v in d.items():
        dd[k] = max(v)
    sorted_dd = sorted(dd.items(), key=itemgetter(1))  # stackoverflow, 613183
    # ...
    # ('283', '2017-10-01T06:48:21Z'),
    # ('379', '2017-10-01T07:43:06Z'),
    # ('464', '2017-10-01T08:03:37Z'),
    # ...
    
    
    with open(outdir + '/pores_active.csv', 'w+') as out:
        total = len(d.keys())  
        # there are at most 512 channels that could be  active
    
        out.write('time,active\n')  # header csv
        out.write('{},{}\n'.format(0, total/512))  
        # 0 .. starting point in hours
    
        for i in sorted_dd:
            _, t = i
            delta = t.datetime() - mn.datetime()
            total -= 1
            out.write('{},{}\n'.format(delta.total_seconds()/3600, total/512))


    if plot:
        # get the path to the R plotting script
        script_path = os.path.dirname(os.path.realpath(peek.__file__))
        # call Rscript and put all non-plot stuff into /dev/null
        command = f'''
            Rscript --vanilla {script_path}/scripts/pores.R {outdir} \
            > /dev/null 2>&1
            '''
        # click.echo(f'\nPlotting, command:\n{command}')
        print('\nPlotting ...')
        os.system(command)


    print('Done. Now go and prosper.')
