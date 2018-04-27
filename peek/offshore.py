#! /usr/bin/env python
'''
Execution script for snakemake workflows.

Shamelessly adopted from Titus Brown, https://github.com/ctb/2018-snakemake-cli
'''

import os.path
from snakemake import snakemake
# import sys
# import pprint
import json

import click


thisdir = os.path.abspath(os.path.dirname(__file__))
# print(__file__)  # stackoverflow.com/questions/779495
# print(thisdir)   # the directory of offshore.py

# default to thisdir + '/workflow/'


@click.option(
    '--snakefile', '-s', required=False, type=click.Path(exists=True), help='lalala', default=os.path.join(thisdir, 'Snakefile'))
@click.option(
    '--workflowfile', '-w', required=True, type=click.Path(exists=True), help='lalala')
@click.option(
    '--paramsfile', '-p', required=True, type=click.Path(exists=True), help='lololo')
@click.option(
    '--dryrun', '-n', default=False, required=False, is_flag=True, help='lululu')
@click.option(
    '--outdir', '-o', default='.', help='results directory')
@click.command()
def offshore(snakefile, workflowfile, paramsfile, dryrun, outdir):
    '''
    The offshore subcommand diverts calls to workflows organisized via
    Snakemake. You can specify the target in the workflow.json like so:

    {'workflow_target': 'deplete'}

    Available workflows:

    - deplete .. remove sequences that map to a certain reference

    Usage:

    \b
    cd peek/workflow/
    peek offshore \\
        --paramsfile params-deplete.json \\
        --workflowfile workflow-deplete.json \\
        --snakefile Snakefile \\
        -o ~/tmp/results
    '''

    with open(workflowfile, 'rt') as fp:
        workflow_info = json.load(fp)

    target = workflow_info['workflow_target']
    # pass these over the command line or in the params file
    # see https://bitbucket.org/snakemake/snakemake/issues/98/allow-passing-variables-to-the-workflow
    # snakemake(myworkflow, config={'var1': 1, 'var2': 'foo', 'var3': [1,2,3]})

    # config = {
    #     'query': 'reads',
    #     'reference': 'hucov229e',
    #     'outfile': out}
    config = {'outdir': outdir}

    print('--------')
    print('details!')
    print('\tsnakefile: {}'.format(snakefile))
    print('\tconfig: {}'.format(workflowfile))
    print('\tparams: {}'.format(paramsfile))
    print('\ttarget: {}'.format(target))
    print('--------')

    # run!!
    status = snakemake(
        snakefile,
        configfile=paramsfile,
        targets=[target],
        printshellcmds=True,
        dryrun=dryrun,
        config=config)

    if status:  # translate 'success' into shell exit code of 0
        return 0
    return 1


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='run snakemake workflows', usage='''run <workflow> <parameters> [<target>]

# Run snakemake workflows, using the given workflow name & parameters file.

# ''')

#     parser.add_argument('workflowfile')
#     parser.add_argument('paramsfile')
#     parser.add_argument('-n', '--dry-run', action='store_true')
#     args = parser.parse_args()

#     sys.exit(bar(args))
