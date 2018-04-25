#! /usr/bin/env python
"""
Execution script for snakemake workflows.

Shamelessly adopted from Titus Brown, https://github.com/ctb/2018-snakemake-cli
"""

import os.path
from snakemake import snakemake
# import sys
# import pprint
import json

import click


thisdir = os.path.abspath(os.path.dirname(__file__))


@click.option(
    '--snakefile', '-s', required=False, type=click.Path(exists=True), help='lalala', default=os.path.join(thisdir, 'Snakefile'))
@click.option(
    '--workflowfile', '-w', required=True, type=click.Path(exists=True), help='lalala')
@click.option(
    '--paramsfile', '-p', required=True, type=click.Path(exists=True), help='lololo')
@click.option(
    '--dryrun', '-n', default=False, required=False, is_flag=True, help='lululu')
@click.command()
def bar(snakefile, workflowfile, paramsfile, dryrun):

    with open(workflowfile, 'rt') as fp:
        workflow_info = json.load(fp)

    target = workflow_info['workflow_target']
    # pass these over the command line or in the params file
    # see https://bitbucket.org/snakemake/snakemake/issues/98/allow-passing-variables-to-the-workflow
    # snakemake(myworkflow, config={"var1": 1, "var2": "foo", "var3": [1,2,3]})
    config = {"query": "reads", "reference": "hucov229e"}

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

    if status:  # translate "success" into shell exit code of 0
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
