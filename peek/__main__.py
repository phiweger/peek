'''
peek command line.
'''


import click
# https://kushaldas.in/posts/building-command-line-tools-in-python-with-click.html


from peek.stats import stats
from peek.normalize import normalize
from peek.identity import identity
from peek.filtlen import filtlen
from peek.chop import chop

# test for snakemake file integration
from peek.foo import bar


@click.group()
def cli():
    pass


cli.add_command(stats)
cli.add_command(normalize)
cli.add_command(identity)
cli.add_command(filtlen)
cli.add_command(chop)

# test for snakemake file integration
cli.add_command(bar)
