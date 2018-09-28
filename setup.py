from setuptools import setup

setup(
    name="peek",
    version='0.1',
    install_requires=[
        'Click>=6.7',
        'biopython>=1.70',
        'maya',
        'numpy<1.14.0',
        'tqdm>=4.7.2',
        'mappy',
        'scikit-bio',
        'pandas<0.23.0',
        'pysam',
    ],
    entry_points={
        'console_scripts': [
            'peek = peek.__main__:cli'
        ]})

# 'snakemake',