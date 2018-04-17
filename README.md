## peek

A couple of small utilities to quickly peek at results from nanopore-based sequencing runs.

### Install

```bash
# create a virtual environment (e.g. using conda)
# cd into a directory where you store source code
git clone https://github.com/phiweger/peek
cd peek 
pip install -r requirements.txt -e .
```

### Usage

List subcommands:

```bash
peek
```

Look at the API of subcommands:

```bash
peek stats --help

# Usage: peek stats [OPTIONS]

#   Return basic stats from basecalled fastq.

#   stats: read name, length, median quality, label

#   Usage:

#   peek stats -f bar.fq -o bar.txt -l bad
#   # 3964e59a-f07d-46ff-bff5-1e2c03e4254a    981 10.0    bad
#   # c0dd59b3-bd29-488f-a78f-fa75d6fc6a3d    987 12.0    bad
#   # 3f6a3cd0-c788-40c2-900a-4ea3a6eb424c    869 13.0    bad

#   # only summary
#   peek stats --summary -f bar.fq > /dev/null

#   # or read from stdin
#   seqtk sample -s42 bar.fq 10000 | peek stats -f - -o bar.txt -l bad
#   cat workspace/pass/barcode01/* | peek stats --summary -f - > /dev/null

#   # concatenate multiple fastq stats and plot using R
#   peek stats --label bc01 -f BC01.fastq --summary > read_len_dist.tsv
#   peek stats --label bc02 -f BC02.fastq --summary >> read_len_dist.tsv
#   R

#   library(ggplot2)
#   library(readr)
#   df <- read_tsv('read_len_dist.tsv', col_names=c('name', 'length', 'quality', 'label'))
#   ggplot(df, aes(x=length, colour=label)) +
#       geom_histogram(bins=100) +
#       scale_x_log10(labels=scales::comma)

# Options:
#   --summary         Write summary to stderr
#   --header          Include header
#   -l, --label TEXT  Label
#   -n INTEGER        How many reads?
#   -o TEXT           Where to write output?
#   -f TEXT           Raw fastq file from e.g. Albacore basecaller  [required]
#   --help            Show this message and exit.
```
