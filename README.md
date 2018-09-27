## peek

A couple of small utilities to quickly peek at results from nanopore-based sequencing runs.


### Install

```bash
# create a virtual environment (e.g. using conda)
# cd into a directory where you store source code
git clone https://github.com/phiweger/peek
cd peek 
pip install -e .
```

For plotting in R there are 3 dependencies:

```r
# R
install.packages('ggplot2')
install.packages('R.devices')
install.packages('readr')
```


### Usage

List subcommands:

```bash
peek
```

Look at the API of subcommands:

```
peek stats --help

Usage: peek stats [OPTIONS]

  Return basic stats from basecalled fastq.

  stats: read name, length, median quality, label

  Usage:

  peek stats -f bar.fq -o bar.txt -l bad
  # 3964e59a-f07d-46ff-bff5-1e2c03e4254a    981 10.0    bad
  # c0dd59b3-bd29-488f-a78f-fa75d6fc6a3d    987 12.0    bad
  # 3f6a3cd0-c788-40c2-900a-4ea3a6eb424c    869 13.0    bad

  # only summary
  peek stats --summary -f bar.fq > /dev/null

  # or read from stdin
  seqtk sample -s42 bar.fq 10000 | peek stats -f - -o bar.txt -l bad
  cat workspace/pass/barcode01/* | peek stats --summary -f - > /dev/null

  # concatenate multiple fastq stats and plot using R
  peek stats --label bc01 -f BC01.fastq --summary > read_len_dist.tsv
  peek stats --label bc02 -f BC02.fastq --summary >> read_len_dist.tsv

  # vis w/ R
  library(ggplot2)
  library(readr)
  df <- read_tsv(
    'read_len_dist.tsv', col_names=c('name', 'length', 'quality', 'label'))
  ggplot(df, aes(x=length, colour=label)) +
      geom_histogram(bins=100) +
      scale_x_log10(labels=scales::comma)

Options:
  --summary         Write summary to stderr
  --header          Include header
  -l, --label TEXT  Label
  -n INTEGER        How many reads?
  -o TEXT           Where to write output?
  -f TEXT           Raw fastq file from e.g. Albacore basecaller  [required]
  --help            Show this message and exit.
```


## peek offshore

The idea here is to not wrap command line tools w/ arbitrary python code, but to use [Snakemake](http://snakemake.readthedocs.io/en/latest/) behind the command line interface of a Python package, as e.g. illustrated [here](https://github.com/ctb/2018-snakemake-cli) and described in more detail [here](http://ivory.idyll.org/blog/2018-workflows-applications.html).

Initially, only _host depletion_ is implemented. This means that we provide a reference genome and only care about reads that are not (or are) in this genome. For example in the following example data set, we only care for viral sequences but are not interested in the human reads present:

- [Ebola data set](https://www.ncbi.nlm.nih.gov/sra/SRX3544109[accn]) -- a mixture of human cells and virus
- [Ebola reference genome](https://www.ncbi.nlm.nih.gov/nuccore/LT605058.1)

```bash
cd peek/workflow/
peek offshore \
    --paramsfile params-deplete.json \
    --workflowfile workflow-deplete.json \
    --snakefile Snakefile \
    -o ~/tmp/results
```



