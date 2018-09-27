#!/usr/bin/env Rscript
# Rscript --vanilla sillyScript.R iris.txt out.txt

options(scipen=999)                 # stackoverflow.com/questions/5352099
options(warn=-1)                    # stackoverflow.com/questions/16194212
library(R.devices, quietly=TRUE)    # stackoverflow.com/questions/17348359
library(ggplot2)
library(readr)

args = commandArgs(trailingOnly=TRUE)
out <- args[1]


theme_minimal <-
    theme_classic() +
    theme(
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.spacing=unit(2, 'lines'),
        strip.background=element_blank(),
        )


df <- read.table(paste0(out, '/pores.csv'), header=T, sep=',')
df <- df[with(df, order(time)),]  # for cumulative sum on fragment length


pore_gel <-
    ggplot(df, aes(x=time, y=channel)) +
    geom_point(size=0.01, alpha=0.5) +
    theme_minimal +
    xlab('time [h]')
suppressGraphics(ggsave(
    paste0(out, '/qc_pore_gel.png'), height=12, width=12, units='cm'))


length_distribution <-
    ggplot(df, aes(x=length)) +
    geom_histogram(bins=20, fill='white', color='black') +
    theme_minimal +
    xlab('length') +
suppressGraphics(ggsave(
    paste0(out, '/qc_length_distribution.pdf'), height=6, width=6, units='cm'))


yield <-
    ggplot(df, aes(x=time, y=cumsum(length/1e6))) +
    geom_line() +
    theme_minimal +
    ylab('yield [mb]') +
    xlab('time [h]')
suppressGraphics(ggsave(
    paste0(out, '/qc_yield.pdf'), height=6, width=6, units='cm'))


df2 <- read.table(paste0(out, '/pores_active.csv'), header=T, sep=',')


active_pores <-
    ggplot(df2, aes(x=time, y=active)) +
    geom_line() +
    theme_minimal +
    scale_y_continuous(limits=c(0, 1)) +
    xlab('time [h]') +
    ylab('active pores, total')
suppressGraphics(ggsave(
    paste0(out, '/qc_pores_active.pdf'), height=6, width=6, units='cm'))

