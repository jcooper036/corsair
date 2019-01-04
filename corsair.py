#!/usr/bin/env python3
import Corsair as cor

sample_ctl = '/Users/Jacob/corsair/Corsair/test_files/ctl_test.txt'

## parse the ctl file
ctl = cor.load_ctl(sample_ctl)

## for the first time
cor.corsair_initialize(ctl)

iso = 'gene1'

## run all the execs on a gene by gene basis - allows for parallelizing this part
cor.corsair_execs(ctl, iso)

## run all the results gathering functions for the whole gene list
cor.corsair_results(ctl, iso)