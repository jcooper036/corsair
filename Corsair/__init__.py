#!/usr/bin/env python3
"""init for the Corsair module"""

## classes
from Corsair.classes.Isoform import Isoform
from Corsair.classes.Ctlog import Ctlog

## wrapper definintions
from Corsair.wrappers.corsair_initialize import corsair_initialize
from Corsair.wrappers.corsair_execs import corsair_execs
from Corsair.wrappers.corsair_results import corsair_results

## definintions
from Corsair.defs.translate import translate
from Corsair.defs.read_fasta import read_fasta
from Corsair.defs.sample_ctl import sample_ctl
from Corsair.defs.load_ctl import load_ctl
from Corsair.defs.initialize_isoforms import initialize_isoforms
from Corsair.defs.load_isoform import load_isoform
from Corsair.defs.save_isoform import save_isoform
