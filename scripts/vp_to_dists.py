#!/usr/bin/env python

from sys import stdout,stderr,argv,exit
from itertools import izip, combinations
from os.path import dirname, basename, join
from glob import glob
import re

from alf import readHomologies, readGenePositions, matrix2families


PAT_SPECIES = re.compile('SE(\d+)')

if __name__ == '__main__':
    if len(argv) != 2:
        print '\tusage: %s <VP MATRIX>' %argv[0]
        exit(1)

    out = stdout

    MATRIX = readHomologies(open(argv[1]))
    FAMS = matrix2families(MATRIX)
    DB_DIR = join(dirname(argv[1]), '..', 'DB')

    for f in glob(join(DB_DIR, '*_aa.fa')):
        sid = int(PAT_SPECIES.match(basename(f)).group(1))
        # read gene positions in ALF's fasta files, we do not need to know the
        # orientation of the genes inside the genome..
        ORIENT, POS, GNAME = izip(*readGenePositions(open(f)))
        POS_rev = [0] * len(FAMS[sid-1])
        for i, p in enumerate(POS[1:]):
            POS_rev[p] = i

        print >> out, '>SE%03i' %sid
        for p in POS_rev:
            print >> out, FAMS[sid-1][p],
        print >> out, ''
