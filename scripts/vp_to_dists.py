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
    DB_DIR = join(dirname(argv[1]), '..', 'DB')

    POS_rev = list()
    GNAME   = list()
    SID     = list()

    for f in glob(join(DB_DIR, '*_aa.fa')):
        sid_f = int(PAT_SPECIES.match(basename(f)).group(1))
        SID.append(sid_f)
        # read gene positions in ALF's fasta files, we do not need to know the
        # orientation of the genes inside the genome..
        ORIENT_f, POS_f, GNAME_f = izip(*readGenePositions(open(f)))
        POS_rev_f = [0] * (max(POS_f)+1)
        for i, p in enumerate(POS_f[1:]):
            POS_rev_f[p] = i
        POS_rev.append(POS_rev_f)
        GNAME.append(GNAME_f)

    FAMS = matrix2families(MATRIX, POS_rev, GNAME)

    for i, POS_rev_i in enumerate(POS_rev):
        sid = SID[i]
        print >> out, '>SE%03i' %sid
        for p in POS_rev_i:
            print >> out, FAMS[sid-1][p],
        print >> out, ''
