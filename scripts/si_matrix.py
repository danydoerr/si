#!/usr/bin/env python

from sys import stdout, stderr, exit, maxint
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import combinations, chain
import logging

import numpy as np

from vp_to_dists import PAT_SPECIES
from sampleGenome import  si


LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def inverseSI(si, lam, n, k):

    x = max(0, 1.-si/(1-2.*k/(n-1)))
    return min(-n/(lam*(3.-(5.*k)/(n-1))) * np.log(x), 6*n)


def readGenomes(data):

    res = list()
    for line in data:
        if line.startswith('>'):
            res.append((line[1:].strip(), list()))
        elif line.strip():
            res[-1][1].extend(map(int, line.strip().split(' ')))

    return res


def constructSIDistMat(genomes, k, lam):

    res = np.zeros((len(genomes), len(genomes)))

    # assume that all genomes have equal size
    n = len(genomes[0])

    for i, j in combinations(xrange(len(genomes)), 2):
        x = si(genomes[i], genomes[j], k)
        val = inverseSI(1-x, lam, n, k)
        res[i, j] = res[j, i] = val
    return res


def writeDistMat(gNames, D, out):

    print >> out, '\t'.join(chain(('', ), gNames))

    for i, name in enumerate(gNames):
        print >> out, '\t'.join(chain((name, ), map(str, D[i, :])))


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('genomes', type=file,
            help = 'Gene order data of all genomes in fasta-like format')
    parser.add_argument('k', type=int,
            help = 'SI parameter \"neighborhood size\" k')
    parser.add_argument('-l', '--lambdda', type=float, default=1.0,
            help = 'SI parameter \"lambda\"')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    LOG.info('read genome data')
    genomeNames, genomes = zip(*readGenomes(args.genomes))

    LOG.info('constructing distance matrix')
    D = constructSIDistMat(genomes, args.k, args.lambdda)
    LOG.info('writing matrix to standard out')
    writeDistMat(genomeNames, D, stdout)
    LOG.info('DONE')


