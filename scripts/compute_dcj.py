#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import combinations, chain
from sys import stderr, stdout, exit
import logging

import networkx as nx
import numpy as np

CHR_CIRCULAR = ')'
CHR_LINEAR = '|'
ORIENT_POSITIVE = '+'
ORIENT_NEGATIVE = '-'
EXTREMITY_TAIL = 't'
EXTREMITY_HEAD = 'h'

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def readGenomes(data, genomesOnly=None):
    """Read genome in UniMoG format
    (https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual)"""

    res = list()

    # helper function for parsing each individual gene
    str2gene = lambda x: x.startswith(ORIENT_NEGATIVE) and (ORIENT_NEGATIVE, \
            x[1:]) or (ORIENT_POSITIVE, x.lstrip(ORIENT_POSITIVE))
    # process each line, assuming that the file is well-formatted
    skip = False
    for line in data:
        line = line.strip()
        if line:
            if line.startswith('>'):
                genomeName = line[1:].strip()
                if genomesOnly == None or genomeName in genomesOnly:
                    skip = False
                    res.append((genomeName, list()))
                elif genomesOnly:
                    skip = True
            elif not skip:
                structure = None
                genes = None
                if line[-1] in (CHR_CIRCULAR, CHR_LINEAR):
                    structure = line[-1]
                    genes = map(str2gene, line[:-1].split())
                else:
                    struture = CHR_LINEAR
                    genes = map(str2gene, line.split())
                res[-1][1].append((structure, list(genes)))
    return res


def constructAdjacencyGraph(genomes):
    """ Breakpoint graph over all genomes.
    Colors: 0 = gene, 1 = 1st genome, 2 = 2nd genome, etc """

    G = nx.MultiGraph()

    gene2ext = lambda x: x[0] == ORIENT_POSITIVE and ((EXTREMITY_TAIL, x[1]), \
            (EXTREMITY_HEAD, x[1])) or ((EXTREMITY_HEAD, x[1]), \
            (EXTREMITY_TAIL, x[1]))


    ext2node = dict()
    for i, (_, chromosomes) in enumerate(genomes):
        for chr_type, genes in chromosomes:
            # find initial extrimity
            prev = None
            if chr_type == CHR_LINEAR:
                # telomere is defined as tuple, corresponding to adjoining
                # extremity
                prev = (gene2ext(genes[0])[0], EXTREMITY_HEAD)
            else:
                prev = gene2ext(genes[-1])[1]
            # extend
            for g in genes:
                ext1, ext2 = gene2ext(g)
                v = (i, ) + tuple(sorted((prev, ext1)))
                G.add_node(v)
                for ext in (prev, ext1):
                    if ext not in ext2node:
                        ext2node[ext] = list()
                    ext2node[ext].append(v)
                prev = ext2

            if chr_type == CHR_LINEAR:
                telomere = (prev, EXTREMITY_HEAD)
                v = (i, ) + tuple(sorted((prev, telomere)))
                G.add_node(v)
                for ext in (prev, telomere):
                    if ext not in ext2node:
                        ext2node[ext] = list()
                    ext2node[ext].append(v)
    for ext, vertices in ext2node.items():
        for u, v in combinations(vertices, 2):
            G.add_edge(u, v, extremity=ext)
    return G


def calcDCJ(G, n):
    c = 0
    i = 0
    for C in nx.connected_components(G):
        degrees = tuple(sorted(set(dict(G.degree(C)).values())))
        l = len(G.edges(C))
        if degrees == (2,):
            c += 1
        elif degrees == (1,2) or degrees == (1, ):
            if l %2:
                i += 1
        else:
            raise Exception('There is an error in the program, the adjacency ' + \
                    'graph has components that are neither simple cycles ' + \
                    'nor simple paths')
    return n - c - i/2


def caldGenomeSize(genome):

    res = 0
    for _, genes in genome:
        res += len(genes)

    return res


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('genomes', type=open,
            help='Genome sequences in UniMoG format')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    genomes = readGenomes(args.genomes)

    if len(genomes) < 2:
        LOG.fatal('Expected at least two genomes, received %s, exiting' %len(genomes))
        exit(1)

    n = caldGenomeSize(genomes[0][1])
    D = np.zeros((len(genomes), len(genomes)), dtype=int)
    for i, j in combinations(range(len(genomes)), 2):
        G = constructAdjacencyGraph((genomes[i], genomes[j]))
        D[i,j] = D[j, i] = calcDCJ(G, n)

    gNames, _ = zip(*genomes)
    print('\t'.join(chain(('', ), gNames)))
    for i, name in enumerate(gNames):
        print('\t'.join(chain((name, ), map(str, D[i, :]))))
