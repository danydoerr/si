#!/usr/bin/env python

from sys import stdout, stderr, exit, maxsize
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from functools import partial
from scipy import random
import logging

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

from matplotlib import pylab as plt
import numpy as np


from sampleGenome import constructGenome, si, hgt

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def calcIndelDistribution(n, t):

    dist = [1]
    for t in range(t):
        tmp = [0] * (t+2)
        for i, p in enumerate(dist):
            p0=1-float(t-i)/n
            p1 = 1-p0 
            tmp[i] += p0 * p
            tmp[i+1] += p1 * p
        yield tmp
        dist = tmp
        
def countHits(A, B, k):
    """ count the number of common genes in the 2k-neighborhood of each gene of
    genome A w.r.t. genome B """

    g2pos = [-1] * (max(max(A), max(B))+1)
    
    for i, g in enumerate(A):
        g2pos[g] = i

    # total count of neighborhood matches
    s = list()

    hit_indel = 0
    for i, g in enumerate(B):
        si = 0
        if g2pos[g] > -1:
            NgB = set()
            # neighborhood in B
            for j in range(i-k, i+k+1):
                NgB.add(B[j % len(B)])

            # neighborhood in A
            for j in range(g2pos[g]-k,g2pos[g]+k+1):
                if A[j % len(A)] in NgB:
                    si += 1
            # remove count for i itself
            si -= 1
        else:
            hit_indel += 1
        s.append(si)

    return s, hit_indel


def runExperimentOneStep(k):
    """ let interval evolve one step and record the SI distance"""

    # data structure that maintains SI distance in experimental runs
    data = [None] * (2*k-1)
    
    A = constructGenome(4*k)
    for pp in range(k+1, 3*k):
        B = list(A)
        # always choose gene at position k
        gene = B.pop(k)
        B.insert(pp, gene)

        LOG.debug('moving gene from position %s to %s' %(k, pp))
        hitdist, _ = countHits(A[:pp+k+1], B[:pp+k+1], k)
        data[pp-k-1] = map(lambda x: 2*k-x, hitdist)
   
    return data

def expectedIndelSI(n, k, indels):
    return 
    

def calculateDistOneStep(n, k):

    res = np.empty((n, ), dtype=float)
    # distance in k region
    for l in range(1, k+1):
        p = (l-1)*2
        res[p:p+2] = 2*l

    for l in range(1, k):
        p = (k+(l-1))*2
        res[p:p+2] = 2*k+l
    res[4*k-2:] = 3*k
    res /= k*n
    return res

def makeStat(interval, k):

    res = np.zeros((2*k+1, )) 

    for v in interval:
        res[v] += 1

    return res

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('n', type=int, help='number of genes')
    parser.add_argument('t', type=int, help='time (number of steps)')
#    parser.add_argument('k', type=int)
#    parser.add_argument('-s', '--samples', type=int, default=1000,
#            help='number of repeats/samples')
    args = parser.parse_args()
     
    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    n = args.n
    dist = np.zeros((args.t, args.n))
    for t, t_dist in enumerate(calcIndelDistribution(n, args.t)):
        dist[t,n-t-1:] = np.array(t_dist[:-1])

#    dist = calculateDistOneStep(args.n, args.k)
#    sample_from_analytic = np.array([np.random.choice(dist) for _ in
#        range(args.samples)])
#   
#    sample_from_simulated = np.empty((args.samples, ))
#    A = constructGenome(args.n)
#    for i in range(args.samples):
#        sample_from_simulated[i] = 1-si(A, hgt(list(A)), args.k)
#    
#    plt.hist([sample_from_analytic, sample_from_simulated])
#    plt.show()
#    data = runExperimentOneStep(args.k)
#    stats = np.array(map(partial(makeStat, k=args.k), data))
#
#    x = np.linspace(0, 2*args.k-1, 2*args.k)
#    plt.plot(np.array(map(sum, data))/(2.0*args.k), 'x')
#    plt.plot(x, yp[0:4*args.k:2]*args.n, 'o')
   
#    plt.figure()
#    x = np.array(map(sum, data) + [6.*args.k])/(2.0*args.k)
#    y = np.ones((2*args.k, ))
#    y[:-1] = 2
#    y[-1] = args.n-(4*args.k-2)
#    plt.bar(x, y)
#
#    plt.show()
