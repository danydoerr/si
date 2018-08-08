#!/usr/bin/env python

from sys import stdout, stderr, exit, maxint
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from scipy import random
import logging

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

from matplotlib import pylab as plt
import numpy as np

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def constructGenome(n):
    """ construct a genome with n genes """
    return list(range(n))


def indel(genome, n, c):
    """ perform an event that alters the gene content """
    #common = [i for i, x in enumerate(genome) if x <= n]
    #p = random.choice(common)
    p = random.randint(0, len(genome))
    new_indel = genome[p] <= n
    genome[p] = c
    c +=1
    return c, p, new_indel


def hgt(genome, c):
    """ perform an Horzontal Gene Transfer (HGT) event """
    # pick a gene
    p = random.randint(0, len(genome))
    gene = genome.pop(p)
    # new insertion spot
    pp = random.randint(0, len(genome)+1)
    genome.insert(pp, gene)
    return c


def __count_hits__(A, B, k):
    """ count the number of common genes in the 2k-neighborhood of each gene of
    genome A w.r.t. genome B """

    g2pos = [-1] * (max(max(A), max(B))+1)
    
    for i, g in enumerate(A):
        g2pos[g] = i

    # total count of neighborhood matches
    s = 0

    hit_indel = 0
    s_comm  = 0.
    s_n = 0.
    for i, g in enumerate(B):
        if g2pos[g] > -1:
            NgB = set()
            # neighborhood in B
            for j in xrange(i-k, i+k+1):
                NgB.add(B[j % len(B)])

            # neighborhood in A
            si = 0
            for j in xrange(g2pos[g]-k,g2pos[g]+k+1):
                if A[j % len(A)] in NgB:
                    si += 1
            # remove count for i itself
            si -= 1
            s_comm = (s_n * s_comm + si)/(s_n+1)
            s_n += 1
            s += si
        else:
            hit_indel += 1

    return s, hit_indel, s_comm

def estimateIndels(n, mu, t):
    """ \rho(t) """
    return 1-(1-float(mu)/n)**t

def si(A, B, k):
    """ compute the Synteny Index distance of two genomes A and B """
    n = len(set(A+B))

    sAB, hit_indelAB, m_AB = __count_hits__(A, B, k)
    sBA, hit_indelBA, m_BA = __count_hits__(B, A, k)

    return (sAB + sBA) /float(4*k*n)


def estimateSI(n, lam, mu, k, t):
    """ \widetilde{SI} """
    it = estimateIndels(n, mu, t)
    return (1-np.exp(-(3*lam+2.*mu)*t/n)) 
    #return (1-np.exp(-(3*lam+2.*mu)*t/n)) * (n-1.-2.*k*(1-it))/(n-1.)

def approximateSI(n, lam, mu, k, t):
    """ \widehat{SI} """ 
    it = estimateIndels(n, mu, t)
    return 1-(1-it)**2


def evolve(genome, hgt_rate, indel_rate, time, c):
    """ evolves a genome along evolutionary time according to the given rates of
    HGT and indel events.
    
    The last parameter "c" denotes the next unused gene family ID count and is
    necessary to ensure that only new gene families are inserted in indel
    events.

    The method yields after every time step. 
    """

    orig = c-1
    for t in xrange(1, time+1):
        indel_events = random.poisson(indel_rate)
        new_indel_events = 0
        ps = list()
        for _ in xrange(indel_events):
            c, p, is_new = indel(genome, orig, c)
            new_indel_events += is_new and 1 or 0
            ps.append(p)
        if indel_events:
            LOG.debug(('%s indels occurring at time step %s at position(s) ' + \
                    '%s') %(indel_events, t, ','.join(map(str, sorted(ps)))))
        hgt_events = random.poisson(hgt_rate)
        for _ in xrange(hgt_events):
            c = hgt(genome, c)
        if hgt_events:
            LOG.debug('%s HGTs occurring at time step %s' %(hgt_events, t))

        yield t, indel_events, new_indel_events


def runExperiment(n, samples, time, hgt_rate, indel_rate, k):
    """ let {samples} many genome evolve along evolutionary time and record the
    SI distance at each time step"""

    # data matrix that maintains changing SI distance along evolutionary time
    data = np.empty((samples, time), dtype=float)
    LOG.info(('sampling %s genomes of size %s over %s mutational ' + \
            'steps...') %(samples, n, time))
    
    new_indels = np.empty(time)
    new_indels[0] = 0
    for s in xrange(samples):
        if not (s % max(1, samples/100)):
            LOG.info('%s%%' %((100*s)/samples))
        A = constructGenome(n)
        B = list(A)
        for t, indel_events, new_indel_events in evolve(B, hgt_rate,
                indel_rate, time, n):
            data[s,t-1] = si(A, B, k)
            if t-1:
                new_indels[t-1] = new_indels[t-2]
            new_indels[t-1] += new_indel_events
   
    return data, new_indels


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('n', type=int, help='number of genes')
    parser.add_argument('-g', '--hgtrate', type=float, default=1,
            help='rate of horzontal gene transfer (HGT) events')
    parser.add_argument('-i', '--indelrate',  type=float, default=0.1,
            help='rate of gene content modification events')
    parser.add_argument('k', type=int)
    parser.add_argument('-s', '--samples', type=int, default=1000,
            help='number of repeats/samples')
    parser.add_argument('-t', '--time', type=int, default=1000,
            help='number of steps of mutations')
    args = parser.parse_args()

     
    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    data, new_indels = runExperiment(args.n, args.samples, args.time,
            args.hgtrate, args.indelrate, args.k)

#    x = np.array(xrange(args.time))
#    plt.plot(x, new_indels, '.')
#    plt.plot(x, estimateIndels(args.n, args.i, x)*args.n)
#    plt.show()
    #it = estimateIndels(args.n, args.indelrate, args.time)
    #plt.plot(x, [(1.-(2.*args.k*(1-it))/(args.n-1.))*(1-it) + it] * args.time, label='upper bound')
    #plt.boxplot(1-data, labels=list(xrange(1, args.time+1)))
    #plt.plot(xrange(1, args.time+1), xrange(1./(args., args.time+1))

    LOG.info('plotting result')
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
    x = np.array(xrange(1, args.time+1))
    plt.plot(x, [1-np.median(data[:, i]) for i in xrange(data.shape[1])],
            label=r'median $SI$ over %s samples' %args.samples)
    plt.plot(x, approximateSI(args.n, args.hgtrate, args.indelrate, args.k, x),
            label=r'$\widetilde{SI}(G_0, G_t)$')
    plt.plot(x, estimateSI(args.n, args.hgtrate, args.indelrate, args.k, x),
            '--', label=r'$\widehat{SI}(G_0, G_t)$')
    title = r''
    if args.hgtrate:
        title += r'HGT $\gamma=%s$' % args.hgtrate

    if args.indelrate:
        if title:
            title += r', '
        title += r'indel $\mu=%s$' % args.indelrate
    plt.title(title)
    plt.xlabel('$t$ (time)')
    plt.ylabel('distance')
    plt.legend()

    plt.savefig(stdout, format='pdf')
    LOG.info('DONE')
