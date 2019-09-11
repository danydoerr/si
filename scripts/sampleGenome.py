#!/usr/bin/env python

from sys import stdout, stderr, exit, maxint
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from scipy.optimize import newton
from scipy import random
import logging

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

from matplotlib import pylab as plt
import numpy as np

NEWTON_MAX_ITER = 100000

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def constructGenome(n):
    """ construct a genome with n genes """
    return list(range(n))


def indel(genome, n, c, random_pos=False):
    """ perform an event that alters the gene content """
    #common = [i for i, x in enumerate(genome) if x <= n]
    #p = random.choice(common)
    p = random.randint(0, len(genome))
    new_indel = genome[p] <= n
    if random_pos:
        gene = genome.pop(p)
        # new insertion spot
        pp = random.randint(0, len(genome)+1)
        genome.insert(pp, c)
    else:
        genome[p] = c
    c +=1
    return c, p, new_indel


def hgt(genome):
    """ perform an Horzontal Gene Transfer (HGT) event """
    # pick a gene
    p = random.randint(0, len(genome))
    gene = genome.pop(p)
    # new insertion spot
    pp = random.randint(0, len(genome)+1)
    genome.insert(pp, gene)
    return genome


def countHits(A, B, k):
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

def si(A, B, k):
    """ compute the Synteny Index distance of two genomes A and B """
    n = len(set(A+B))

    sAB, hit_indelAB, m_AB = countHits(A, B, k)
    #sBA, hit_indelBA, m_BA = countHits(B, A, k)

    return sAB/(2.*k*n)

def estimateSI(n, mu, k, t):
    """ \widetilde{SI} """
    return (1-np.exp(-3./n*t)) * (1-2.*k/(n-1)*(np.exp(-3./n*mu*t)))
    #return (1-np.exp(-(3*lam+2.*mu)*t/n))


def derivativeSI(n, mu, k, t):
    return np.exp(-3./n*mu*t)*(3./n*mu)*(1-np.exp(-3./n*t)) + \
            np.exp(-3./n*t)*(3./n*t)* \
            ((n-1.)/(2.*k)-np.exp(-3./n*mu*t))

def sndDerivativeSI(n, lam, mu, k, t):

    return 9.*((lam+mu)**2/(2.*k*n)*np.exp(3./n*mu*t) + (mu/float(n))**2 * \
            np.exp(3./n*(lam+mu)*t) - ((lam+2.*mu)/n)**2) * \
            np.exp(-3./n*(lam+2.*mu)*t)

def inverseSI(si, mu, n, k):

    fix_der = lambda t: derivativeSI(n, mu, k, t)
    fix_snd = lambda t: sndDerivativeSI(n, mu, k, t)
    fix_est = lambda t: estimateSI(n, mu, k, t) - si

#    try:
#        t = newton(fix_est, x0=n*si, fprime=fix_der, fprime2=fix_snd,
#                maxiter=NEWTON_MAX_ITER)
#    except RuntimeError:
#        # try without derivative
    t = newton(fix_est, n*si, maxiter=NEWTON_MAX_ITER)
    return t

def inverseSIHgtOnly(si, lam, n, k):

    x = max(0, 1.-si/(1-2.*k/(n-1)))
    return min(-n/((1-lam)*(3.-(5.*k)/(n-1))) * np.log(x), 6*n)

def evolve(genome, indel_ratio, time, c):
    """ evolves a genome along evolutionary time according to the given rates of
    HGT and indel events.

    The last parameter "c" denotes the next unused gene family ID count and is
    necessary to ensure that only new gene families are inserted in indel
    events.

    The method yields after every time step.
    """

    orig = c-1
    for t in xrange(1, time+1):
        indel_events = random.poisson(indel_ratio)
        new_indel_events = 0
        ps = list()
        for _ in xrange(indel_events):
            c, p, is_new = indel(genome, orig, c)
            new_indel_events += is_new and 1 or 0
            ps.append(p)
        if indel_events:
            LOG.debug(('%s indels occurring at time step %s at position(s) ' + \
                    '%s') %(indel_events, t, ','.join(map(str, sorted(ps)))))
        hgt_events = random.poisson(1-indel_ratio)
        for _ in xrange(hgt_events):
            hgt(genome)
        if hgt_events:
            LOG.debug('%s HGTs occurring at time step %s' %(hgt_events, t))

        yield t, indel_events, new_indel_events


def runExperiment(n, samples, time, indel_ratio, k):
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
        for t, indel_events, new_indel_events in evolve(B, indel_ratio, time,
                n):
            data[s,t-1] = si(A, B, k)
            if t-1:
                new_indels[t-1] = new_indels[t-2]
            new_indels[t-1] += new_indel_events

    return data, new_indels


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('n', type=int, help='number of genes')
    parser.add_argument('-i', '--indelratio', type=float, default=0,
            help='ratio of indels to horzontal gene transfer (HGT) events')
    parser.add_argument('k', nargs='+', type=int)
    parser.add_argument('-s', '--samples', type=int, default=1000,
            help='number of repeats/samples')
    parser.add_argument('-t', '--time', type=int, default=1000,
            help='number of steps of mutations')
    parser.add_argument('-d', '--estimate_distance', action='store_true',
            help='whether or not the evolutionary distance should be estimated')
    args = parser.parse_args()


    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
    plt.figure()
    title = r'$n=%s$' %args.n
    x = np.array(xrange(1, args.time+1))

    data = list()
    for i, k in enumerate(args.k):
        LOG.info(('running simulation with n = %s, k=%s, indelratio = %s, ' + \
                'and %s samples') %(args.n, k, args.indelratio, args.samples))

        data_i, new_indels = runExperiment(args.n, args.samples, args.time,
                args.indelratio, k)

        LOG.info('plotting result')
        ext = len(args.k) > 1 and ' for $k = %s$' %k or ''
        plt.plot(x, [1-np.median(data_i[:, z]) for z in xrange(data_i.shape[1])],
                color = 'C%s' %i, label=r'median SI%s' %ext)
        plt.plot(x, estimateSI(args.n, args.indelratio, k, x), '--',
                color='C%s' %i, label=r'$si%s(t)$%s' %(args.indelratio and '\''
                    or '', ext))

        data.append(data_i)

    if len(args.k) == 1:
        title += r', $k=%s$' %args.k[0]
    if args.indelratio:
        title += r', $\mu=%s$' % args.indelratio
    title += r', %s samples' %args.samples

    plt.title(title)
    plt.xlabel('$t$ (time)')
    plt.ylabel('$d_{\overline{SI}_k}$')
    plt.legend(loc='lower right')
    plt.savefig(stdout, format='pdf')
    plt.close()

    if args.estimate_distance:

        plt.figure()
        for i, k in enumerate(args.k):
            LOG.info(('estimate evolutionary distance from simulated SI ' + \
                    'values for k = %s') %k)
            est_d = np.empty(data[i].shape)
            for s in xrange(data[i].shape[0]):
                if not (s % max(1, data[i].shape[0]/100)):
                    LOG.info('%s%%' %((100*s)/data[i].shape[0]))
                for t in xrange(data[i].shape[1]):
                    iv = args.indelratio and inverseSI or inverseSIHgtOnly
                    est_d[s,t] = iv(1-data[i][s,t], args.indelratio, args.n, k)

            medians = [np.median(est_d[np.isfinite(est_d[:, z]), z]) for z in xrange(est_d.shape[1])]
            upper_q = [np.quantile(est_d[np.isfinite(est_d[:, z]), z], 0.95) for z in xrange(est_d.shape[1])]
            lower_q = [np.quantile(est_d[np.isfinite(est_d[:, z]), z], 0.05) for z in xrange(est_d.shape[1])]

            ext = len(args.k) > 1 and ' for $k = %s$' %k or ''
            plt.plot(x, medians, color='C%s' %i, label = r'median $\hat t$%s' %ext)
            plt.plot(x, upper_q, '-.', color = 'C%s' %i, label = r'$0.95\%%$ quantile%s' %ext)
            plt.plot(x, lower_q, ':', color = 'C%s' %i, label = r'$0.05\%%$ quantile%s' %ext)

        plt.plot(x, x, color='C%s' %len(args.k), label='true distance')
        plt.title(title)
        plt.legend(loc='upper right')
        plt.xlabel('$t$ (time)')
        plt.ylabel('inferred evolutionary distance')
        est_d_flat = est_d.flatten()
        plt.ylim([0, np.min((np.max(est_d_flat[np.isfinite(est_d_flat)]), args.time, 10*args.n))])
        plt.savefig('distance_plot.pdf', format='pdf')

    LOG.info('DONE')
