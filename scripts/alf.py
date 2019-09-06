#!/usr/bin/env python

from collections import deque
from Bio import SeqIO
import re

FASTA_HEADER_PAT = re.compile('^(G\d+_SE\d+), sequence type: (.*), locus: (-?\d+)$')
TOKEN_PAT = re.compile('\[|\]|\d+')

DIRECTION_CRICK_STRAND = '+'
DIRECTION_WATSON_STRAND = '-'



def readHomologies(vp_data):
    res = None

    stack = deque()
    for line in vp_data:
        if line.startswith('#'):
            continue
        start = line.find(':=')
        for token in TOKEN_PAT.findall(line[start+2:]):
            # get current list from stack
            cur_elem = None
            if stack:
                cur_elem = stack.pop()

            # parse token
            if token == '[':
                new_l = list()
                if cur_elem != None:
                    cur_elem.append(new_l)
                else:
                    # add root node twice, so it will be on top of stack after
                    # all elements have been processed
                    cur_elem = new_l

                stack.append(cur_elem)
                cur_elem = new_l

            elif token == ']':
                cur_elem = None
            elif token.isdigit():
                cur_elem.append(int(token))

            if cur_elem != None:
                stack.append(cur_elem)
    return stack.pop()

def matrix2families(mtrx):

    res = [[0] * len(mtrx[i][0]) for i in xrange(len(mtrx))]

    c = 1
    for i in xrange(len(mtrx)):
        for p in xrange(len(mtrx[i][0])):
            for j in xrange(len(mtrx)):
                if mtrx[i][j][p] and res[j][mtrx[i][j][p][0]-1]:
                    res[i][p] = res[j][mtrx[i][j][p][0]-1]
                    break
            if not res[i][p]:
                res[i][p] = c
                c += 1
    return res

def readGenePositions(fastadata):
    res = list()

    i = 0
    for record in SeqIO.parse(fastadata, 'fasta'):
        m = FASTA_HEADER_PAT.match(record.description)
        if not m:
            print >> stderr, 'Unable to parse FASTA sequence header "%s". Exiting' %(record.description)
            exit(1)
        sid, _, locus = m.groups()

        locus = int(locus)
        orient = DIRECTION_CRICK_STRAND
        if locus < 0:
            orient = DIRECTION_WATSON_STRAND
            locus = -locus
        res.append((locus, orient, i, sid))
        i += 1

    res.sort()
    return [(DIRECTION_CRICK_STRAND, -1, 'TELOMERE_START')] + map(lambda x: x[1:4], res)

