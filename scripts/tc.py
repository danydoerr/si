#!/usr/bin/env python3
from sys import argv

from dendropy import Tree, TaxonNamespace
from dendropy.calculate import treecompare

tns = TaxonNamespace()
T0 = Tree.get(file=open(argv[1]), schema='Newick',
        taxon_namespace=tns)
T1 = Tree.get(file=open(argv[2]), schema='Newick',
        taxon_namespace=tns)
#fp, fn = treecompare.false_positives_and_negatives(T0, T1)
print(treecompare.symmetric_difference(T0, T1))
