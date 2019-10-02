#!/usr/bin/env python2

from sys import argv, stdout, stderr, exit
from matplotlib import pylab as plt
import numpy as np


if __name__ == '__main__':

    data = np.loadtxt(argv[1])

    xs = sorted(set(data[:, 1]))
    ys = [np.sum(data[data[:, 1] == x, 2:4], axis=1) for x in xs]

    bp = plt.boxplot(ys, labels=map(int, xs), positions = xs, notch=0, sym='+',
            vert=1, whis=1.5, widths=20)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+', markersize=4)

    plt.ylabel('total number of insertions+deletions')
    plt.xlabel('PAM')
    plt.xlim([-10, 550])

    plt.show()

