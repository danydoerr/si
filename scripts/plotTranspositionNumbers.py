#!/usr/bin/env python2

from sys import argv, stdout, stderr, exit
from matplotlib import pylab as plt
import numpy as np


if __name__ == '__main__':

    data = np.loadtxt(argv[1])

    xs = sorted(set(data[:, 1]))
    ys = [data[data[:, 1] == x, 2] for x in xs]

    bp = plt.boxplot(ys, labels=map(int, xs), positions = xs, notch=0, sym='+',
            vert=1, whis=1.5, widths=20)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+', markersize=4)

    plt.ylabel('total number of transpositions')
    plt.xlabel('PAM')
    plt.xlim([-10, 850])

    plt.show()

