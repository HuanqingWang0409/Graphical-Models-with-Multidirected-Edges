import math
import numpy as np
import scipy.stats as ss
import scipy
from scipy.stats import ttest_1samp


def all_partition(collection):
    if len(collection) == 1:
        yield [collection]
        return

    first = collection[0]
    for smaller in all_partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[first] + subset] + smaller[n + 1:]
        # put `first` in its own subset
        yield [[first]] + smaller


def moment(data):
    sample_size = len(data[0, :])
    summand = 0
    for n in range(sample_size):
        summand += np.prod(data[:, n])
    return float(summand) / sample_size


def cumulant(data, index):
    # index is a set of the variable indices that we want to compute.
    # index starts at 0
    # Define a function to find all partitions of this index.
    sumofcum = 0
    # For each partition, compute the summand for cumulant.
    for temp in all_partition(index):
        par = sorted(temp)
        l = len(par)
        prod_moment = 1
        for i in par:
            tofindmom = data[i, :]
            prod_moment = prod_moment * moment(tofindmom)
        toadd = pow(-1, l - 1) * math.factorial(l - 1) * prod_moment
        sumofcum += toadd
    return sumofcum


def find_multiedge(data, graph, R, P, X):
    # graph is a bow_acyclic_graph.
    # P, R, and X are lists of vertices. They should be subsets of graph.vertices
    if len(P) == 0 & len(X) == 0:
        yield R
    for v in P:
        pos = [v.order]
        for s in R:
            pos.append(s.order)

        cumu = cumulant(data, pos)
        if R.issubset(v.sib) & cumu != 0:  # May need some tolerance.
            find_multiedge(data, R.append(v), list(set(P).intersection(v.sib)), list(set(X).intersection(v.sib)))
        P = list(set(P).difference({v}))
        X = X.append(v)


