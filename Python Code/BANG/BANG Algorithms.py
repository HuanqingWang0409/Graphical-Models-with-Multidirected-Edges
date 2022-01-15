import graph
import math
import numpy as np
import scipy.stats
import itertools
import random
import copy

def algorithm1(S,Y):

    #line 2
    p = S.shape[0]
    G = graph.bow_acyclic_graph(p)
    G.initialize_all_sibs()
    D = np.zeros([p,p])
    l = 1
    #line 4-10
    while G.max_num_sib() >= l:
        all_vertices = copy.deepcopy(G.vertices)
        i=0
        for v1 in all_vertices:
            v1.sib = G.vertices[i].sib.copy()
            i+=1
        D1=D.copy()
        for v in G.vertices:
            # v = G.vertices[all_vertices.index(v1)]
            v1 = all_vertices[G.vertices.index(v)]
            algorithm2(v, v.sib, Y, D, G.vertices)
            algorithm3(v, v.par, v.sib, D,D, S, Y, l, G.vertices)
        if (D1==D).all():
            l+=1
        else:
            l=1
    #line 11
    R = D
    algorithm4(G,D,S,Y,l)
    #line 12
    print(D)
    return G,D


def algorithm2(v,sibv, Y,D,vertices):
    #line 2
    gamma = Y-np.matmul(D,Y)
    sib_v = sibv.copy()
    for u in sib_v:
        gamma_u = np.array(gamma)[vertices.index(u)]
        gamma_v = np.array(gamma)[vertices.index(v)]
        if graph.indep_var(gamma_u,gamma_v):
            v.sib.discard(u)
            u.sib.discard(v)


def algorithm3(v, pav, sibv, D1,D, S, Y, l, vertices):
    C_1=set()
    if (len(sibv)>=l):
        for C in findsubsets(sibv,l):
            i = True
            for c in C:
                gamma_c = np.array(Y - np.matmul(D1, Y))[vertices.index(c)]
                r = graph.an_D(C.union(v.par), D1, vertices)
                gamma_v = np.array(graph.gamma_v_set(Y, C.union(v.par), graph.an_D(C.union(v.par), D1, vertices), S, D1, vertices, v))[0]
                indep_or_not = graph.indep_var(gamma_c,gamma_v)
                if (indep_or_not == False):
                    i = False
                    break
            if i==True:
                C_1 = C_1.union(C)
        v.par = v.par.union(C_1)
        r = graph.trans_index(v.par,vertices)
        a = graph.trans_index(graph.an_D(v.par, D1, vertices),vertices)
        delta_v = graph.delta_v(v.par, graph.an_D(v.par, D1, vertices), S, D1, vertices, v)
        delta_v_transposed = np.array(np.transpose(delta_v))[0]
        par_indx = graph.trans_index(v.par,vertices)
        if len(delta_v)>0:
            graph.replace_sub_mat(D, delta_v_transposed, par_indx, vertices.index(v))
        v.sib = v.sib.difference(v.par)

        for s in v.par:
            s.sib = s.sib.difference({v})


def algorithm4(G,D,S,Y,l):
    top_order = (topo_order(D))
    top_order.reverse()
    for i in range(len(top_order)):
        v = G.vertices[top_order[i]]
        par_v = v.par.copy()
        for s in par_v:
            k = par_v.difference({s})
            gamma_s = np.array(Y-np.matmul(D,Y))[G.vertices.index(s)]
            gamma_v = np.array(graph.gamma_v_set(Y, par_v.difference({s}), graph.an_D(par_v.difference({s}), D, G.vertices), S, D, G.vertices,v))[0]
            indep_or_not = graph.indep_var(gamma_s,gamma_v)
            if indep_or_not:
                a = G.vertices.index(v)
                b = G.vertices.index(s)
                v.par.discard(s)
                D[G.vertices.index(v)][G.vertices.index(s)]=0



def findsubsets(S,m):
    if (m<= len(S)):
        return list(map(set, itertools.combinations(S,m)))

def last_non_zero_indx(vec):
    i=0
    for x in reversed(range(len(vec))):
        if vec[x]!=0:
            break
        i+=1
    return i

def topo_order(D):
    I = np.identity(D.shape[0])
    D_inv = np.linalg.pinv(I - D)
    D_inv[np.abs(D_inv)<0.01]=0
    list_to_test = []
    for i in range(D.shape[0]):
        test = np.array(D[i])
        cord = (i, last_non_zero_indx(np.array(D_inv[i])))
        list_to_test.append(cord)
    list_to_test.sort(key = lambda tup:tup[1],reverse=True)
    list_to_return = []
    for i in range(D.shape[0]):
        list_to_return.append((list_to_test[i])[0])
    return list_to_return

def random_list_unif(num, start, stop):
    result = []
    for i in range(num):
        result.append(round(random.uniform(start,stop),2))
    return result




