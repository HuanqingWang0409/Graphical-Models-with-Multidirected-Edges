import math
import numpy as np
from scipy.stats import chi2_contingency
from bisect import bisect_right
import scipy.stats


class bow_acyclic_graph:

    def __init__(self, p):
        self.vertices = []

        # initialize vertices in the graph
        for i in range(1, p + 1):
            v = vertex(i)
            self.vertices.append(v)
        self.directed_edges = set()
        self.bidirected_edges = set()

    # set all sibling of a vertex to V\{v}
    def initialize_all_sibs(self):
        length = len(self.vertices)
        for i in range(length):
            for j in range(length):
                if i!=j:
                    (self.vertices[i]).sib.add(self.vertices[j])

    def check_num_sib(self):
        n = 0
        vertices_totest = []
        for v in self.vertices:
            if len(v.sib) > n:
                n = len(v.sib)
            vertices_totest.append(v)
        return (n, vertices_totest)

    def max_num_sib(self):
        return self.check_num_sib()[0]

    def ver_num_sib(self):
        return self.check_num_sib()[1]


class vertex:
    def __init__(self,i):
        self.order=i-1 #Start from 0
        self.par = set()
        self.sib = set()


# v1,v2 are lists of observations.
def indep_var(v1, v2):
    maximum=max(v1+v2)
    minimum=min(v1+v2)
    bin = 8
    interval=float(maximum-minimum)/bin
    freq1=[bisect_right(v1,interval+minimum)]
    freq2=[bisect_right(v2,interval+minimum)]
    for i in range(1,bin):
        temp1=bisect_right(v1,(i+1)*interval+minimum)
        temp2=bisect_right(v2, (i + 1) * interval + minimum)
        freq1.append(temp1-sum(freq1))
        freq2.append(temp2 - sum(freq2))
    obs=[freq1,freq2]
    teststat,pval,dof,expected=chi2_contingency(obs)

    if pval <0.01:
        #We reject the null hypothesis that the variables are independent
        return False
    else:
        return True


def delta_v(C, A, S, D, vertices, v):
    C1 = trans_index(C, vertices)
    A1 = trans_index(A, vertices)
    n = np.shape(D)[0]
    I = np.identity(n)
    first_matrix = np.linalg.pinv(np.matmul((I - D)[np.ix_(C1, A1)], S[np.ix_(A1, C1)]))
    second_matrix = np.matmul((I - D)[np.ix_(C1, A1)], S[np.ix_(A1, [vertices.index(v)])])
    return np.matmul(first_matrix, second_matrix)


def trans_index(C, vertices):
    C_list = []
    for v in C:
        C_list.append(vertices.index(v))
    C_list.sort()
    return C_list


def gamma_v_set(Y, C, A, S, D, vertices, v):
    indx = vertices.index(v)
    C1 = trans_index(C, vertices)
    a = delta_v(C, A, S, D, vertices, v)
    return Y[indx] - np.matmul(np.transpose(delta_v(C, A, S, D, vertices, v)), Y[np.ix_(C1)])


def trans_nonzero_indx(vec):
    nonzero_list = []
    for v in vec:
        if v!=0:
            nonzero_list.append(vec.index(v))
    return nonzero_list

def trans_indx_to_vertex(vec,vertices):
    set_to_return = set()
    for v in vec:
        set_to_return.add(vertices[v])
    return set_to_return

def an_D(C,D,vertices):
    set_to_return = set()
    n = np.shape(D)[0]
    C1 = trans_index(list(C),vertices)
    inv_matrix = np.linalg.pinv(np.identity(n)-D)
    inv_matrix[np.abs(inv_matrix)<0.01]=0
    for c in C1:
        list_to_test = trans_nonzero_indx((inv_matrix).tolist()[c])
        set_to_return = set_to_return.union(trans_indx_to_vertex(list_to_test,vertices))
    return set_to_return

#A is to be replaced, B is 1xN list, C is 1xN list of indices
def replace_sub_mat(A,B,C,v):
    i=0
    for c in C:
        A[v][c]=B[i]
        i+=1



