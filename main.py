import matplotlib.pyplot as plt
import networkx as nx
import random
import numpy as np
import itertools as it
import scipy
from scipy import stats
from matplotlib import rc
import math



def infections_caused_matrix(P_k, beta, x=1):
    dist = np.zeros((len(P_k), len(P_k)))
    for k in range(len(P_k)):
        for l in range(k):
            dist[l][k] = P_k[k]*p_l_infected(k, l, beta)*(x**l)
    return dist

# def p_l_infected(k, l, beta):
#     return scipy.special.binom(k, l)*(beta**l)*((1-beta)**(k-l))

def little_test():
    degree_data = np.random.poisson(10, 50)
    P_k = pdf_of(degree_data)
    plt.plot(P_k)
    plt.show()
    g_0_matrix = infections_caused_matrix(P_k, .15)
    avg_degree = z1_of(P_k)
    z2 = z1_of(g1_of(g1_of(P_k)))   # This might be wrong
    beta = .3
    r_0 = beta*z2/avg_degree
    g_1 = g1_of(P_k)
    g_1_matrix = infections_caused_matrix(g_1, .15)
    # myplot = plt.plot(g_0_matrix)
    fig, ax = plt.subplots()
    fig1 = ax.imshow(g_0_matrix, cmap='plasma')
    plt.show()
    fig, ax = plt.subplots()
    fig2 = ax.imshow(g_1_matrix, cmap='plasma')
    plt.show()
    Psi_sm = phase_space(P_k, g_1)
    fig, ax = plt.subplots()
    fig3 = ax.imshow(Psi_sm[0], cmap='plasma')
    plt.show()
    return g_0_matrix

def pdf_of(degree_list):
    g0 = np.zeros(len(degree_list))
    for i in range(len(g0)):
        g0[i] = degree_list[i]/ np.sum(degree_list)
    return g0

def z1_of(g_0):
    z1=0
    for k in range(len(g_0)):
        z1+=(k * g_0[k])
    return z1

def z2_of(g_0):
    z2=0
    for k in range(len(g_0)-2):
        z2+=(k+1)*(k+2)*g_0[k+2]
    return z2


def g1_of(g_0):
    g_1 = np.zeros(len(g_0))
    for k in range(len(g_0)-1):
        g_1[k] = (k+1)*g_0[k+1]
    return g_1/(z1_of(g_0))

def phase_space(g_0, g_1, g=10):
    Psi_sm = np.zeros((10, 100,100))
    # Initial condition:
    Psi_sm[0][1][1] = 1
    return Psi_sm





def formalism():
    # Given a degree distribution for G0 (or the degree distribution of the entire network).
    T = 0.5 # Transmissability
    r0 = 3 # mean
    a =  0.2 # dispersion/clustering coefficient

    maxk = 100
    p_k = np.empty(maxk)
    degreeDist = np.empty(maxk)
    p_LK = []
    for k in range(maxk):
        p_k[k] = math.gamma(r0 + k)/(math.factorial(k)* math.gamma(r0))*(a/(r0+a))**(a) * (a/(r0+a))**(k) # make vector
        p_LgivenK = []
        for l in range(k):
            p_LgivenK.append(math.gamma(k + 1) / (math.gamma(l + 1) * math.gamma(k - l + 1)) * T**(l) * (1 - T)**(k - l))
        p_LK.append(p_LgivenK)
        #p_l[k] = np.sum(p_LK[k]) # need to sum these the other way


    p_LK =  [elem[::-1] for elem in p_LK] # Need to normalize this by the columns
    b = np.zeros([len(p_LK), len(max(p_LK, key=lambda x: len(x)))])
    for i, j in enumerate(p_LK):
        b[i][0:len(j)] = j
    p_l = b.sum(axis = 0)
    for r in range(maxk-1):
        degreeDist[r] = p_k[r] * p_l[r]
    start_G0 = pdf_of(degreeDist)

    start_G1 = g1_of(start_G0)

    G1_func = start_G1

    G1_func[1] = G1_func[1] - 1

    fun = np.poly1d(np.flip(start_G1))
    roots = np.roots(fun)
    u = roots[(roots > 0) & (roots < 1)]

    # Outbreak size, What is going on here with the imaginary numbers
    if len(u) == 0:
        S = 1
    else:
        S = 1 - np.polyval(np.flip(start_G0), u[1])
    #print(S)
    return start_G1, start_G0

def constructMatrixM():
    g1, g0 = formalism()
    #N_0 = len(g0)
    N_1 = len(g1)
    # Need to make matrix that does all the powers of G_g-1 (collect like terms)
    #M_0 = np.empty([N_0, N_0])
    M_1 = np.zeros((N_1, N_1))
    # For loop to feed new dist into the convolve again

    M_0 = g0
    newDist = g1
    M_1[1] = newDist
    for row in range(2,N_1):
        M_1[row] = convolve_dists(newDist, g1)
        newDist = M_1[row]

    return M_1

def computeLittlePsi(s, m, prevGenPsi, M):
    s_prime = s - m
    newPsi = prevGenPsi[s_prime].dot(M[:,m])
    return newPsi


def layeredPsi(g0, g1, initProb, gen, s_count, m_count, M_0, M_1):
    # The number of people that are infective is important for the k values of the matrix
    # The matrix should by and s By m, so the k values should line up with the s values
    allPsi = np.zeros(((gen, s_count, m_count)))
    #onePsiMat = np.zeros((s_count, m_count))
    allPsi[0][1][1] = initProb
    for s in range(s_count):
        for m in range(m_count):
            #probMat[s, m] = initProb*g0*(g1)**(gen-1) INCORRECT
            allPsi[1][s][m] = computeLittlePsi(s, m, allPsi[0], M_0)
            #probMat[s,m] = initProb*
        #Figure out what is going on with the sequences here, should we be accounting for k value somewhere?
                # This is where the implementation for the convolution needs to come into play

    for g in range(gen):
        for s in range(s_count):
            for m in range(m_count):
                # probMat[s, m] = initProb*g0*(g1)**(gen-1) INCORRECT
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g-1], M_1)
    return allPsi

def phaseSpace(gen, s, m):
    # need to construct the generating function for psi gen g (prob of having s infected by the end of gen g of which m became infected during gen g
    g0, g1 = formalism()
    initProb = 1
    mat = psiRecursion(g0, g1, initProb, gen, s, m)
    return mat

# How to structure this code.
# The phase space formalism needs
        # 1) Need the G_{g-1} formalism (This includes the double sum formula) **DONE
        # 2) Need to put together the matrix that Andrea is drawing,
            # this matrix will have rows of PGFs, each row is a exponentiated G_{g-1}
        # 3) Convolve the rows of this matrix to get one value associated with [G_{g-1}]^m'
        # 4) The value from the convolution will then lead to the Psi matrix

# This is the convolution code
from scipy.stats import binom


def find_pairs(m, len_dist_1, len_dist_2):
    ### must have all three args nat num valued
    ### must be that m <= len_dist_1 + len_dist_2
    pairs = []
    if (m <= len_dist_1 and m <= len_dist_2):
        for i in np.arange(0,m+1,1):
            pairs.append([i,m-i])
    elif (m <= len_dist_1 and m > len_dist_2):
        for i in np.arange(m-len_dist_2,m+1,1):
            pairs.append([i,m-i])
    elif (m > len_dist_1 and m <= len_dist_2):
        for i in np.arange(0, len_dist_1+1,1):
            pairs.append([i,m-i])
    else:
        for i in np.arange(m-len_dist_2 ,len_dist_1+1,1):
            pairs.append([i,m-i])
    return pairs

def convolve_dists(X,Y):
    ### X and Y should be vectors of probabilities (For our problem we need to incorporate the G0 and G1 probabilities in here.)
    ### each giving a distribution on a finite subset of the naturals
    len_X = len(X)
    len_Y = len(Y)
    new_dist_len = len_X + len_Y - 1
    new_dist = np.zeros(new_dist_len)
    for m in np.arange(new_dist_len):
        new_prob = 0
        pairs = find_pairs(m, len_X-1, len_Y-1)
        for l_pair in pairs:
            new_prob = new_prob + X[l_pair[0]]*Y[l_pair[1]]
        new_dist[m] = new_prob
    return new_dist

# def gen_ic_probs(sens, vec):
#     ### sens is test sensitivity should be in [0,1]
#     ### vec is a vector of pairs (n,p) where n is natural num p is a prob
#     ### n represents number of incoming students from a county
#     ### p is estimated prevalence in county/prob of having COVID if from county
#     ic_dist = binom.pmf(np.arange(vec[0][0]+1), vec[0][0], vec[0][1]*(1-sens))
#     for v in vec[1:]:
#         ic_dist = convolve_dists(ic_dist, binom.pmf(np.arange(v[0]+1), v[0], v[1]*(1-sens)))
#     return ic_dist / sum(ic_dist)




if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')
    #little_test()
    #formalism()
    probMat = phaseSpace(2, 4, 2)