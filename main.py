import matplotlib.pyplot as plt
import numpy as np
import itertools as it
import scipy
from matplotlib import rc
import math
import matplotlib as m


import SIR_sims

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
    z2 = 0                       #Is this wrong?
    for k in range(len(g_0)-2):
        z2 += (k+1)*(k+2)*g_0[k+2]
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
    T = 0.8 # Transmissability
    r0 = 3 # mean
    a =  0.2 # dispersion/clustering coefficient

    maxk = 100
    p_k = np.empty(maxk)
    degreeDist = np.empty(maxk)
    p_LK = []
    p_k[0] = 0
    for k in range(1, maxk):
        p_k[k] = (k ** (-2)) * (math.e ** (-k / 5))
        # p_k[k] = math.gamma(r0 + k)/(math.factorial(k)* math.gamma(r0))*(a/(r0+a))**(a) * (a/(r0+a))**(k) # make vector
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

    # G1_func = start_G1
    #
    # G1_func[1] = G1_func[1] - 1
    #
    # fun = np.poly1d(np.flip(start_G1))
    # roots = np.roots(fun)
    # u = roots[(roots > 0) & (roots < 1)]

    # Outbreak size, What is going on here with the imaginary numbers
    # if len(u) == 0:
    #     S = 1
    # else:
    #     S = 1 - np.polyval(np.flip(start_G0), u[1])
    #print(S)
    return start_G1, start_G0

def constructMatrixM(g_0, g_1):
    N_0 = len(g_0)
    N_1 = len(g_1)
    # Need to make matrix that does all the powers of G_g-1 (collect like terms)
    M_0 = np.zeros((N_0, N_0))
    M_1 = np.zeros((N_1, N_1))
    # For loop to feed new dist into the convolve again

    M_0[1] = g_0
    newDist = g_1
    M_1[1] = newDist
    for row in range(2,N_1):
        convol = convolve_dists(newDist, g_1)
        M_1[row] = convol
        newDist = M_1[row]

    return (M_0, M_1)

def computeLittlePsi(s, m, prevGenPsi, M):
    s_prime = s - m
    newPsi = prevGenPsi[s_prime].dot(M[:,m])
    return newPsi


def layeredPsi(initProb, num_gens, s_count, m_count, M_0, M_1):
    # The number of people that are infective is important for the k values of the matrix
    # The matrix should by and s By m, so the k values should line up with the s values
    allPsi = np.zeros(((num_gens, s_count, m_count)))
    #onePsiMat = np.zeros((s_count, m_count))
    allPsi[0][1][1] = initProb
    for s in range(s_count):
        for m in range(m_count):
            #probMat[s, m] = initProb*g0*(g1)**(num_gens-1) INCORRECT
            allPsi[1][s][m] = computeLittlePsi(s, m, allPsi[0], M_0)
            #probMat[s,m] = initProb*
        #Figure out what is going on with the sequences here, should we be accounting for k value somewhere?
                # This is where the implementation for the convolution needs to come into play

    for g in range(2, num_gens):
        for s in range(s_count):
            for m in range(m_count):
                # probMat[s, m] = initProb*g0*(g1)**(gen-1) INCORRECT
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g-1], M_1)
    return allPsi

def phaseSpace(num_gens):
    cdict = {
        'red': ((0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green': ((0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue': ((0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }

    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    # need to construct the generating function for psi gen g (prob of having s infected by the end of gen g of which m became infected during gen g
    g0, g1 = formalism()
    initProb = 1
    M = constructMatrixM(g0, g1)
    all_psi_results = layeredPsi(initProb, num_gens, len(g0), len(g0), M[0], M[1])
    fig, ax = plt.subplots()
    inverted_s_m = all_psi_results[5].T
    ax.imshow(inverted_s_m[:30][:,:30], cmap=cm) #gen 5
    ax.invert_yaxis()
    plt.title('Phase Space at Generation 5 of Power Law Network')
    plt.ylabel('$m$')
    plt.xlabel('$s$')
    plt.savefig('draft_phase_space.png')
    plt.show()
    return all_psi_results

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
    new_dist_len = len_X + len_Y - 1 # Don't think it has to be this long
    new_dist_len = len_X
    new_dist = np.zeros(new_dist_len)
    for m in np.arange(new_dist_len):
        new_prob = 0
        pairs = find_pairs(m, len_X-1, len_Y-1)
        for l_pair in pairs:
            new_prob = new_prob + X[l_pair[0]]*Y[l_pair[1]]
        new_dist[m] = new_prob
    return new_dist


if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')
    #little_test()
    #formalism()
    probMat = phaseSpace(20)
    SIR_sims.run()
