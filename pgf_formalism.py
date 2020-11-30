import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib as m


def pdf_of(degree_list):
    g0 = np.zeros(len(degree_list))
    for i in range(len(g0)):
        g0[i] = degree_list[i] / np.sum(degree_list)
    return g0


def z1_of(g_0):
    z1 = 0
    for k in range(len(g_0)):
        z1 += (k * g_0[k])
    return z1


def z2_of(g_0):
    z2 = 0
    for k in range(len(g_0) - 2):
        z2 += (k + 1) * (k + 2) * g_0[k + 2]
    return z2


def g1_of(g_0):
    g_1 = np.zeros(len(g_0))
    for k in range(len(g_0) - 1):
        g_1[k] = (k + 1) * g_0[k + 1]
    return g_1 / (z1_of(g_0))


def phase_space(g_0, g_1, g=10):
    Psi_sm = np.zeros((10, 100, 100))
    # Initial condition:
    Psi_sm[0][1][1] = 1
    return Psi_sm


def power_law_degree_distrb(maxk):
    p_k = np.empty(maxk)
    p_k[0] = 0
    for k in range(1, maxk):
        p_k[k] = (k ** (-2)) * (math.e ** (-k / 5))
    p_k = p_k / np.sum(p_k)
    return p_k


def generating_function_metrics(gen_func_g0, gen_func_g1):
    # Placeholder function for computing outbreak size and other metrics on generating functions
    G1_func = gen_func_g1
    #
    G1_func[1] = G1_func[1] - 1
    #
    fun = np.poly1d(np.flip(gen_func_g1))
    roots = np.roots(fun)
    u = roots[(roots > 0) & (roots < 1)]

    # Outbreak size, What is going on here with the imaginary numbers
    if len(u) == 0:
        S = 1
    else:
        S = 1 - np.polyval(np.flip(gen_func_g0), u[1])
    print(S)


def gen_functions_with_transmissibility(degree_distrb, T):
    # Given a degree distribution for G0 (or the degree distribution of the entire network).
    # Given transmissibility T
    maxk = 100
    p_k = degree_distrb
    p_LK = np.zeros((maxk, maxk))

    for k in range(0, maxk):
        for l in range(0, k + 1):
            p_LgivenK = p_k[k] * (
                    math.gamma(k + 1) / (math.gamma(l + 1) * math.gamma(k - l + 1)) * T ** (l) * (1 - T) ** (k - l))
            p_LK[k][l] = p_LgivenK
    p_l = np.sum(p_LK, axis=0)
    p_l = p_l / (np.sum(p_l))

    G0_with_T = pdf_of(p_l)
    G1_with_T = g1_of(
        G0_with_T)  # Question for Mariah: Is this correct, or do we need to re-compute p_L USING G1, and then use that?
    return G1_with_T, G0_with_T


def constructMatrixM(g_0, g_1):
    N_0 = len(g_0)
    N_1 = len(g_1)
    M_0 = np.zeros((N_0, N_0))
    M_1 = np.zeros((N_1, N_1))

    M_0[1] = g_0
    newDist = g_1
    M_1[1] = newDist
    for row in range(2, N_1):
        convol = convolve_dists(newDist, g_1)
        M_1[row] = convol
        newDist = M_1[row]

    return (M_0, M_1)


def computeLittlePsi(s, m, prevGenPsi, M):
    s_prime = s - m
    newPsi = prevGenPsi[s_prime].dot(M[:, m])
    return newPsi


def Psi(initProb, num_gens, max_s, max_m, intervention_gen=-1):
    # The number of people that are infective is important for the k values of the matrix
    # The matrix should by and s By m, so the k values should line up with the s values
    allPsi = np.zeros(((num_gens, max_s, max_m)))
    allPsi[0][1][1] = initProb

    original_degree_distrb = power_law_degree_distrb(100)
    initial_T = .8
    g1, g0 = gen_functions_with_transmissibility(original_degree_distrb,
                                                 initial_T)  # this g0 and g1 is for the G(1-(xy+1)T) in terms of the l's
    M_0, M_1 = constructMatrixM(g0, g1)
    for s_g1 in range(max_s):
        for m_g1 in range(max_m):
            allPsi[1][s_g1][m_g1] = computeLittlePsi(s_g1, m_g1, allPsi[0], M_0)

    for g in range(2, num_gens):
        # If g is intervention, re-call g0_l's, g1_l's, M0, M1 etc
        if g == intervention_gen:
            new_T = 0.7
            new_g1, new_g0 = gen_functions_with_transmissibility(original_degree_distrb, new_T)
            new_M = constructMatrixM(new_g0, new_g1)
            M_1 = new_M[1]
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M_1)
    return allPsi


def phaseSpace(num_gens, num_nodes):
    cdict = {
        'red': ((0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green': ((0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue': ((0.0, 0.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }

    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    # need to construct the generating function for psi gen g (prob of having s infected by the end of gen g of which m became infected during gen g
    initProb = 1
    num_nodes = 100
    all_psi_results = Psi(initProb, num_gens, num_nodes, num_nodes)
    fig, ax = plt.subplots()
    inverted_s_m = all_psi_results[5].T  # example for gen 5
    ax.imshow(inverted_s_m[:180][:, :180], cmap=cm)  # gen 5
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
        for i in np.arange(0, m + 1, 1):
            pairs.append([i, m - i])
    elif (m <= len_dist_1 and m > len_dist_2):
        for i in np.arange(m - len_dist_2, m + 1, 1):
            pairs.append([i, m - i])
    elif (m > len_dist_1 and m <= len_dist_2):
        for i in np.arange(0, len_dist_1 + 1, 1):
            pairs.append([i, m - i])
    else:
        for i in np.arange(m - len_dist_2, len_dist_1 + 1, 1):
            pairs.append([i, m - i])
    return pairs


def convolve_dists(X, Y):
    ### X and Y should be vectors of probabilities (For our problem we need to incorporate the G0 and G1 probabilities in here.)
    ### each giving a distribution on a finite subset of the naturals
    len_X = len(X)
    len_Y = len(Y)
    new_dist_len = len_X + len_Y - 1  # Don't think it has to be this long
    new_dist_len = len_X
    new_dist = np.zeros(new_dist_len)
    for m in np.arange(new_dist_len):
        new_prob = 0
        pairs = find_pairs(m, len_X - 1, len_Y - 1)
        for l_pair in pairs:
            new_prob = new_prob + X[l_pair[0]] * Y[l_pair[1]]
        new_dist[m] = new_prob
    return new_dist
