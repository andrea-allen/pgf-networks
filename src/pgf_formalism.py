import numpy as np
import math
import scipy.special


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


def gen_functions_with_transmissibility(degree_distrb, T):
    # Given a degree distribution for G0 (or the degree distribution of the entire network).
    # Given transmissibility T
    # TODO this is what will change depending on the intervention
    # TODO don't modify this one but make new function for each vaccination strategy
    maxk = len(degree_distrb)
    p_k = degree_distrb

    # This is the matrix of k x l resulting p l given k in each cell
    p_LK = np.zeros((maxk, maxk))

    # Generates pgf in variable l as probabilities of infection of l neighbors given the original degree distribution and transmission prob T
    for k in range(0, maxk):
        # somewhere here have another matrix that's l by j
        # construct that for each p_j_given_l (whatever the order is)
        for l in range(0, k + 1):
            try:
                # this will be some operation of this with the vector for the particular l
                # sum over all j for this particular l
                p_LgivenK = p_k[k] * (
                        math.gamma(k + 1) / (math.gamma(l + 1) * math.gamma(k - l + 1)) * T ** (l) * (1 - T) ** (k - l))
                p_LK[k][l] = p_LgivenK
            except OverflowError:
                p_LK[k][l] = 0
    p_l = np.sum(p_LK, axis=0)
    p_l = p_l / (np.sum(p_l))

    G0_with_T = pdf_of(p_l)
    G1_with_T = g1_of(G0_with_T)
    return G1_with_T, G0_with_T


def random_vacc_distribution(degree_distrb, T, V):
    maxk = len(degree_distrb)

    # Important! Only feed this function the oringinal degree distribution
    # it will be then transformed to g1, assuming this is happening
    # at a generation later than 0
    # TODO Make sure this gets normalized?
    q_k = g1_of(degree_distrb)

    # we ultimately want:
    # vector of powerks of k
    # matrix of j and k (j behaves like the oriignal k in the paper)
    # matrix of j and l


    P_lk = np.zeros((maxk, maxk))

    for k in range(0, maxk):
        P_jk = np.zeros(maxk)
        for j in range(0, k + 1):
            try:
                p_j_given_k = q_k[k] * (math.gamma(k + 1) / (math.gamma(j + 1) * math.gamma(k - j + 1))
                                        * ((1 - V) ** (j)) * (V ** (k - j)))
                P_jk[j] = p_j_given_k
            except OverflowError:
                P_jk[j] = 0

            for l in range(0, j + 1):
                try:
                    p_l_given_j = (math.gamma(j + 1) / (math.gamma(l + 1) * math.gamma(j - l + 1))) * ((1 - T) ** (j - l)) * (T ** (l))
                    P_lk[k][l] = P_jk[j] * p_l_given_j
                except OverflowError:
                    P_lk[k][l] = 0

    P_l_distrb = np.sum(P_lk, axis=0)
    P_l_distrb = P_l_distrb / (np.sum(P_l_distrb))

    # This return value is in the form of an excess degree distribution
    # Question for Weds: do we need to transform again into excess, or did
    #the original G1 take care of that?
    return P_l_distrb



def constructMatrixM(g_0, g_1):
    # Constructs the matrix of pgfs for G0_with transmission convolved to every mth power
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

    M_1[0][0] = 1
    return (M_0, M_1)


def computeLittlePsi(s, m, prevGenPsi, M):
    s_prime = s - m
    newPsi = prevGenPsi[s_prime].dot(M[:, m])
    return newPsi


# COMPUTATION STARTS HERE
def Psi(degree_distrb, initProb, num_gens, max_s, max_m, initial_T, intervention_gen=-1, intervention_T=0.5, prop_vacc=0.5):
    # 3-d matrix with one matrix per generation of Psi_g
    allPsi = np.zeros(((num_gens, max_s, max_m)))
    allPsi[0][1][1] = initProb

    # Assign initial degree distribution here
    original_degree_distrb = degree_distrb
    g1, g0 = gen_functions_with_transmissibility(original_degree_distrb,
                                                 initial_T)  # this g0 and g1 is for the G(1-(xy+1)T) in terms of the l's
    M_0, M_1 = constructMatrixM(g0, g1)
    for s_g1 in range(max_s):
        for m_g1 in range(max_m):
            allPsi[1][s_g1][m_g1] = computeLittlePsi(s_g1, m_g1, allPsi[0], M_0)

    for g in range(2, num_gens):
        print('working on gen ' + str(g))
        # If g is intervention, re-call g0_l's, g1_l's, M0, M1 etc
        if g == intervention_gen:
            new_T = intervention_T
            # THIS IS WHERE WE:D MODIFY, make a version of this function that's the intervention
            # new_g1, new_g0 = gen_functions_with_transmissibility(original_degree_distrb, new_T)
            #when using the random one, feed it the old T
            new_g1 = random_vacc_distribution(degree_distrb, initial_T, prop_vacc)
            new_M = constructMatrixM(g0, new_g1) # TODO Is this correct? to have the oriignal g0 here?
            M_1 = new_M[1]
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M_1)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g
    return allPsi


# START HERE WHEN CALLING
def phaseSpace(num_gens, num_nodes, degree_distribution, transmissibility,
               save_results=True, gens_to_save=None,
               file_fmt_to_save='phase_space/generation_{0}', intervention_gen=-1, intervention_trans=None, vacc_pop=.5,
               do_non_interv = True, do_interv = True):
    # the generating function for psi gen g (prob of having s infected by the end of gen g of which m became infected during gen g
    initProb = 1

    if do_non_interv:
        all_psi_results = Psi(degree_distribution, initProb, num_gens, num_nodes, num_nodes, transmissibility)

        if save_results:
            try:
                for gen in gens_to_save:
                    if gen < num_gens:
                        np.savetxt(file_fmt_to_save.format(gen) + '.txt', all_psi_results[gen], delimiter=',')
            except Exception:
                print('Must provide gens_to_save in arguments as list')
    else:
        all_psi_results = np.zeros((2,2))

    if do_interv:
        all_psi_results_with_intervention = Psi(degree_distribution, initProb, num_gens, num_nodes, num_nodes,
                                                transmissibility, intervention_gen, intervention_trans, vacc_pop)
        if save_results:
            try:
                for gen in gens_to_save:
                    np.savetxt(file_fmt_to_save.format(gen) + '_intv.txt', all_psi_results_with_intervention[gen],
                               delimiter=',')
            except Exception:
                print('Must provide gens_to_save in arguments as list')
    else:
        all_psi_results_with_intervention = np.zeros((2, 2))
    return all_psi_results, all_psi_results_with_intervention


# Convolution code below:
# Used in methods generating the phase space pgf's
# The phase space formalism needs
# 1) Need the G_{g-1} formalism (This includes the double sum formula) **DONE
# 2) Need to put together the matrix that Andrea is drawing,
# this matrix will have rows of PGFs, each row is a exponentiated G_{g-1}
# 3) Convolve the rows of this matrix to get one value associated with [G_{g-1}]^m'
# 4) The value from the convolution will then lead to the Psi matrix

# This is the convolution code


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
