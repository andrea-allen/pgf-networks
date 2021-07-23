import math

import numpy as np

from src import gen_extinct_prob


def compute_extinct_prob_all(deg_dist=None, T=1.0, n_gens=20, renorm=True, fft = True, custom_g0=None, custom_g1=None):
    if deg_dist is not None:
        psi = Psi(deg_dist, initProb=1, num_gens=n_gens, max_s=len(deg_dist), max_m=len(deg_dist), initial_T=T)
    else:
        psi = Psi(deg_dist, initProb=1, num_gens=n_gens, initial_T=T, custom_g0=custom_g0, custom_g1=custom_g1,
                  max_m=len(custom_g1), max_s=len(custom_g1))
# =============================================================================
#     for g in range(n_gens):
#         psi[g][:,0] = np.zeros(psi.shape[0])
#         psi[g] = psi[g]/np.sum(psi[g])
# =============================================================================
    if renorm:
        for g in range(n_gens):
            for s in range(psi.shape[1]):
                if np.sum(psi[g][s][:]) > 0:
                    psi[g][s,0] = 0
                    psi[g][s,:] = psi[g][s,:] / np.sum(psi[g][s,:])
    if deg_dist is None:
        deg_dist = custom_g1 # this is just for computing the extinction prob
        extct_array = gen_extinct_prob.gen_ext_prob_array(psi, deg_dist, T, fft=fft, custom=True)
    else:
        extct_array = gen_extinct_prob.gen_ext_prob_array(psi, deg_dist, T, fft=fft)
    return [extct_array, psi]

def expected_num_infected(deg_dist, T):
    psi = Psi(deg_dist, initProb=1, num_gens=50, max_s=400, max_m=400, initial_T=T)
    expected_cum_array = np.zeros(50)
    for gen in range(50):
        psi[gen][:,0] = np.zeros(400)
        psi[gen] = psi[gen]/np.sum(psi[gen])
        inverted_s_m = psi[gen]
        ps_g_analytical = np.sum(inverted_s_m, axis=0)
        ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)
        expected_cum_array[gen] = np.sum(ps_g_analytical * np.arange(len(ps_g_analytical)))
    np.savetxt(f'./../data/expected_cum_{T}_m.txt', expected_cum_array)
    print(expected_cum_array[:10])

# CALL FROM OUTSIDE CLASS CALLS HERE
def compute_phase_space(num_gens, num_nodes, degree_distribution, transmissibility,
                        save_results=True, gens_to_save=None,
                        file_fmt_to_save='phase_space/generation_{0}', intervention_gen=-1, intervention_trans=None,
                        vacc_pop=.5, rollout_dict=None,
                        do_non_interv=True, do_interv=True, intervention_type="none"):
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
        all_psi_results = np.zeros((2, 2))

    if do_interv:
        all_psi_results_with_intervention = Psi(degree_distribution, initProb, num_gens, num_nodes, num_nodes,
                                                transmissibility, intervention_gen, intervention_trans, vacc_pop, rollout_dict, intervention_type)
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

    # Important! Only feed this function the original degree distribution
    # it will be then transformed to g1, assuming this is happening
    # at a generation later than 0
    # TODO Make sure this gets normalized?
    q_k = g1_of(degree_distrb)

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
                    p_l_given_j = (math.gamma(j + 1) / (math.gamma(l + 1) * math.gamma(j - l + 1))) * (
                            (1 - T) ** (j - l)) * (T ** (l))
                    P_lk[k][l] = P_jk[j] * p_l_given_j
                except OverflowError:
                    P_lk[k][l] = 0

    P_l_distrb = np.sum(P_lk, axis=0)
    P_l_distrb = P_l_distrb / (np.sum(P_l_distrb))
    return P_l_distrb


def critical_degree_calc(prop, degree_d):
    k_crit = len(degree_d) #default for the critical k value
    temp_prop = prop
    while temp_prop > 0:
        temp_prop = temp_prop - degree_d[k_crit-1]
        k_crit -= 1

    return k_crit

def targeted_vacc_distribution(degree_distrb, T, V):
    maxk = len(degree_distrb)

    # Important! Only feed this function the original degree distribution
    # it will be then transformed to g1, assuming this is happening
    # at a generation later than 0

    q_k = g1_of(degree_distrb)
    crit_value = critical_degree_calc(V, degree_distrb)

    num_H = 0
    for c in range(crit_value, maxk-1):
        num_H += c*q_k[c]

    denom_H = 0
    for f in range(maxk-1):
        denom_H += f*q_k[f]

    H = num_H / denom_H

    P_lkTarget = np.zeros((maxk, maxk))

    # need to have an if statement in the second for loop that adjusts the transmissibility
    for k in range(0, maxk):
        P_jkTarget = np.zeros(maxk)
        for j in range(0, k + 1):

            try:
                p_j_given_k_tar = q_k[k] * (math.gamma(k + 1) / (math.gamma(j + 1) * math.gamma(k - j + 1))
                                        * ((1 - H) ** (j)) * (H ** (k - j)))
                P_jkTarget[j] = p_j_given_k_tar
            except OverflowError:
                P_jkTarget[j] = 0

            for l in range(0, j + 1):

                #Determine what k is to change T_k accordingly
                T_k = T
                if k >= crit_value:
                    T_k = 0

                try:
                    p_l_given_j_tar = (math.gamma(j + 1) / (math.gamma(l + 1) * math.gamma(j - l + 1))) * (
                            (1 - T_k) ** (j - l)) * (T_k ** (l))
                    P_lkTarget[k][l] = P_jkTarget[j] * p_l_given_j_tar
                except OverflowError:
                    P_lkTarget[k][l] = 0

    P_l_tar_distrb = np.sum(P_lkTarget, axis=0)
    P_l_tar_distrb = P_l_tar_distrb / (np.sum(P_l_tar_distrb))
    return P_l_tar_distrb


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
        second_convol = convolve_with_ifft(newDist, g_1)
        M_1[row] = second_convol
        newDist = M_1[row]

    M_1[0][0] = 1 #nice
    return (M_0, M_1)


def computeLittlePsi(s, m, prevGenPsi, M):
    s_prime = s - m
    newPsi = prevGenPsi[s_prime, :].dot(M[:, m])
    return newPsi

def offspring_dists(r0, k, p0, length):
    a = 1/k
    g1 = np.zeros(length)
    g0 = np.zeros(length)
    for i in range(len(g1)):
        try:
            # fact = np.sqrt(2 * math.pi * i) * (i / np.exp(1)) ** i
            g1[i] = (math.gamma(i + k) / (math.factorial(i) * math.gamma(k))) * ((a * r0) / (1 + a * r0)) ** (i) * (
                        1 / (1 + a * r0)) ** (k)
        except OverflowError:
            g1[i] = 0
    g1 = g1 / np.sum(g1)
    g0 = compute_g0_from_offspring(g1, p0)
    g0 = g0[:-1]
    return g0, g1

def compute_g0_from_offspring(g1, p0):
    g0 = np.zeros(len(g1)+1)
    for i in range(1, len(g0)):
        g0[i] = g1[i-1]/i
    g0 = g0 / np.sum(g0)
    g0[0] = p0
    g0[1:] = g0[1:]*(1-p0)
    return g0



# COMPUTATION STARTS HERE
def Psi(degree_distrb=None, initProb=1, num_gens=400, max_s=400, max_m=400, initial_T=0.8,
        intervention_gen=-1, intervention_T=0.5,
        prop_vacc=0.5, rollout_dict=None, intervention_type="none",
        custom_g0=None, custom_g1=None):
    # 3-d matrix with one matrix per generation of Psi_g
    allPsi = np.zeros(((num_gens, max_s, max_m)))
    allPsi[0][1][1] = initProb

    # Assign initial degree distribution here
    if custom_g0 is None and custom_g1 is None:
        original_degree_distrb = degree_distrb
        g1, g0 = gen_functions_with_transmissibility(original_degree_distrb,
                                                     initial_T)  # this g0 and g1 is for the G(1-(xy+1)T) in terms of the l's
    elif custom_g0 is None or custom_g1 is None:
        print('PLEASE PROVIDE BOTH CUSTOM G1 AND CUSTOM G0')
    else:
        g0 = custom_g0
        g1 = custom_g1

    M_0, M_1 = constructMatrixM(g0, g1)
    for s_g1 in range(max_s):
        for m_g1 in range(s_g1):
            allPsi[1][s_g1][m_g1] = computeLittlePsi(s_g1, m_g1, allPsi[0], M_0) #s are the rows, m are the columns

    if intervention_type=="none":
        allPsi = baseline(num_gens, max_s, max_m, allPsi, M_1)

    if intervention_type=="universal_intervention":
        allPsi = universal_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, intervention_T, allPsi, g0, M_1)

    elif intervention_type=="random_rollout":
        allPsi = random_rollout_intervention(num_gens, max_s, max_m, original_degree_distrb, initial_T, allPsi, g0, M_1, rollout_dict)

    elif intervention_type=="random":
        allPsi = random_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, prop_vacc, initial_T, allPsi, g0, M_1)

    return allPsi

def baseline(num_gens, max_s, max_m, allPsi, M_1):
    for g in range(2, num_gens):
        print('working on gen ' + str(g))
        for s in range(max_s):
            for m in range(0, s):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M_1)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g
    return allPsi

def universal_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, intervention_T, allPsi, g0, M_1):
    for g in range(2, num_gens):
        print('working on gen ' + str(g))
        if g == intervention_gen:
            new_T = intervention_T
            new_g1, new_g0 = gen_functions_with_transmissibility(original_degree_distrb, new_T)
            new_M = constructMatrixM(g0, new_g1)
            M_1 = new_M[1]
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M_1)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g
    return allPsi

def random_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, prop_vacc, initial_T, allPsi, g0, M_1):
    for g in range(2, num_gens):
        print('working on gen ' + str(g))
        if g == intervention_gen:
            new_g1 = random_vacc_distribution(original_degree_distrb, initial_T, prop_vacc)
            new_M = constructMatrixM(g0, new_g1)
            M_1 = new_M[1]
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M_1)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g

    return allPsi

def random_rollout_intervention(num_gens, max_s, max_m, original_degree_distrb, initial_T, allPsi, g0, M_1, rollout_dict):
    intervention_gen_keys = list(rollout_dict.keys())
    current_gen_idx = 0
    next_up_intervention_gen = intervention_gen_keys[current_gen_idx]
    for g in range(2, num_gens):
        print('working on gen ' + str(g))
        if g == next_up_intervention_gen:
            new_g1 = random_vacc_distribution(original_degree_distrb, initial_T, rollout_dict[next_up_intervention_gen])
            new_M = constructMatrixM(g0, new_g1)
            M_1 = new_M[1]
            if current_gen_idx < len(intervention_gen_keys) - 1:
                current_gen_idx += 1
                next_up_intervention_gen = intervention_gen_keys[current_gen_idx]
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M_1)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g

    return allPsi

def targeted_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, prop_vacc, initial_T, allPsi, g0, M_1):
    print('Working on targeted intervention')
    #TODO need to do the target vaccination distribution

    for g in range(2, num_gens):
        print('working on gen ' + str(g))
        if g == intervention_gen:
            new_g1 = targeted_vacc_distribution(original_degree_distrb, initial_T, prop_vacc)
            new_M = constructMatrixM(g0, new_g1)
            M_1 = new_M[1]
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M_1)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g

    return allPsi

def targeted_rollout_intervention(num_gens, max_s, max_m, original_degree_distrb, initial_T, allPsi, g0, M_1, rollout_dict):
    print('Working on targeted rollout intervention')
    #TODO need to do the target vaccination distribution with rollout



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


def convolve_dists(X, Y, static_length=None):
    ### X and Y should be vectors of probabilities (For our problem we need to incorporate the G0 and G1 probabilities in here.)
    ### each giving a distribution on a finite subset of the naturals
    if static_length is None:
        len_X = len(X)
        len_Y = len(Y)
    elif static_length is not None:
        len_X = static_length
        len_Y = static_length
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

def convolve_with_ifft(firstdist, seconddist):
    result_length = len(firstdist)

    # Copy each array into a 2d array of the appropriate shape.
    rows = np.zeros((2, result_length))
    for i, array in enumerate([firstdist, seconddist]):
        rows[i, :len(array)] = array

    # Transform, take the product, and do the inverse transform
    # to get the convolution.
    fft_of_rows = np.fft.fft(rows)
    fft_of_convolution = fft_of_rows.prod(axis=0)
    convolution = np.fft.ifft(fft_of_convolution)

    # Assuming real inputs, the imaginary part of the output can
    # be ignored.
    return convolution.real


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
