import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib as m
import matplotlib.colors as colors


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

def binomial_degree_distb(N):
    degree_dist = np.zeros(40)
    p = 6/N
    for k in range(0, len(degree_dist)):
        p_k = (p**k)*((1-p)**(N-k))*math.comb(N, k)
        degree_dist[k] = p_k
    return degree_dist





def gen_functions_with_transmissibility(degree_distrb, T):
    # Given a degree distribution for G0 (or the degree distribution of the entire network).
    # Given transmissibility T
    maxk = len(degree_distrb)
    p_k = degree_distrb
    p_LK = np.zeros((400, 400))

    # Generates pgf in variable l as probabilities of infection of l neighbors given the original degree distribution and transmission prob T
    for k in range(0, maxk):
        for l in range(0, k + 1):
            try:
                p_LgivenK = p_k[k] * (
                    math.gamma(k + 1) / (math.gamma(l + 1) * math.gamma(k - l + 1)) * T ** (l) * (1 - T) ** (k - l))
                p_LK[k][l] = p_LgivenK
            except OverflowError:
                p_LK[k][l] = 0
    p_l = np.sum(p_LK, axis=0)
    p_l = p_l / (np.sum(p_l))

    G0_with_T = pdf_of(p_l)
    G1_with_T = g1_of(
        G0_with_T)
    return G1_with_T, G0_with_T


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

    M_1[0][0]=1
    return (M_0, M_1)

def computeLittlePsi(s, m, prevGenPsi, M):
    s_prime = s - m
    newPsi = prevGenPsi[s_prime].dot(M[:, m])
    return newPsi


def Psi(degree_distrb, initProb, num_gens, max_s, max_m, initial_T, intervention_gen=-1, intervention_T=0.5):
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
        print('working on gen '+str(g))
        # If g is intervention, re-call g0_l's, g1_l's, M0, M1 etc
        if g == intervention_gen:
            new_T = intervention_T
            new_g1, new_g0 = gen_functions_with_transmissibility(original_degree_distrb, new_T)
            new_M = constructMatrixM(new_g0, new_g1)
            M_1 = new_M[1]
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M_1)
        psi_g = allPsi[g]
        psi_g = psi_g/np.sum(psi_g)
        allPsi[g] = psi_g
    return allPsi

def newFigure():
    # TODO move into "do" in this file, "visualize" in next file so it saves results locally and then can be easily played with
    s = 400
    g = 100
    T = 0.8
    heatmap_m = np.zeros((g, s))
    heatmap_s = np.zeros((g, s))
    # x-axis gens g
    # y axis: s, number cumulative infections
    # each (s,g) square corresponds to: probability that s are infected at generation g, as a number
    # instead of a curve you get a heat map column
    # TODO
    # Plot the whole psi matrix for each generation
    # compute the s-marginal (using the code for the slice distribution)
    # using each distribution which is 1 g, over all s, then make that one row of the matrix for row 'g'
    initProb = 1
    power_law_degree_dist = power_law_degree_distrb(s)
    # binomial_degree_dist = binomial_degree_distb(1000)
    all_psi_results = Psi(power_law_degree_dist, initProb, g, s, s, T)
    # Specify intervention parameters for gen_intervene and T_intervene after initial T:
    # all_psi_results_with_intervention = Psi(power_law_degree_dist, initProb, num_gens, num_nodes, num_nodes, T, 3, 0.4)
    color_key = {2: 'blue', 6: 'red', 11: 'orange', 18: 'black'}
    for gen in range(1, g):
        inverted_s_m = all_psi_results[gen].T
        s_marginal = np.sum(inverted_s_m, axis=0)
        m_marginal = np.sum(inverted_s_m, axis=1)
        s_marginal = s_marginal/np.sum(s_marginal) #normalize
        m_marginal = m_marginal/np.sum(m_marginal) #normalize
        heatmap_m[gen] = m_marginal
        heatmap_s[gen] = s_marginal
        # label='$g='+str(gen)+'$'
        # color = color_key[gen]
        # plt.plot(ps_g_analytical[2:], label=label, color=color, linestyle='-')
        # plt.title('No intervention')
    # plt.show()



    cmap = plt.cm.hot(np.linspace(1, 0, 100000))
    cmap = m.colors.ListedColormap(cmap[:, :-1])

    fig, ax = plt.subplots()
    ax.axis('equal')
    # ax.imshow(psi_g[:60][:, :100], cmap=cmap, norm=colors.PowerNorm(gamma=0.05, vmin=0, vmax=max(psi_g[0])), label='$gen='+str(gen)+'$')  # gen 5
    ax.imshow(heatmap_m.T, cmap=cmap, norm = plt.cm.colors.SymLogNorm(linthresh=0.00005, vmax=0.4, vmin=0.000), aspect='auto')
    # ax.imshow(psi_g[:80][:, :150], cmap=cm)  # gen 5
    ax.invert_yaxis()
    # plt.title('Phase Space at Generation '+str(gen)+' of '+str(title_label))
    plt.ylabel('$m$', fontsize=16)
    plt.xlabel('$gen$', fontsize=16)
    plt.xticks(np.arange(0, g, 10), np.arange(0, g, 10))
    # plt.legend(loc='upper right')
    # plt.savefig(str(title_label)+str(gen)+'_phase_space.png')
    plt.show()

    fig, ax = plt.subplots()
    ax.axis('equal')
    # ax.imshow(psi_g[:60][:, :100], cmap=cmap, norm=colors.PowerNorm(gamma=0.05, vmin=0, vmax=max(psi_g[0])), label='$gen='+str(gen)+'$')  # gen 5
    ax.imshow(heatmap_s.T, cmap=cmap, norm = plt.cm.colors.SymLogNorm(linthresh=0.00005, vmax=0.4, vmin=0.000), aspect='auto')
    # ax.imshow(psi_g[:80][:, :150], cmap=cm)  # gen 5
    ax.invert_yaxis()
    # plt.title('Phase Space at Generation '+str(gen)+' of '+str(title_label))
    plt.ylabel('$s$', fontsize=16)
    plt.xlabel('$gen$', fontsize=16)
    plt.xticks(np.arange(0, g, 10), np.arange(0, g, 10))
    # plt.legend(loc='upper right')
    # plt.savefig(str(title_label)+str(gen)+'_phase_space.png')
    plt.show()

def phaseSpace(num_gens, num_nodes):
    # TODO toggle saving results
    # the generating function for psi gen g (prob of having s infected by the end of gen g of which m became infected during gen g
    initProb = 1
    binom = binomial_degree_distb(1000)
    all_psi_results = Psi(binom, initProb, num_gens, num_nodes, num_nodes, 0.2)
    all_psi_results_with_intervention = Psi(binom, initProb, num_gens, num_nodes, num_nodes, 0.2, 3, 0.1)

    # np.savetxt('../pgf-nets-data/binom_allPsiT8_2_int.txt', all_psi_results_with_intervention[2], delimiter=',')
    # np.savetxt('../pgf-nets-data/binom_allPsiT8_6_int.txt', all_psi_results_with_intervention[6], delimiter=',')
    # np.savetxt('../pgf-nets-data/binom_allPsiT8_11_int.txt', all_psi_results_with_intervention[11], delimiter=',')
    # np.savetxt('../pgf-nets-data/binom_allPsiT8_18_int.txt', all_psi_results_with_intervention[18], delimiter=',')
    # np.savetxt('../pgf-nets-data/binom_allPsiT8_2_reg.txt', all_psi_results[2], delimiter=',')
    # np.savetxt('../pgf-nets-data/binom_allPsiT8_6_reg.txt', all_psi_results[6], delimiter=',')
    # np.savetxt('../pgf-nets-data/binom_allPsiT8_11_reg.txt', all_psi_results[11], delimiter=',')
    # np.savetxt('../pgf-nets-data/binom_allPsiT8_18_reg.txt', all_psi_results[18], delimiter=',')


    initProb = 1
    power_law = power_law_degree_distrb(400)
    all_psi_results = Psi(power_law, initProb, num_gens, num_nodes, num_nodes, 0.8)
    all_psi_results_with_intervention = Psi(power_law, initProb, num_gens, num_nodes, num_nodes, 0.8, 3, 0.4)
    # TODO save all the results and return

    return all_psi_results


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
