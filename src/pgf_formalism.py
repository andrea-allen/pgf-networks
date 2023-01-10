import math
import numpy as np
from src import gen_extinct_prob, manatee

"""
For use, call compute_phase_space from outside file with associated parameters to obtain and save
matrix files of the s,m phase space for each generation specified.
"""

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
    d_len = len(deg_dist)
    psi = Psi(deg_dist, initProb=1, num_gens=50, max_s=d_len, max_m=d_len, initial_T=T)
    expected_cum_array = np.zeros(50)
    for gen in range(50):
        psi[gen][:,0] = np.zeros(d_len)
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
                        do_non_interv=True, do_interv=True, intervention_type="none", pre_vax_correction=False):
    """
    Computes and saves a matrix for each generation of the s,m phase space
    :param num_gens: Number of generations to compute til
    :param num_nodes: Number of nodes in the network. Should match the length of the passed degree_distribution
    :param degree_distribution: Vector or list of length num_nodes.
    :param transmissibility: Original transmission paramter T
    :param save_results: boolean, will save files to directory and format in file_fmt_to_save. Must specify dyamic arg param for generation, {0}.
    :param gens_to_save: Which generations to save a file for
    :param file_fmt_to_save: Directory and file name pattern, with arg for gen number, such as my_files/generation_{0}
    :param intervention_gen: FOR SINGLE, UNIVERSAL INTERVENTION ONLY: Generation of intervention
    :param intervention_trans: FOR SINGLE, UNIVERSAL INTERVENTION ONLY: Transmissiblity change at intervention
    :param vacc_pop: FOR SINGLE, UNIVERSAL INTERVENTION ONLY: population vaccinated at intervention
    :param rollout_dict: Non-cumulative dict of generations and their vaccination percentages, e.g. {3: .05, 4: .01, 5:.02}
    :param do_non_interv: boolean. Will run non-intervention case on same degree distribution and transmissibility.
    :param do_interv: boolean. Will run intervention model.
    :param intervention_type: string. Specify universal, random_rollout, or targeted_rollout
    :param pre_vax_correction: IGNORE, DO NOT USE
    :return:
    """
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
                                                transmissibility, intervention_gen, intervention_trans, vacc_pop,
                                                rollout_dict, intervention_type, pre_vax_correction=pre_vax_correction)
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
    if np.sum(degree_list) == 0:
        return np.array(degree_list)
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
    if np.sum(g_0) == 0:
        return g_1
    return g_1 / (z1_of(g_0))

def g_g_of(g_gminus1, deltas, g):
    g_g = np.zeros(len(g_gminus1))
    denom = np.zeros(len(g_gminus1))

    for k in range(len(g_gminus1) - 1):
        g_g[k] = k * (k+1) * (1 - deltas[g][k+1]) * g_gminus1[k + 1]
        denom[k] = (k + 1) * (1 - deltas[g][k+1]) * g_gminus1[k + 1]
    if np.sum(denom) == 0
        return g_g
    return g_g / np.sum(denom)



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
    if np.sum(p_l) != 0:
        p_l = p_l / (np.sum(p_l))

    G0_with_T = pdf_of(p_l)
    G1_with_T = g1_of(G0_with_T)
    return G1_with_T, G0_with_T


def critical_degree_calc(prop, degree_d, delta_k, g): # This function is the change
    k_crit = len(degree_d) #default for the critical k value
    temp_prop = prop
    while temp_prop > 0:
        temp_prop = temp_prop - degree_d[k_crit-1]
        delta_k[g][k_crit-1] = 1
        k_crit -= 1

    delta_k[g][k_crit] = (1 + temp_prop) / degree_d[k_crit]

    return k_crit, -temp_prop, delta_k

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
        custom_g0=None, custom_g1=None, pre_vax_correction=False):
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
        g0 = custom_g0 # Custom g0 and g1 set if using a distribution for secondary cases (such as a negative binomial) and not the original network degree distribution
        g1 = custom_g1

    M_0, M_1 = constructMatrixM(g0, g1)

    if intervention_type=="none":
        allPsi = baseline(num_gens, max_s, max_m, allPsi, M_1, M_0)

    if intervention_type=="universal":
        allPsi = universal_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, intervention_T, initial_T, allPsi, g0, M_1, M_0, pre_vax_correction)

    elif intervention_type=="random_rollout":
        allPsi = random_rollout_intervention(num_gens, max_s, max_m, original_degree_distrb, initial_T, allPsi, g0, M_1, M_0, rollout_dict, pre_vax_correction)

    elif intervention_type=="random":
        print('PLEASE USE RANDOM ROLLOUT')
        # allPsi = random_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, prop_vacc, initial_T, allPsi, g0, M_1, M_0, pre_vax_correction)

    elif intervention_type=="targeted":
        print('PLEASE USE TARGETED ROLLOUT')
        # allPsi = targeted_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, prop_vacc, initial_T, allPsi, g0, M_1, M_0, pre_vax_correction)

    elif intervention_type=="targeted_rollout":
        allPsi = targeted_rollout_intervention(num_gens, max_s, max_m, original_degree_distrb, initial_T, allPsi, g0, M_1, M_0, rollout_dict, pre_vax_correction)

    else:
        print(f"Specified an unrecognized intervention type of {intervention_type}, please specify one of universal, random, targeted, random_rollout or targeted_rollout")

    return allPsi

def baseline(num_gens, max_s, max_m, allPsi, M_1, M_0):
    for g in range(1, num_gens):
        M = M_1
        if g == 1:
            M = M_0
        print('working on gen ' + str(g))
        for s in range(max_s):
            for m in range(0, s):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g
    return allPsi

def universal_intervention(num_gens, max_s, max_m, original_degree_distrb, intervention_gen, intervention_T, initial_T, allPsi, g0, M_1, M_0, pre_vax_correction):
    for g in range(1, num_gens):
        if pre_vax_correction:
            if g < intervention_gen:
                # Re-compute T_g based on beta, gamma, g, v, and g_int
                beta_1 = .5  # TODO is it ok that this is arbitrary? Ans: Yes.
                gamma_1 = (beta_1 - beta_1 * initial_T) / initial_T
                q_1 = z1_of(g1_of(original_degree_distrb))
                beta_2 = .5
                gamma_2 = (beta_2 - beta_2 * intervention_T) / intervention_T

                T_g_i = T_pre_vax_fancy_2(beta_1=beta_1, beta_2=beta_2, gamma_1=gamma_1, gamma_2=gamma_2, q_1=q_1, gen_i=intervention_gen-g,
                                  v=1) #should v be 1 or 0?
                # Use new T_g to compute: new secondary degree distributions with T G(x;T):
                g1_T, g0_T = gen_functions_with_transmissibility(original_degree_distrb, T_g_i)
                #                                                      new T_g)
                # Re-compute M0 and M1 from new g1, g0
                M_0_g, M_1_g = constructMatrixM(g0_T, g1_T)
                M = M_1_g
                if g == 1:
                    M = M_0_g
        else:
            if g < intervention_gen:
                M = M_1
                if g == 1:
                    M = M_0
        print('working on gen ' + str(g))
        if g == intervention_gen:
            new_T = intervention_T
            new_g1, new_g0 = gen_functions_with_transmissibility(original_degree_distrb, new_T)
            new_M = constructMatrixM(g0, new_g1)
            M_1 = new_M[1]
            M = M_1
            if g == 1:
                M = M_0
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g
    return allPsi


def random_rollout_intervention(num_gens, max_s, max_m, original_degree_distrb, initial_T, allPsi, g0, M_1, M_0,
                                rollout_dict, pre_vax_correction):
    intervention_gen_keys = list(rollout_dict.keys())
    inter_list = [0]
    inter_list.extend(intervention_gen_keys)
    beta = 0.8
    gamma = (beta - beta*initial_T) / initial_T
    q_1 = z1_of(g1_of(original_degree_distrb))
    v_rollout_cumu = [0]
    v_rollout_cumu.extend(np.cumsum(list(rollout_dict.values())))
    vacc_vector = np.zeros(num_gens)
    for i in range(num_gens):
        if i in inter_list and i!=0:
            vacc_vector[i] = rollout_dict[i]+vacc_vector[i-1]
        else:
            vacc_vector[i] = vacc_vector[i-1]

    for g in range(1, num_gens):
        # change num_gens to a big big number in computation and then see what the plot looks like
        T_g = manatee.t_of_g(betas=np.full(num_gens, beta), gammas=np.full(num_gens, gamma), qs=np.full(num_gens, q_1),
                             vaccs=vacc_vector, g=g)
        print(g, T_g) # this makes sense
        g1_T, g0_T = gen_functions_with_transmissibility(original_degree_distrb, T_g) # changing it to modify the new G1 with the new Tg?
        M_0_g, M_1_g = constructMatrixM(g0_T, g1_T) # this now might be wiping away the changes made with the random vax dist
        M = M_1_g

        if g == 1:
            M = M_0
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)  # normalize
        allPsi[g] = psi_g

    return allPsi

def modify_g0(G0, k_crit, replacement):
    G0_x = np.zeros((len(G0)))
    for k in range(0, k_crit):
        G0_x[k] = G0[k]
    G0_x[k_crit] = replacement
    if np.sum(G0_x) > 0:
        G0_x = G0_x/np.sum(G0_x)
    return G0_x



def targeted_rollout_intervention(num_gens, max_s, max_m, original_degree_distrb,
                                  initial_T, allPsi, g0, M_1, M_0, rollout_dict,
                                  pre_vax_correction):
    # First, need to create and find H, q and T for each generation g in the rollout:
    intervention_gen_keys = list(rollout_dict.keys())
    inter_list = [0]
    inter_list.extend(intervention_gen_keys)
    v_rollout_cumu = [0]
    v_rollout_cumu.extend(np.cumsum(list(rollout_dict.values())))
    vacc_vector = np.zeros(num_gens)
    for i in range(num_gens):
        if i in inter_list and i!=0:
            vacc_vector[i] = rollout_dict[i]+vacc_vector[i-1]
        else:
            vacc_vector[i] = vacc_vector[i-1]
    beta = 0.8 # arbitrary, via choice of initial_T
    gamma = (beta - beta*initial_T) / initial_T
    maxk = len(original_degree_distrb)
    g1_orig = g1_of(original_degree_distrb)
    origin_q = np.sum([k*g1_orig[k] for k in range(len(g1_orig))])
    dynamic_q = np.zeros((num_gens))
    first_gen_inter = intervention_gen_keys[0]
    dynamic_q[:first_gen_inter] = origin_q
    dynamic_H = np.zeros((num_gens))
    store_Ggs = np.zeros((maxk, maxk))
    store_G0s = np.zeros((maxk, maxk))
    store_Ggs[:first_gen_inter] = g1_orig
    store_G0s[:first_gen_inter] = original_degree_distrb
    delta_g_k = np.ones((num_gens,maxk))

    ##### finding the q's and H's for each intervention:
    for gen, V in rollout_dict.items():
        print(vacc_vector[gen])
        crit_value, replace_prob, delta_g_k = critical_degree_calc(vacc_vector[gen], original_degree_distrb, delta_g_k, gen) # check if original makes sense
        print(crit_value)
        #H = np.sum([(k) * original_degree_distrb[k] for k in range(crit_value, len(original_degree_distrb))])\
        #    /np.sum([(k) * original_degree_distrb[k] for k in range(0, len(original_degree_distrb))])
        H_g = np.sum([(k+1) * delta_g_k[gen][k] * original_degree_distrb[k+1] for k in range(0, len(original_degree_distrb))])\
                \np.sum([(k+1) * original_degree_distrb[k+1] for k in range(0, len(original_degree_distrb))])
        print(f'H:{H_g}')
        mod_G0_H = modify_g0(original_degree_distrb, crit_value, replace_prob)
        print(mod_G0_H)
        mod_Gg_H = g1_of(mod_G0_H)
        store_Ggs[gen:] = mod_Gg_H
        store_G0s[gen:] = mod_G0_H
        q_g = (1-H_g)* np.sum(mod_Gg_H) #change to the sum since the derivations are different for normalizing
        if np.isnan(q_g):
            q_g = 0
        print(q_g)
        dynamic_q[gen:] = q_g
        dynamic_H[gen:] = H_g

    ######

    for g in range(1, num_gens):
        T_g = manatee.t_of_g(betas=np.full(num_gens, beta), gammas=np.full(num_gens, gamma), qs=dynamic_q,
                             vaccs=dynamic_H, g=g)
        print(g, T_g)
        g0_x = store_G0s[g]
        g1_T, g0_T = gen_functions_with_transmissibility(g0_x, T_g)
        M_0_g, M_1_g = constructMatrixM(g0_T, g1_T)
        M = M_1_g

        if g == 1:
            M = M_0
        for s in range(max_s):
            for m in range(max_m):
                allPsi[g][s][m] = computeLittlePsi(s, m, allPsi[g - 1], M)
        psi_g = allPsi[g]
        psi_g = psi_g / np.sum(psi_g)
        allPsi[g] = psi_g

    return allPsi

def T_pre_vax(beta_1, beta_2, gamma_1, gamma_2, q_1, gen_i, v):
    T_gen_i = beta_1 / (beta_1 + gamma_1 + ((q_1*beta_1) / gen_i)) \
              + (((q_1*beta_1) / gen_i) / (beta_1 + gamma_1 + ((q_1*beta_1) / gen_i))
                 * (beta_2/ (beta_2+gamma_2)) * (1-v))
    return T_gen_i

def T_pre_vax_fancy_2(beta_1, beta_2, gamma_1, gamma_2, q_1, gen_i, v):
    T_gen_i = beta_1 / (beta_1 + gamma_1 + ((q_1*beta_1) / (gen_i - .5))) \
              + (((q_1*beta_1) / (gen_i - .5)) / (beta_1 + gamma_1 + ((q_1*beta_1) / (gen_i-0.5)))
                 * (beta_2/ (beta_2+gamma_2)) * (1-v))
    return T_gen_i

def T_pre_vax_rollout(beta_1, beta_2, gamma_1, gamma_2, q_1, current_g, rollout_dict):
    T_modified = 0
    vax_prop_cum = 0
    for g in rollout_dict.keys():
        vax_prop_cum += rollout_dict[g]
        i_w = g - current_g
        if (T_modified == 0) and (g > current_g): #i.e. first term has not been added yet, and found next soonest intervention, then add 1st term:
            T_modified += beta_1 / (beta_1+gamma_1 + (q_1*beta_1/i_w))
        if g > current_g:
            T_modified += ((beta_1 * (1-vax_prop_cum)) / (beta_1 + gamma_1 + (q_1*beta_1/(i_w)))) * ((q_1*beta_1/i_w)/(beta_1+gamma_1+(q_1*beta_1/i_w)))
    return T_modified


def T_pre_vax_fancy(beta_1, beta_2, gamma_1, gamma_2, q_1, gen_i, v):
    term_1 = -math.e ** (-1 * gamma_1 * gen_i /(q_1 * beta_1)) + ((gamma_1/(beta_1+gamma_1))*math.e**(-(beta_1+gamma_1)*gen_i/(q_1*beta_1))) + (beta_1/(beta_1+gamma_1))
    term_2 = math.e ** (-gamma_2*gen_i/(q_1*beta_1)) - ((gamma_2/(gamma_2+beta_2)) * math.e**(-(beta_2+gamma_2)*gen_i/(q_1*beta_1)))
    T_gen_i = term_1 + (1-v) * term_2
    return T_gen_i



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
