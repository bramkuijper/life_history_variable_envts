#!/usr/bin/env python

import numpy as np

mu_a_tau = [ 0, 0.02 ]
mu_a_eps = [ 0, 0.02 ]

sdmu = 0.01

m_r = [ 0.05, 0.1 ]
m_n = [ 0.05 ]

b = [ 2,3,4 ]

clutch_multiplier = [ 1.0, 2.0 ]

min_clutch = [ 0, 1.0 ]

ampl = [ 0, 0.1, 0.5, 1.0 ]

stoch = [ 0, 0.1, 0.5, 1.0 ]

exe = "./xseasonality_thr"

sd_envt = 1.0

nrep = 5

ctr = 0

for mu_a_tau_i in mu_a_tau:
    for mu_a_eps_i in mu_a_eps:
        for m_r_i in m_r:
            for m_n_i in m_n:
                for b_i in b:
                    for clutch_multiplier_i in clutch_multiplier:
                        for min_clutch_i in min_clutch:
                            for ampl_i in ampl:
                                for stoch_i in stoch:
                                    for nrep_i in range(0,nrep):

                                        ctr += 1
                                        print("echo " + str(ctr))
                                        print(exe + " 0 " + str(mu_a_tau_i) + " " + str(mu_a_eps_i) + " " 
                                                + str(sdmu) + " " +  str(m_r_i) + " " + str(m_n_i) + " " 
                                                + str(b_i) + " " + str(clutch_multiplier_i) + " " 
                                                + str(min_clutch_i) + " " + str(ampl_i) + " " + str(stoch_i) + " "
                                                + str(sd_envt))
