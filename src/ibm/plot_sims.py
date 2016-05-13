#!/usr/bin/env python3

import sys, re, os.path
import pandas as pd
import matplotlib

import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# should we plot each nth generation? 
skip = 1

# get the filename from the command line
filename = sys.argv[1]

# read the file and look for parline
f = open(filename);
fl = f.readlines();
f.close()


parline = -1

for idx, line in enumerate(fl):
    if re.match("^mu_a_eps.*",line) != None:
        parline = idx - 1;
        break;

parameters = None

# read in the csv file
if parline > 0:
    histdat = pd.read_csv(filename, nrows=parline-2, sep=";")
    parameters = pd.read_csv(filename, skiprows=parline, sep=";", header=None)
    parameters.columns = ['name','value','x']
else:
    histdat = pd.read_csv(filename, sep=";")


# only take every tenth generation, otherwise too much data....
histdat = histdat[histdat["time"] % skip == 0]


# generate the figure

# initialize and specify size 
fig = plt.figure(figsize=(10,20))

num_rows = 8

# prob of becoming a breeder
plt.subplot(num_rows,1,1)
plt.plot(histdat["time"],histdat["environment"],'b',
        histdat["time"],histdat["mean_p"],'r',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Breeding prob, $\bar{p}_{t}$')
plt.legend((r'$\varepsilon_{t}$',r'$\bar{p}_{t}$'))

# mean evolving values
plt.subplot(num_rows,1,2)
plt.plot(
        histdat["time"],histdat["mean_a_tau"],'c',
        histdat["time"],histdat["mean_a_eps"],'m',
        histdat["time"],histdat["mean_b_tau"],'y',
        histdat["time"],histdat["mean_b_eps"],'k',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'evolving loci')
plt.legend((r'$\bar{a}_{\tau,t}$',r'$\bar{a}_{\varepsilon,t}$',r'$\bar{b}_{\tau,t}$',r'$\bar{b}_{\varepsilon,t}$'))

# variances
plt.subplot(num_rows,1,3)
plt.plot(
        histdat["time"],histdat["var_a_tau"],'c',
        histdat["time"],histdat["var_a_eps"],'m',
        histdat["time"],histdat["var_b_tau"],'y',
        histdat["time"],histdat["var_b_eps"],'k',
        histdat["time"],histdat["var_p"],'g',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Variances')
plt.legend((r'$\sigma_{a_{\tau},t}^2$',r'$\sigma_{a_{\varepsilon},t}^2$',r'$\sigma_{b_{\tau},t}^2$',r'$\sigma_{b_{\varepsilon},t}^2$',r'$\sigma_{p,t}^2$'))

# numbers of individuals
plt.subplot(num_rows,1,4)
plt.plot(histdat["time"],histdat["Nrep"],'#129aff',
        histdat["time"],histdat["Nnonrep"],'#a60090',
        histdat["time"],histdat["Noff"],'#c74c4c',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='on')
plt.ylabel(r'$n_{i}$')
plt.legend((r'$N_{r}$',r'$N_{n}$',r'$N_{\mathrm{off}}$',r'$N_{\mathrm{recov}}$'))

nrow_histdat = histdat.shape[0]

histdat_sub = histdat.iloc[range(nrow_histdat-100,nrow_histdat,1),:]

# zoomed in bit

# prob of becoming a breeder
plt.subplot(num_rows,1,5)
plt.plot(histdat_sub["time"],histdat_sub["environment"],'b',
        histdat_sub["time"],histdat_sub["mean_ntime"],'g',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='on')
plt.ylabel(r'Breeding prob, $\bar{p}_{t}$')
plt.legend((r'$\varepsilon_{t}$',r'$\bar{\tau}_{t}$'))

# prob of becoming a breeder
plt.subplot(num_rows,1,6)
plt.plot(
        histdat_sub["time"],histdat_sub["mean_p"],'r',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='on')
plt.ylabel(r'Breeding prob, $\bar{p}_{t}$')
plt.ylim((0,1))

plt.subplot(num_rows,1,7)
plt.plot(histdat_sub["time"],histdat_sub["Nrep"],'#129aff',
        histdat_sub["time"],histdat_sub["Nnonrep"],'#a60090',
        histdat_sub["time"],histdat_sub["Noff"],'#c74c4c',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='on')
plt.ylabel(r'$n_{i}$')
plt.legend((r'$N_{r}$',r'$N_{n}$',r'$N_{\mathrm{off}}$',r'$N_{\mathrm{recov}}$'))

def thr(h, envt):
    if envt > h:
        return 1

    return -1


if parameters is not None:
    # make sample timeseries
    dat = pd.DataFrame(columns=['t','envt','p'])
    dat['t'] = pd.Series(list(np.arange(0,100,1)))

    def fg(row):
        global parameters

        ampl = double(parameters[parameters['name'] == "ampl"]["value"])
        stoch = double(parameters[parameters['name'] == "stoch"]["value"])
        sd_envt = double(parameters[parameters['name'] == "sd_envt"]["value"])
        eps = ampl * math.sin(2 * math.pi / 12 * row['t']) + stoch * np.random.normal(0, sd_envt)
        
        nrow = histdat.shape[0] - 1

        # get values for the evolving traits
        a_tau = histdat.ix[nrow,"mean_a_tau"]
        a_eps = histdat.ix[nrow,"mean_a_eps"]
        b_tau = histdat.ix[nrow,"mean_b_tau"]
        b_eps = histdat.ix[nrow,"mean_b_eps"]

        p = 1.0 / (1.0 + exp(-(a_tau + b_tau * row['t']) - (a_eps + b_eps * eps)))

        return(pd.Series(dict(envt=eps,p=p)))

    dat[['envt','p']] = dat.apply(fg,axis=1)

    # sample environmental response
    plt.subplot(num_rows,1,8)
    plt.plot(dat["t"],dat["p"],'g',
            dat["t"],dat["envt"],'b',
            linewidth=1)
    plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
    plt.ylabel(r'$p_{t}$')
    plt.legend((r'$p_{t}$',r'$\varepsilon_{t}$'))

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
