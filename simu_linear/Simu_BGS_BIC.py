# Importing standard Qiskit libraries and configuring account
from qiskit import QuantumCircuit, execute, Aer, IBMQ
from qiskit.compiler import transpile, assemble
from qiskit.visualization import *
from qiskit import IBMQ, QuantumRegister, ClassicalRegister
import numpy as np
from qiskit import BasicAer
from qiskit.tools.visualization import plot_histogram
import time
from qiskit.tools.monitor import job_monitor
from datetime import datetime
from utility import *
import multiprocessing
import pandas as pd
from sklearn import linear_model
from sklearn.model_selection import train_val_split
from sklearn.preprocessing import scale
import pickle
from itertools import product

#### backend setting for quantum system ####

IBMQ.load_account()
provider = IBMQ.get_provider(hub = "ibm-q", group='open', project='main')
backend = Aer.get_backend('aer_simulator')#('ibmq_manila')
backend.configuration()


#### data generating ####
N=1000
nrep = 100
snr=1
rho=0.7

for p in range(6,16):
    s=round(p/2)
    B = np.zeros((p,))
    B[0:s]=1
    Sigma = np.zeros((p,p))
    for i in range(p):
        for j in range(p):
            Sigma[i,j] = rho**abs(i-j)
    mu = np.zeros((p,))
    Input_list = []
    test_list = []
    for rep in range(nrep):
        np.random.seed(101+rep)
        X_train = np.random.multivariate_normal(mu, Sigma, N)+1
        epsion = np.random.normal(0,1,N)
        Y_train  = X_train.dot(B)+np.sqrt(B.dot(Sigma).dot(B)/snr)*epsion
        initial_candidate = [] # candidates of initial subset
        m = 5
        for i in range(m):
            initial_candidate.append([i for i in np.random.binomial(1, 0.5, p)])
        Input = dict({"Y_train":Y_train,"X_train":X_train,\
                    "initial_index":initial_candidate,"criterion": "Rsquared"})
        Input_list.append(Input)

    with open("./data/simu_Input_p%d"%p, "wb") as fp1:
        pickle.dump(Input_list, fp1)

#### experiment running for BGS ####
for p in range(6,16):
    print(p)
    Input_list = []
    with open("./data/simu_Input_p%d"%p, 'rb') as fr3:
        Input_list.extend(pickle.load(fr3))
    pool = multiprocessing.Pool(processes=36)
    experiment_res = pool.map(BGS_linearReg_BIC,Input_list)
    with open("./output/res_simu_BGS_p%d"%p, "wb") as fp4:
        pickle.dump(experiment_res,fp4)

