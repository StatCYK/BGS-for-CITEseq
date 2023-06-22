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
import sys
import os
from utility import *
import multiprocessing
import pandas as pd
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
import pickle
from itertools import product

from multiprocessing import get_context

os.chdir("./CBMC analysis/code")

backend = Aer.get_backend('aer_simulator')
backend.configuration()


nrep = 100
Input_list = []

save_path = "../results/"

for num in range(1,14):
    print(num)   
    res_test_accu = np.array(pd.read_csv(save_path + "num%d_test_MSE.csv"%num))
    res_BIC = pd.read_csv(save_path + "num%d_BIC.csv"%num)
    if num == 1:
        BIC_mat = np.array(res_BIC)
        test_accu_mat = res_test_accu
        subset_names = list(res_BIC.columns)
    else:
        BIC_mat = np.concatenate((BIC_mat,np.array(res_BIC)),axis=1)
        test_accu_mat = np.concatenate((test_accu_mat,res_test_accu),axis=1)
        subset_names.extend(list(res_BIC.columns))

m = 8
np.random.seed(1234)

for irep in range(nrep):
    Input_list = []
    for irep in range(nrep):
        BIC_irep = np.concatenate((BIC_mat[irep,:],10**40*np.ones(2**p-len(subset_names))))
        initial_idx = np.random.choice(len(subset_names),m)
        Input = dict({"initial_index":initial_idx,"criterion": BIC_irep})
        Input_list.append(Input)

pool = multiprocessing.Pool(processes=38)
experiment_res = pool.map(BGS_BIC,Input_list)


selected_test_accu =  [test_accu_mat[i, experiment_res[i][0] ] for i in range(len(experiment_res))]
selected_subset = [subset_names[ experiment_res[i][0]]   for i in range(len(experiment_res))]

np.save(save_path+"BGS_select_accu.npy",selected_test_accu)
with open(save_path+"selected_subset_BGS", "wb") as fp1:
   pickle.dump(selected_subset, fp1)



experiment_res_BSS  = pool.map(BSS_BIC,Input_list)
print(selected_subset_BSS)

selected_test_accu_BSS =  [test_accu_mat[i, experiment_res_BSS[i][0] ] for i in range(len(experiment_res_BSS))]
selected_subset_BSS = [subset_names[ experiment_res_BSS[i][0]]   for i in range(len(experiment_res_BSS))]
np.save(save_path+"BSS_select_accu.npy",selected_test_accu_BSS)
with open(save_path+"selected_subset_BSS", "wb") as fp2:
    pickle.dump(selected_subset_BSS, fp2)

