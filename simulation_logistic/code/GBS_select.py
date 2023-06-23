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
os.chdir('./simu_logistic/code')
from utility import *
import multiprocessing
import pandas as pd
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
import pickle
from itertools import product
import os

backend = Aer.get_backend('aer_simulator')
backend.configuration()


res = pd.DataFrame({})
for p in range(3,14):
    save_path = "./select_result/dim_%d"%p
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    nrep = 100
    m = 8
    np.random.seed(1234)
    Input_list = []
    for irep in range(1,101):
        res_BIC = pd.read_csv("./BIC_result/dim%d_sim%d.csv"%(p,irep),index_col=0)
        subset_names = list(res_BIC.subset)
        BIC_values = np.array(res_BIC.BIC)
        BIC_irep = np.concatenate((BIC_values,10**40*np.ones(2**p-len(subset_names))))
        initial_idx = np.random.choice(len(subset_names),m)
        Input = dict({"initial_index":initial_idx,"criterion": BIC_irep})
        Input_list.append(Input)
    
    with open(save_path+"/Input", "wb") as fp1:
        pickle.dump(Input_list, fp1)

    pool = multiprocessing.Pool(processes=38)
    experiment_res = pool.map(BGS_BIC,Input_list)

    with open(save_path+ "/res_BGS", "wb") as fp3:
        pickle.dump(experiment_res,fp3)
     
    selected_test_accu =  [experiment_res[i][3] for i in range(len(experiment_res))]
    np.save(save_path+"/BGS_select_accu.npy",selected_test_accu)
    
    res["dim%d"%p] = [experiment_res[i][0] for i in range(len(experiment_res))]

res.to_csv("./select_result_all_size/BGS_select.csv")
