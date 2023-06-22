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
sys.path.insert(0, '~/Projects/Quantum/pbmc_HIV')
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


CellType_dirs = os.listdir("./BIC_result")

for Dir in CellType_dirs:
    save_path = "./select_result/"+Dir
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    for num in range(1,10):
        print(num)
        res_BIC = pd.read_csv("./BIC_result/"+Dir+ "/num%d_BIC.csv"%num)
        res_values = np.array(res_BIC)
        nrep = res_values.shape[0]
        subset_names = res_BIC.columns
        BIC_mat = np.array(res_BIC)
        p = int(np.ceil(np.log2(len(subset_names))))
        m = 8
        np.random.seed(1234)
        Input_list = []
        for irep in range(nrep):
            BIC_irep = np.concatenate((BIC_mat[irep,:],10**40*np.ones(2**p-len(subset_names))))
            initial_idx = np.random.choice(len(subset_names),m)
            Input = dict({"initial_index":initial_idx,"criterion": BIC_irep})
            Input_list.append(Input)
        
        with open(save_path+"/Input", "wb") as fp1:
            pickle.dump(Input_list, fp1)
        
        pool = multiprocessing.Pool(processes=38)
        experiment_res = pool.map(BGS_BIC,Input_list)

        with open(save_path+ "/res_CiteSeq_BGS_p_%d"%num, "wb") as fp3:
            pickle.dump(experiment_res,fp3)
        print(sum([res[3] for res in experiment_res]))

    res = pd.DataFrame({})

    for num in range(1,10):
        experiment_res  = []
        with open(save_path+"/res_CiteSeq_BGS_p_%d"%num, 'rb') as fr3:
            experiment_res.extend(pickle.load(fr3))
        res_test_accu = np.array(pd.read_csv("./BIC_result/"+Dir+ "/num%d_test_auc.csv"%num))
        selected_test_accu =  [res_test_accu[i, experiment_res[i][0]] for i in range(len(experiment_res))]
        np.save(save_path+"/BGS_select_num%d_auc.npy"%num,selected_test_accu)
        res["num%d"%num] = selected_test_accu

    res.to_csv(save_path+"/BGS_select_auc.csv")