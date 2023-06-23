# Importing standard Qiskit libraries and configuring account
from qiskit import QuantumCircuit, execute, Aer, IBMQ
from qiskit.compiler import transpile, assemble
#from qiskit.tools.jupyter import *
#from qiskit.visualization import *
from qiskit import IBMQ, QuantumRegister, ClassicalRegister, execute, QuantumCircuit
# Loading your IBM Q account(s)
#import pylab
import numpy as np
from qiskit import BasicAer
from qiskit.tools.visualization import plot_histogram
import time
from qiskit.tools.monitor import job_monitor
from datetime import datetime
from sklearn.linear_model import LinearRegression
from qiskit.quantum_info.operators import Operator, Pauli
from heapq import nlargest

######## backend setting ###############

#IBMQ.load_account()
#provider = IBMQ.get_provider(hub = "ibm-q-research-2", group='uni-georgia-1', project='main')
#provider = IBMQ.get_provider(hub = "ibm-q", group='open', project='main')
#backend = provider.get_backend('ibmq_qasm_simulator')#('ibmq_manila')
#backend.configuration()
backend = Aer.get_backend('aer_simulator')#provider.get_backend('ibmq_qasm_simulator')#('ibmq_manila')
backend.configuration()

def encode(state):
    state_index = int(state,2)
    state = ["0" for i in range(len(initial_oracle))]
    state[state_index]="1"
    encoded_state = ''.join(state)
    return(encoded_state)

def State2Bin(state):
    p = int(np.log2(len(state)))
    state_index = state.index('1')
    Bin = bin(state_index).replace('0b', '')
    return Bin.rjust(p,'0') 

def Bin2State(Bin):
    state_index = int(Bin,2)
    state = ["0" for i in range(2**len(Bin))]
    state[state_index]="1"
    encoded_state = ''.join(state)
    return(encoded_state)

# score function for linear regression
def score(Bin,Y_train,X_train,Y_test,X_test,criterion = "MSE"):
    predictor_index = [i for i in range(len(Bin)) if Bin[i]=="1"]
    Z_train = X_train[:,predictor_index]
    Z_test = X_test[:,predictor_index]
    if len(predictor_index)==0:
        Z_train = np.ones((X_train.shape[0],1))
        Z_test = np.ones((X_test.shape[0],1))
    reg = LinearRegression().fit(Z_train, Y_train)
    if criterion == "MSE":
        Score = np.mean((reg.predict(Z_test)-Y_test)**2)
    elif criterion == "Rsquared":
        Score = np.mean((reg.predict(Z_test)-Y_test)**2)/np.mean((np.mean(Y_test)-Y_test)**2)
    else:
        print("used method not found! The choices are: MSE, Rsquared.")
    return(Score)

def Oracle(qc, nqubits, oracle_index):
    oracle_diag = np.ones(2**nqubits,)
    oracle_diag[oracle_index] = -1
    oracle_op = Operator(np.diag(oracle_diag))
    return(oracle_op)

def initialize_s(qc, qubits):
    """Apply a H-gate to 'qubits' in qc"""
    for q in qubits:
        qc.h(q)
    return qc

def Grover(oracle_index ,nqubits,iterations,backend):
    qc = QuantumCircuit(nqubits)
    qc = initialize_s(qc, [k for k in range(nqubits)])
    for iter in range(iterations):
        ##### apply oracle gate
        qc.append(Oracle(qc, nqubits, oracle_index),list(range(nqubits)))
        ##### apply applitude gate
        for qubit in range(nqubits):
            qc.h(qubit)
        # Apply transformation |00..0> -> |11..1> (X-gates)
        for qubit in range(nqubits):
            qc.x(qubit)
        # Do multi-controlled-Z gate
        qc.h(nqubits-1)
        qc.mct(list(range(nqubits-1)), nqubits-1)  # multi-controlled-toffoli
        qc.h(nqubits-1)
        # Apply transformation |11..1> -> |00..0>
        for qubit in range(nqubits):
            qc.x(qubit)
        # Apply transformation |00..0> -> |s>
        for qubit in range(nqubits):
            qc.h(qubit)
    qc.measure_all()
    transpiled_grover_circuit = transpile(qc, backend, optimization_level=1)
    job = backend.run(transpiled_grover_circuit,shots=128)
    job_monitor(job, interval=2)
    results = job.result()
    answer = results.get_counts()
    #time = job.time_per_step()['COMPLETED'] - job.time_per_step()['RUNNING']
    time = 999
    #new_Bin = max(answer,key=answer.get)
    new_Bins = answer.keys()#nlargest(16,answer,key=answer.get)
    return(new_Bins,time)

def BIC(Bin,Y_train,X_train):
    predictor_index = [i for i in range(len(Bin)) if Bin[i]=="1"]
    Z_train = X_train[:,predictor_index]
    if len(predictor_index)==0:
        Z_train = np.ones((X_train.shape[0],1))
    p = Z_train.shape[1]
    n = X_train.shape[0]
    reg = LinearRegression().fit(Z_train, Y_train)
    Y_est = reg.predict(Z_train)
    loss = np.sum((Y_est-Y_train)**2)
    BIC = (p+2)*np.log(n)+n*np.log(loss/n)
    return(BIC)

def BGS(Input):
    score= Input["score"]
    initial_index= Input["initial_index"]
    p = Input["StateNum"]
    state = ["1" for i in range(initial_index+1)]
    state.extend(["0" for i in range(2**p-initial_index-1)])
    selected_index = initial_index
    score_initial = score[initial_index]
    Selected_Scores =[]
    time_used = 0
    while int(score_initial) != 0:
        iter_time = int(0.25*np.pi*np.sqrt(2**p/(state_index+1) ) )+1
        oracles_bin = [list(bin(k).replace('0b', '').rjust(p,'0') ) for k in range(initial_index)]
        new_Bin,time = Grover(oracles_bin, p, iter_time, backend)
        time_used = time_used+time.total_seconds()
        state_index = int(new_Bin,2)
        score_new = score[state_index]
        if score_initial> score_new:
            selected_index = state_index
            score_initial = score_new
        Selected_Scores.append(score_initial)
    return([selected_index,Selected_Scores,time_used])

def BGS_linearReg(Input):
    Y_train = Input["Y_train"]
    X_train = Input["X_train"]
    Y_test = Input["Y_val"]
    X_test = Input["X_val"]
    criterion = Input["criterion"]
    initial_index = Input["initial_index"]
    p = X_train.shape[1]
    max_iter = p*10
    initial_Bin = ["0" for i in range(p)]
    for index in initial_index:
        initial_Bin[index]= '1'
    score_initial = score(initial_Bin,Y_train,X_train,Y_test,X_test,criterion)
    ScoresList = []
    Binlist = []
    for i in range(2**p):
        Bin = bin(i).replace('0b', '')
        Binlist.append(list(Bin.rjust(p,'0')))
        ScoresList.append(score(Bin.rjust(p,'0'),Y_train,X_train,Y_test,X_test,criterion))
    ScoresList = np.array(ScoresList)
    Selected_Scores =[]
    time_used = 0
    nums_iter = 0
    selected_index = int(''.join(initial_Bin),2)
    while sum(score_initial>ScoresList) != 0 and nums_iter< max_iter:
        nums_iter = nums_iter+1
        num_oracle = sum(score_initial>=ScoresList)
        iter_time = round(0.25*np.pi/np.arcsin(np.sqrt(num_oracle/2**p))-0.5 )
        OracleState_index = np.where(ScoresList<score_initial)[0]
        oracles_bin = [Binlist[k] for k in OracleState_index]
        new_Bin,time = Grover(OracleState_index,p,iter_time,backend)
        state_index = int(new_Bin,2)
        score_new = ScoresList[state_index]
        if score_initial > score_new:
            selected_index = state_index
            score_initial = score_new
        Selected_Scores.append(score_initial)
    return([selected_index,Selected_Scores,time_used])



def BGS_linearReg_BIC(Input):
    Y_train = Input["Y_train"]
    X_train = Input["X_train"]
    p = X_train.shape[1]
    # initialize
    score_candidate = []
    for initial_index in Input["initial_index"]:
        initial_Bin = ["0" for i in range(p)]
        for index in initial_index:
            initial_Bin[index]= '1'
        score_candidate.append(BIC(initial_Bin,Y_train,X_train))
    best_idx = np.argmin(np.array(score_candidate))
    initial_index = Input["initial_index"][best_idx]
    initial_Bin = ["0" for i in range(p)]
    for index in initial_index:
        initial_Bin[index]= '1'
    score_initial = BIC(initial_Bin,Y_train,X_train)

    ScoresList = []
    Binlist = []
    for i in range(2**p):
        Bin = bin(i).replace('0b', '')
        Binlist.append(list(Bin.rjust(p,'0')))
        ScoresList.append(BIC(Bin.rjust(p,'0'),Y_train,X_train))
    ScoresList = np.array(ScoresList)
    BSS_selected_idx = np.argmin(ScoresList)
    Selected_Scores =[]
    time_used = 0
    nums_iter = 0
    selected_index = int(''.join(initial_Bin),2)
 
    while score_initial != min(ScoresList) and nums_iter < 200:
        num_oracle = sum(score_initial>ScoresList)
        iter_time = round(0.25*np.pi/np.arcsin(np.sqrt(num_oracle/2**p))-0.5)
        OracleState_index = np.where(ScoresList<score_initial)[0]
        oracles_bin = [Binlist[k] for k in OracleState_index]
        new_Bins,time = Grover(OracleState_index,p,iter_time,backend)
        state_indexs = [int(new_Bin,2) for new_Bin in new_Bins]
        best_idx = np.argmin(ScoresList[state_indexs])
        state_index = state_indexs[best_idx]
        #new_Bin = new_Bins[best_idx]
        score_new = ScoresList[state_index]
        if score_initial > score_new:
            selected_index = state_index
            score_initial = score_new
        Selected_Scores.append(score_initial)
        nums_iter = nums_iter+iter_time*16
        # break the loop if 10 times of iteration get the same res
        if len(Selected_Scores)>=5 and min(Selected_Scores[-5:])==max(Selected_Scores[-5:]):
            break
    return([selected_index,Selected_Scores,nums_iter,sum(score_initial>=ScoresList),BSS_selected_idx])


def BGS_BIC(Input):
    BICs = Input["criterion"]
    p = int(np.log2(len(BICs)))
    m = sum(BICs<10**40)
    # initialize
    score_candidate = []
    for initial_index in Input["initial_index"]:
        score_candidate.append(BICs[initial_index])

    best_idx = np.argmin(np.array(score_candidate))

    initial_index = Input["initial_index"][best_idx]
    selected_index = initial_index
    score_initial = BICs[best_idx]
    Selected_Scores = [score_initial]
    nums_iter = 0
    while score_initial != min( BICs) and nums_iter < 200:
        num_oracle = sum(score_initial>BICs)
        iter_time = round(0.25*np.pi/np.arcsin(np.sqrt(num_oracle/2**p))-0.5)
        OracleState_index = np.where(BICs<score_initial)[0]
        new_Bins,time = Grover(OracleState_index,p,iter_time,backend)
        state_indexs = [int(new_Bin,2) for new_Bin in new_Bins]
        best_idx = np.argmin(BICs[state_indexs])
        state_index = state_indexs[best_idx]
        #new_Bin = new_Bins[best_idx]
        score_new = BICs[state_index]
        if score_initial > score_new:
            selected_index = state_index
            score_initial = score_new
        Selected_Scores.append(score_initial)
        nums_iter = nums_iter+iter_time*16
        # break the loop if 10 times of iteration get the same res
        if len(Selected_Scores)>=5 and min(Selected_Scores[-5:])==max(Selected_Scores[-5:]):
            break

    if selected_index>=m:
        selected_index = np.random.choice(m,1)[0]
    
    return([selected_index,Selected_Scores,nums_iter,sum(score_initial>=BICs)])



def BSS_linearReg(Input):
    Y_train = Input["Y_train"]
    X_train = Input["X_train"]
    Y_test = Input["Y_test"]
    X_test = Input["X_test"]
    initial_index = Input["initial_index"]
    p = X_train.shape[1]
    initial_Bin = ["0" for i in range(p)]
    for index in initial_index:
        initial_Bin[index]= '1'
    ScoresList = []
    Binlist = []
    for i in range(2**p):
        Bin = bin(i).replace('0b', '')
        Binlist.append(list(Bin.rjust(p,'0')))
        ScoresList.append(BIC(Bin.rjust(p,'0'),Y_train,X_train))
    ScoresList = np.array(ScoresList)
    best_idx = np.argmin(ScoresList)
    best_select = np.array([i == '1' for i in Binlist[best_idx]])+0
    return(best_select)

def BGS_evalute(data):
    Input = data["Input"]
    test = data["test"]
    res = data["res"]
    criterion = data["criterion"]
    Y_train = Input["Y_train"]
    X_train = Input["X_train"]
    Y_test = test["Y_test"]
    X_test = test["X_test"]
    p = X_train.shape[1]
    BGS_select = bin(res[0]).replace('0b','')
    BGS_score = score(BGS_select.rjust(p,'0'),Y_train,X_train,Y_test,X_test,criterion)
    return(BGS_score)


def BSS_evalute_old(data):
    Input = data["Input"]
    test = data["test"]
    criterion = data["criterion"]
    Y_test = test["Y_test"]
    X_test = test["X_test"]
    Y_train = Input["Y_train"]
    X_train = Input["X_train"]
    p = X_train.shape[1]
    ScoresList = []
    Binlist = []
    for i in range(2**p):
        Bin = bin(i).replace('0b', '')
        Binlist.append(Bin.rjust(p,'0'))
        ScoresList.append(BIC(Bin.rjust(p,'0'),Y_train,X_train))
    ScoresList = np.array(ScoresList)
    best_idx = np.argmin(ScoresList)
    BSS_select = Binlist[best_idx]
    BSS_score = score(BSS_select.rjust(p,'0'),Y_train,X_train,Y_test,X_test,criterion)
    return([BSS_select, BSS_score])

def BSS_evalute(data):
    Input = data["Input"]
    test = data["test"]
    criterion = data["criterion"]
    Y_test = test["Y_test"]
    X_test = test["X_test"]
    Y_train = Input["Y_train"]
    X_train = Input["X_train"]
    p = X_test.shape[1]
    best_idx = data["res"][4]
    Bin = bin(best_idx).replace('0b', '')
    BSS_select = Bin.rjust(p,'0')
    BSS_score = score(BSS_select.rjust(p,'0'),Y_train,X_train,Y_test,X_test,criterion)
    return([BSS_select,BSS_score])

def BGS_evalute(data):
    Input = data["Input"]
    test = data["test"]
    criterion = data["criterion"]
    Y_test = test["Y_test"]
    X_test = test["X_test"]
    Y_train = Input["Y_train"]
    X_train = Input["X_train"]
    p = X_test.shape[1]
    best_idx = data["res"][0]
    Bin = bin(best_idx).replace('0b', '')
    BGS_select = Bin.rjust(p,'0')
    BGS_score = score(BGS_select.rjust(p,'0'),Y_train,X_train,Y_test,X_test,criterion)
    return([BGS_select,BGS_score])