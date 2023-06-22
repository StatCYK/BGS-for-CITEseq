# Importing standard Qiskit libraries and configuring account
from qiskit import QuantumCircuit, execute, Aer, IBMQ
from qiskit.compiler import transpile, assemble
from qiskit.tools.jupyter import *
from qiskit.visualization import *
from qiskit import IBMQ, QuantumRegister, ClassicalRegister, execute, QuantumCircuit
# Loading your IBM Q account(s)
import pylab
import numpy as np
from qiskit import BasicAer
from qiskit.tools.visualization import plot_histogram
import time
from qiskit.tools.monitor import job_monitor
from datetime import datetime
from sklearn.linear_model import LinearRegression


######## backend setting ###############

IBMQ.load_account()
#provider = IBMQ.get_provider(hub = "ibm-q-research-2", group='uni-georgia-1', project='main')
provider = IBMQ.get_provider(hub = "ibm-q", group='open', project='main')
backend = provider.get_backend('ibmq_qasm_simulator')#('ibmq_manila')
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
def score(Bin,Y_train,X_train,Y_test,X_test):
    predictor_index = [i for i in range(len(Bin)) if Bin[i]=="1"]
    Z_train = X_train[:,predictor_index]
    Z_test = X_test[:,predictor_index]
    if len(predictor_index)==0:
        Z_train = np.ones((X_train.shape[0],1))
        Z_test = np.ones((X_test.shape[0],1))
    reg = LinearRegression().fit(Z_train, Y_train)
    Score = np.mean((reg.predict(Z_test)-Y_test)**2)
    return(Score)

def Oracle(qc, oracle):
    nqubits = len(oracle)
    for qubit in range(nqubits):
        if oracle[qubit] =='0':
            qc.x(qubit)
    qc.h(nqubits-1)
    qc.mct(list(range(nqubits-1)), nqubits-1)  # multi-controlled-toffoli
    qc.h(nqubits-1)
    for qubit in range(nqubits):
        if oracle[qubit] =='0':
            qc.x(qubit)
    return(qc)

def initialize_s(qc, qubits):
    """Apply a H-gate to 'qubits' in qc"""
    for q in qubits:
        qc.h(q)
    return qc

def Grover(oracles_bin ,nqubits,iterations,backend):
    qc = QuantumCircuit(nqubits)
    qc = initialize_s(qc, [k for k in range(nqubits)])
    for iter in range(iterations):
        ##### apply oracle gate
        for oracle in oracles_bin:
            oracle.reverse()
            qc = Oracle(qc,oracle)
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
    transpiled_grover_circuit = transpile(qc, backend, optimization_level=3)
    job = backend.run(transpiled_grover_circuit,shots=1)
    job_monitor(job, interval=2)
    results = job.result()
    answer = results.get_counts()
    time = job.time_per_step()['COMPLETED'] - job.time_per_step()['RUNNING']
    new_Bin = max(answer,key=answer.get)
    return(new_Bin,time)

def QAS(Input):
    score= Input["score"]
    learn_rate=Input["learn_rate"]
    initial_index= Input["initial_index"]
    m = 1
    p = Input["StateNum"]
    
    state = ["1" for i in range(initial_index+1)]
    state.extend(["0" for i in range(2**p-initial_index-1)])
    selected_index = initial_index
    score_initial = score[initial_index]
    
    Selected_Scores =[]
    
    time_used = 0
    while (m < 15*learn_rate*np.log(p)) and int(score_initial) != 0:
        iter_time = int(0.25*np.pi*learn_rate**(-m/2))+1
        m = m+1
        oracles_bin = [list(bin(k).replace('0b', '').rjust(p,'0') ) for k in range(initial_index)]
        new_Bin,time = Grover(oracles_bin ,p,iter_time,backend)
        time_used = time_used+time.total_seconds()
        state_index = int(new_Bin,2)
        score_new = score[state_index]
        if score_initial> score_new:
            OracleState_index = (score < score_new)+0
            selected_index = state_index
            score_initial = score_new
            m = 1
        Selected_Scores.append(score_initial)
    print("propose state is %d"%selected_index)
    return([selected_index,Selected_Scores,time_used])


def QAS_linearReg(Input):
    Y_train= Input["Y_train"]
    X_train= Input["X_train"]
    Y_test= Input["Y_test"]
    X_test= Input["X_test"]
    initial_index= Input["initial_index"]
    learn_rate=Input["learn_rate"]
    m = 1
    Grover_iters = 1
    p = X_train.shape[1]
    
    initial_Bin = "0"*p
    for index in initial_index:
        initial_Bin[index]=1
    score_initial = score(initial_Bin,Y_train,X_train,Y_test,X_test)
    
    ScoresList = []
    Binlist = []
    for i in range(2**p):
        Bin = bin(i).replace('0b', '')
        Binlist.append(list(Bin.rjust(p,'0')))
        ScoresList.append(score(Bin.rjust(p,'0'),Y_train,X_train,Y_test,X_test))
    ScoresList = np.array(ScoresList)
    Selected_Scores =[]

    time_used = 0
    while (m < 15*learn_rate*np.log(p)) and sum(score_initial>ScoresList) != 0:
        iter_time = int(0.25*np.pi*learn_rate**(-m/2))+1
        m = m+1
        OracleState_index = np.where(ScoresList<score_initial)[0]
        oracles_bin = [Binlist[k] for k in OracleState_index]

        new_Bin,time = Grover(oracles_bin,p,iter_time,backend)
        time_used = time_used+time.total_seconds()
        state_index = int(new_Bin,2)
        score_new = ScoresList[state_index]
        if score_initial> score_new:
            selected_index = state_index
            score_initial = score_new
            m = 1
        Selected_Scores.append(score_initial)
    print("propose state is %d"%selected_index)
    return([selected_index,Selected_Scores,time_used])