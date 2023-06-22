# PNAS Bisection Grover’s Search Algorithm  and Its Application in Analyzing CITE-seq Data

<img src="circuit_example.png" width="1080" height="220" />

## Requirements
qiskit=0.36.2

numpy

itertools

functools

## Run experiments

### experiment on IBM Quantum system

Upload `QAS_IBM_tunePar.ipynb` ,`QAS_IBM_linearReg.ipynb` and `utility.py` to https://lab.quantum-computing.ibm.com/. You need to change the provider information in the code `utility.py`.

For example, setting `IBMQ.get_provider(hub = "ibm-q-research-2", group='uni-georgia-1', project='main')`

You can change the Quantum system service. For example, setting `provider.get_backend('ibmq_manila')`.

Some services do not provide public access. A recommended approach to test the code is to use the 'ibmq_qasm_simulator' service. Be noticed that 'ibmq_qasm_simulator' isn't real quantum computer thought.

A list of different Quantum systems is available in https://quantum-computing.ibm.com/services?services=systems.

The `QAS_IBM_tunePar.ipynb` is a demo for parameter selection and `QAS_IBM_linearReg.ipynb` is a demo for linear regression experiments.

### simulated experiment on classical computer

Folder `simu experiments_classic computer/python/` contains the python code of the proposed QAS, generating simulated data and the real experiments.

Folder `simu experiments_classic computer/R/` contains the code of the compared methods and the results.

## Acknowledgement

Detailed tutorials of constructing Quantum circuit and implementation of Quantum algorithms can be found in https://qiskit.org/documentation/tutorials.html

## Contact us

Website: https://bigdata.uga.edu/
