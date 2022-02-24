import numpy as np
from qiskit import QuantumCircuit
from qiskit.extensions import UnitaryGate
from qiskit.quantum_info.operators import Operator

def a_dag(Lmax,i,verbose=False):
    """
    TO DO: implement equation 7 and use a binary qubit mapping
    will be the backbone of the algorithm
    Args:
        Lmax (int): _description_
    """
    num_qubits=len(format(Lmax,'b')) #we find the number of qubits needed from the maximum number of binary digits we need
    qc=QuantumCircuit(num_qubits) #it is easiest to create a circuit that represents this operator, then use qiskit to export a unitary matrix 
    sig_01=Operator([[0,1],[0,0]]) #sigma minus
    sig_10=Operator([[0,0],[1,0]]) #sigma plus
    sig_00=Operator([[1,0],[0,0]])
    sig_11=Operator([[0,0],[0,1]])
    
    a_dag=0
    for l in range(1,Lmax): #NB start at 1 to prevent zero error
        l_b_string=format(l,'b') #converts the l value to a binary string
        l_b_string=b_string.zfill(num_qubits) #pads out this number with enough zeros to fill in values for each qubit
        
        l_minus_b_string=format(l-1,'b') #converts the l-1 value to a binary string
        l_minus_b_string=l_minus_b_string.zfill(num_qubits) #pads out this number with enough zeros to fill in values for each qubit
        
        qc=QuantumCircuit(num_qubits)

        for index,digit in enumerate(l_b_string):
            if digit==l_minus_b_string[index]==0:
                qc.unitary(sig_00)
            elif digit==l_minus_b_string[index]==1:
                qc.unitary(sig_11)
            elif digit==0 and l_minus_b_string[index]==1:
                qc.unitary(sig_01)
            elif digit==1 and l_minus_b_string[index]==1:
                qc.unitary(sig_10)

        if verbose==True:
            qc.draw()   

        vib_op=
        a_dag+=vib_op

        continue
    return a_dag

def b_dag(OhmA,OhmB,a,delta,i):
    """
    COMPLETE (hopefully)
    Takes the Ohm matricies of the two PES and uses them to transform the inital
    ladder operators to ones for the new PES. NB must be repeated for each element of the omega vectors
    This is equivelent to equation S10
    N.B can be called once and h.c can be taken to find the other operatore

    Args:
        OhmA (matrix): property of ground state PES
        OhmB (matrix): property of excited state PES
        a (vector): ground state anhilation operator
        delta (vector): displacement vector

    Returns:
        b (vector): excited state creation operator
    """
    a=a_dag(Lmax,i).getH()
    J=OhmB @ S @ np.linalg.inv(OhmA)
    ans=0.5*(J-np.linalg.inv(np.transpose(J))) @ a
    + 0.5*(J+np.linalg.inv(np.transpose(J))) @ a.getH()
    +(1/np.sqrt(2))*delta
    return ans

def HB(omegaB):
    """
    TO DO: work out how to get inital size of matrix
    This constructs the hamiltonian for PES B. The technique use is the second one
    mentioned in supplemntary material part S1, just below equation S10
    Args:
        omegaB (_type_): _description_

    Returns:
        _type_: _description_
    """    
    HB=np.zeros()
    for i in range(len(omegaB)):
        b_dagi=b_dag(OhmA,OhmB,a,delta,i)
        HB+=omegaB[i]*(b_dagi @ b_dagi.getH()+0.5)
    return HB