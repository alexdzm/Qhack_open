import numpy as np
import math
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

    sig_01=np.array([[0,1],[0,0]])#sigma minus
    sig_10=np.array([[0,0],[1,0]])#sigma plus
    sig_00=np.array([[1,0],[0,0]])
    sig_11=np.array([[0,0],[0,1]])
    
    
    for l in range(1,Lmax): #NB start at 1 to prevent zero error
        l_b_string=format(l,'b') #converts the l value to a binary string
        l_b_string=l_b_string.zfill(num_qubits) #pads out this number with enough zeros to fill in values for each qubit
        
        l_minus_b_string=format(l-1,'b') #converts the l-1 value to a binary string
        l_minus_b_string=l_minus_b_string.zfill(num_qubits) #pads out this number with enough zeros to fill in values for each qubit

        for index,digit in enumerate(l_b_string):
            if index==0:
                
                if digit=='0' and l_minus_b_string[index]=='0':
                    a_dag=sig_00
                elif digit=='1' and l_minus_b_string[index]=='1':
                    a_dag=sig_11
                elif digit=='0' and l_minus_b_string[index]=='1':
                    a_dag=sig_01
                elif digit=='1' and l_minus_b_string[index]=='1':
                    a_dag=sig_10
            else:
                if digit=='0' and l_minus_b_string[index]=='0':
                    a_dag=np.tensordot(a_dag,sig_00,axes=0)
                elif digit=='1' and l_minus_b_string[index]=='1':
                    a_dag=np.tensordot(a_dag,sig_11,axes=0)
                elif digit=='0' and l_minus_b_string[index]=='1':
                    a_dag=np.tensordot(a_dag,sig_01,axes=0)
                elif digit=='1' and l_minus_b_string[index]=='1':
                    a_dag=np.tensordot(a_dag,sig_10,axes=0)
                    
            
        if l==1:
            a_dag_tot=a_dag
        else:
            a_dag_tot+=a_dag

    if verbose==True:
            print(a_dag_tot)
    return a_dag_tot

def b_dag(S,OhmA,OhmB,delta,Lmax,i):
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
    a=a_dag(Lmax,i).conj().T
    J=np.matmul(OhmB,np.matmul(S,np.linalg.inv(OhmA)))
    ans=0.5*np.matmul((J-np.linalg.inv(np.transpose(J))),a)
    + 0.5*np.matmul((J+np.linalg.inv(np.transpose(J))),a.conj().T)
    +(1/math.sqrt(2))*np.array(delta)
    return ans

def HB(S,omegaB,OhmA,OhmB,delta,Lmax):
    """
    TO DO: why does this only work for Lmax= powers of 2?
    
    This constructs the hamiltonian for PES B. The technique use is the second one
    mentioned in supplemntary material part S1, just below equation S10
    Args:
        omegaB (_type_): _description_

    Returns:
        _type_: _description_
    """    

    for i in range(len(omegaB)):
        b_dag_i=b_dag(S,OhmA,OhmB,delta,Lmax,i)
        if i==0:
            HB=omegaB[i]*(np.matmul(b_dag_i,b_dag_i.conj().T+0.5))
        else:
            HB+=omegaB[i]*(np.matmul(b_dag_i,b_dag_i.conj().T+0.5))
    return HB