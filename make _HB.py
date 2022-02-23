import numpy as np

def a_dag(Lmax,i):
    """
    TO DO: implement equation 7 and use a binary qubit mapping
    will be the backbone of the algorithm
    Args:
        Lmax (int): _description_
    """    
    a_dag=0
    for l in range(Lmax):
        continue
    return a_dag

def b_dag(OhmA,OhmB,a,delta,i):
    """
    COMPLETE (hopefully)
    Takes the Ohm matricies of the two PES and uses them to transform the inital
    ladder operators to ones for the new PES. NB must be repeated for each element of the omega vectors
    
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