import numpy as np
from .cpu_main import add_dt
def generate_spike_train(N, T, sparsity):
    '''
    Creates an array of spike trains of size (N,T) at a given spike density (sparsity, or total spikes/length of signal)

    returns spike_train (np.array)
    '''
    
    # Check if sparsity is between 0 and 1
    if not (0 <= sparsity <= 1):
        raise ValueError("Spike density must be between 0 and 1.")
    if sparsity == 0:
        return np.zeros((N, T), dtype=int)
    elif sparsity == 1:
        return np.ones((N, T), dtype=int)
        
    num_non_zero = int(N * T * sparsity)
    indices = np.random.choice(N * T, num_non_zero, replace=False)
    spike_train = np.zeros((N, T), dtype=int)
    spike_train.flat[indices] = 1

    return spike_train

def spike_times_to_train():
    pass
    

def STTC_formula(TA, TB, PA, PB):
    '''
    Formula of STTC from the Cutts and Eglen (2014) Paper
    '''
    x = (PA - TB)/(1 - (PA*TB))
    y = (PB - TA)/(1 - (PB*TA))
    STTC = 0.5*(x + y)
    return STTC


def STTC_single(A, B, dt=0, ret_all=True):
    '''
    STTC calculation across two input signals without any loops/matrix. I call it the 'Direct Method'
    '''

    if A.ndim>1 or B.ndim>1:
        raise ValueError("Input must be 1-d signals.")

    if dt == 0:
        TA, tileA = A.mean(), A
        TB, tileB = B.mean(), B
    else:
        tileA = add_dt(A, dt=dt)
        TA = tileA.mean()

        tileB = add_dt(B, dt=dt)
        TB = tileB.mean()

    n_spikes_A = A.sum()
    n_spikes_B = B.sum()

    if n_spikes_A == 0:
        PA = 0
    else:
        PA = A.dot(tileB).sum()/n_spikes_A

    if n_spikes_B == 0:
        PB = 0
    else:
        PB = B.dot(tileA).sum()/n_spikes_B

    STTC = STTC_formula(TA, TB, PA, PB)
    
    if ret_all:
        return dict(zip(['STTC', 'TA', 'TB', 'PA', 'PB'], [STTC, TA, TB, PA, PB]))
    else:
        return STTC