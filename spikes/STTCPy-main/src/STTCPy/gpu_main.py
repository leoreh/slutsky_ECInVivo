import torch
from .helpers import STTC_formula

def add_dt_gpu(signal, dt=0, device='cuda'):
    '''
    Create tiled signal. Make 0->1 in the surrounding (+/- dt) values around each spike (1)
    This is executed using a 1-d convolution with a rectangle window across each row in the signal.

    For example:
    for dt = 1,
    signal = [1 0 0 0 1 0 0 0 0 1 1 0 0]
    output = [1 1 0 1 1 1 0 0 1 1 1 1 0]

    Inputs
    - signal : (torch.Tensor). Can be 1-d or 2-d array of signals. Output is the same type of Tensor as input. 
                Only put in signals in the format (n_cells, n_samples in time)
    - dt : (int, default = 0)

    Output
    - convolved_signal (torch.Tensor) in the same type and shape as input Tensor.
    '''
    if not torch.is_tensor(signal):
        raise ValueError("Input must be a PyTorch tensor.")

    
    if dt == 0:
        return signal
    
    change_ndim = False
    if signal.ndim==1:
        signal = signal.view(1,-1)
        change_ndim = True


    # Pad the signal
    padded_signal = torch.nn.functional.pad(signal, (dt, dt), value=0)

    # Create a window with ones of size (2*dt+1)
    N,T = signal.shape
    window = torch.ones(N, 2 * dt + 1, device=device, dtype=signal.dtype)

    # Convolve the padded signal with the window
    convolved_signal = torch.nn.functional.conv1d(padded_signal, 
                                                  window.unsqueeze(1), padding=dt, groups=N)

    if change_ndim:
        return ((convolved_signal > 0)*1)[0, dt:-dt].type(signal.dtype)
    return ((convolved_signal > 0)*1)[:, dt:-dt].type(signal.dtype)



def STTC_gpu(signals, dt=0, ret_all = False):
    '''
    Calculates pairwise STTC of all rows in the signals array.

    Inputs
    - signals : (torch.Tensor). 2-d array of signals.
                Only put in signals in the format (n_cells, n_samples in time)
    - dt : (int, default = 0)
    - ret_all : (boolean, default = False) if set to True, returns dictionary with the intermediate calculations of 
                                            TA, TB, PA and PB along with STTC

    Outputs
    - STTC : (torch.Tensor) retruns a (n_cells x n_cells) matrix with each cell giving the paired STTC calculated at set dt.
    - TA, TB, PA, PB : intermediate calculated values that go into the STTC.
    '''
    if not torch.is_tensor(signals):
        raise ValueError("Input must be a PyTorch tensor.")

    if signals.dtype!=torch.float:
        raise ValueError("Input must be a float tensor, dtype=torch.float, torch.float32 or torch.float64
        
    N,L = signals.shape
    dt_signals = add_dt(signals, dt=dt)
    NB = torch.tile(signals.sum(axis=1), (N,1))
    TB = torch.tile(dt_signals.mean(axis=1), (N,1))
    TA = TB.T
    
    DAB = signals @ dt_signals.T

    PA = DAB/(NB.T)
    PB = DAB.T/(NB)
    
    STTC = STTC_formula(TA, TB, PA, PB)
    
    if ret_all:
        return dict(zip(['STTC', 'TA', 'TB', 'PA', 'PB'], [STTC, TA, TB, PA, PB]))
    
    return STTC

def shift_signal_gpu(signals, shifts=None, shuffle=False):
    '''
    shuffle or circular shifts signal according to the second dimension
    '''
    if shuffle:
        return signals[:, torch.randperm(signals.size(1))]
    return torch.roll(signals, shifts=shifts, dims=1)

def null_gpu(signals, n_shifts, dt=0, shuffle=False):
    '''
    Calculates null where each row's (cell) STTC is calculated against a circular-shifted array of every other cell. 
    This is done for (ns = 1 to n_shifts + 1) values

    If shuffle=True, instead of circular shifting, the signals are instead randomly shuffled 'n_shift' times.

    The output is shaped as (ns, n_cell, n_cell).
    Each cell can be calculated as STTC(A, torch.roll(B, shifts=ns)) (if shuffle=False)

    Inputs
    - signals : (torch.Tensor). 2-d array of signals.
                Only put in signals in the format (n_cells, n_samples in time)
    - n_shifts : (int). Maximal value of shifts to be calculated.
    - dt : (int, default = 0)
    - ret_all : (boolean, default = False) if set to True, returns dictionary with the intermediate calculations of 
                                            TA, TB, PA and PB along with STTC

    Outputs
    - null : (torch.Tensor) retruns a (ns x n_cells x n_cells) matrix with each cell giving the 
                            paired STTC calculated at set dt with the j-index being shifted by ns.

    '''
    if not torch.is_tensor(signals):
        raise ValueError("Input must be a PyTorch tensor.")

    if signals.dtype!=torch.float:
        raise ValueError("Input must be a float tensor, dtype=torch.float, torch.float32 or torch.float64.")
        
    N,L = signals.shape

    null = []

    dt_signals = add_dt(signals, dt=dt)
    NB = torch.tile(signals.sum(axis=1), (N,1))
    TA = torch.tile(dt_signals.mean(axis=1), (N,1)).T

    for ns in range(n_shifts):
        shifted = shift_signal(signals, shifts=ns+1, shuffle=shuffle)
        dt_shifted = add_dt(shifted, dt=dt)
        TB = torch.tile(dt_shifted.mean(axis=1), (N,1))

        DAB = signals @ dt_shifted.T
        DBA = dt_signals @ shifted.T

        PA = DAB/NB.T
        PB = DBA/NB

        null.append(STTC_formula(TA, TB, PA, PB))
        
    null = torch.stack(null, dim=0)
    return null


