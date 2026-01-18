# from .helpers import generate_spike_train, STTC_formula, STTC_single
# try:
#     import torch
#     from .gpu_main import STTC_gpu, null_gpu, add_dt_gpu, shift_signal_gpu  # Torch version with GPU functionality
#     print("Loaded with PyTorch backend.")
# except ImportError:
#     from .cpu_main import STTC, null, add_dt, shift_signal  # CPU version
#     print("Loaded.")



from .helpers import generate_spike_train, STTC_formula, STTC_single
from .cpu_main import STTC, null, add_dt, shift_signal

try:
    import torch
    from .gpu_main import STTC_gpu, null_gpu, add_dt_gpu, shift_signal_gpu
    print("GPU functions available on PyTorch backend.")
except ImportError:
    STTC_gpu = None
    null_gpu = None
    add_dt_gpu = None
    shift_signal_gpu = None
    print("PyTorch not found. GPU functions are unavailable.")
