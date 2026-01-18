# **STTCPy Documentation**

## **Overview**
`STTCPy` is a Python package designed to efficiently compute the Spike Time Tiling Coefficient (STTC), a metric introduced by Cutts and Eglen et al. to assess synchrony between neuronal spike trains ([Cutts and Eglen, 2014](#citation)). While existing STTC implementations in both Python and MATLAB are limited to one-to-one calculations (i.e., STTC between a single pair of signals), STTCPy vectorizes these computations, enabling much faster processing of large datasets.

This approach significantly accelerates STTC calculations, making it suitable for high-volume data analysis. The package also includes functionality for calculating null distributions for significance testing and offers GPU acceleration using `PyTorch` backend for even greater performance on large-scale analyses.

Key features:
- Efficient, vectorized STTC computation for large datasets.
- Scalable GPU-based functionality with PyTorch and CUDA support.
- Helper functions for pre-processing, analysis, and significance testing via null calculations.

---

## **Installation**

### Requirements:
- Python 3.8 or newer.
- Dependencies: `numpy`, `scipy`
 
**Note**: To enable GPU functionality, PyTorch must be installed with the appropriate CUDA version for your system. For installation instructions, visit the [PyTorch installation guide](https://pytorch.org/get-started/locally/).


### Installation:
Install the package via pip:
```bash
pip install STTCPy
```

Alternatively, clone the repository and install manually:
```bash
git clone https://github.com/aditiaravind/STTCPy.git
cd STTCPy
pip install .
```

---

## **Usage**

### Import the Package:
```python
from STTCPy import STTC, null
```


### Using CPU-only version:
#### Compute STTC:
```python
from STTCPy import STTC, null, generate_spike_train

N = 3 # Number of spike trains
L = 1000 # Length of each spike train
sparsity = 0.2 # Density of spikes within the train

spikes = generate_spike_train(N, L, sparsity)

dt = 2  # Number of frames to pad when calculating STTC

spikes_STTC = STTC(spikes, dt)
```

You can then plot the STTC values as follows
```python
plt.imshow(spikes_STTC, vmax=1, vmin=-1, cmap='bwr')
plt.colorbar()
plt.xlabel('STTC Value')
plt.ylabel('STTC Value')
plt.title(f'STTC across {N} spike trains')
plt.show()
```



#### Null Model for Statistical Testing:
Calculate STTC when the signals are shuffled or shifted 'n' times. Average this to compare against the actualy STTC values to check for significance. 

```python
n_shifts = 5

spikes_null = null(spikes, n_shifts, dt)
print(f"Null STTC Shape: {spikes_null.shape}")

# Output : (5, 10, 10)
```



### Using GPU:

Here, you have to have torch installed. 

```python
import torch
from STTCPy import STTC_gpu as STTC, null_gpu as null

N = 3 # Number of spike trains
L = 1000 $ Length of each spike train
sparsity = 0.2 # Density of spikes within the train

spikes_numpy = generate_spike_train(N, L, sparsity)
spikes = torch.tensor(spikes_numpy, device='cuda', dtype=torch.float32) #Specify dtype as float for best results.

dt = 2  # Number of frames to pad when calculating STTC
n_shifts = 5

spikes_STTC = STTC(spikes, dt)
null_STTC = null(spikes, n_shifts, dt)

```
For full documentation visit the [Github repository](https://github.com/aditiaravind/STTCPy) 

---

## **Acknowledgments**
I would like to acknowledge the following paper for introducing the Spike Time Tiling Coefficient (STTC) metric:

Cutts CS, Eglen SJ. Detecting pairwise correlations in spike trains: an objective comparison of methods and application to the study of retinal waves. *J Neurosci*. 2014 Oct 22;34(43):14288-303. https://doi.org/10.1523/jneurosci.2767-14.2014. PMID: 25339742; PMCID: PMC4205553.

---

## **License**  
STTCPy is distributed under the GNU General Public License v3 (GPLv3). For more details, refer to the [GNU GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.en.html).  


---

## **Contributing**

I welcome contributions! Here’s how you can help:
1. Fork the repository.
2. Create a new branch for your feature or bug fix:
   ```bash
   git checkout -b feature-name
   ```
3. Commit your changes and push to your branch:
   ```bash
   git commit -m "Description of changes"
   git push origin feature-name
   ```
4. Submit a pull request with a detailed description of your changes.

For issues or feature requests, please use the [GitHub issues](https://github.com/aditiaravind/STTCPy/issues) page.

---
