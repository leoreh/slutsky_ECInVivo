# Welcome to the Slutsky lab repository!<br>

This repository involves MATLAB© code for analyzing neurophysiological signals, mainly tetrode recordings from the rodent hippocampus. However, most of the code can be easily applied to other brain regions and/or recording modalities.<br>

Everything here is continuously developed and improved by the lab. Accordingly, most functions and pipelines (wrappers) are far from polished. If you have any questions please feel free to contact heimleore@gmail.com.

### Dependencies
Many other repositories served as invaluable references for the code presented here, including:
1. https://github.com/michael-zugaro/FMAToolbox
2. https://github.com/open-ephys/analysis-tools
3. https://github.com/buzsakilab/buzcode
4. https://github.com/petersenpeter/CellExplorer
5. https://github.com/MouseLand/Kilosort
6. https://github.com/zekebarger/AccuSleep
7. https://github.com/IoSR-Surrey/MatlabToolbox

Currently, the only external dependencies for some of the code here is CellExplorer (by P. Peterson from the Buzsáki lab) and the Open Ephys analysis tools.

### Data format
Our recording rigs consist of hardware by Tucker-Davis Technologies, Intan Technologies, and A-M systems.<br>
Typically, we first convert the raw recordings (TDT tanks or winWCP files) to flat binaries before further processing.

We adopted the file organization standard from György Buzsáki's lab (NYU). Each recording is in its own directory *basepath* and the name of the directory *basename* is used as the prefix for all other files. For example, if *'/user/../m1_d2'* is the *basepath* input to the function *spktimesWh.m*, the function will search the folder *'/user/../m1_d2'* for the data file *'m1_d2.dat'* and save the output as *m1_d2.spktimes.mat*.<br>
Typically, *basename* is the name of the mouse followed by the date-time string in the format *'yyMMdd_HHmmss'*. For exmaple, *'lh96_211202_070500'*.

### General pipeline

See FunctionList.m for example calls of the most frequently used functions.

##### File conversion
Convert the raw data files to flat binaries with *tdt2dat.m* (for TDT tanks) or *preprocOE.m* (for Open Epyhs), or to .mat files with *getLFP.m* (for .wcp or .abf files).

##### Spike sorting
"Whiten" the raw data by common referencing adjusted to the cross correlation between channels (adapted from Kilosort). Extract the spikes from the whitened data with *spktimesWh.m*. Convert the spikes to neurosuite files (fet, res, and spk) with *spktimes2ns.m*.<br>
On a linux distribution (or Windows wsl) - cluster the spikes with klustakwik version 1 and manually inspect the clusters with klusters (neurosuite).<br>
In Matalb - Use *fixSpkAndRes.m* to realign spikes in a cluster after merging. Use *cleanCluByFet.m* to remove spikes from a cluster based on their distance in feature space.

##### Spikes
Load and visualize the clusters with *Cell Explorer*. Use *spktimesMetrics.m* and *spkwvMetrics.m* to analyze the temporal firing patterns and waveform parameters of the clusters, respectively. These metrices will be added to the *cell_metrics* struct of Cell Explorer so they can be visualized in their GUI.<br>
Note: the .spk files contain the waveform from the whitened data. the *spikes* struct from Cell Explorer contains the waveforms from the raw data (detrended).

##### LFP
Create an LFP flat binary from the raw data with *LFPfromDat.m*. Analyze burst-suppression (*getBS.m*), inter-ictal spikes (*getIIS.m*), ripples (*getRipples.m*), etc.<br> Analyze sleep states with a neural network from AccuSleep.
