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
