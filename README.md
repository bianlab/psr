# Plug-and-play pixel super-resolution phase retrieval for digital holography

This repository contains the code for a pixel super-resolution methods. For more information, please contact Liheng Bian (bian at bit dot edu dot cn).

## About pixel super-resolution and the reported methods
In order to increase signal-to-noise ratio in optical imaging, most detectors sacrifice resolution to increase pixel size in a confined area, which impedes further development of high throughput holographic imaging. Although the pixel super-resolution technique (PSR) enables resolution enhancement, it suffers from the trade-off between reconstruction quality and super-resolution ratio. In this work, we report a high-fidelity PSR phase retrieval method with plug-and-play optimization, termed PNP-PSR. It decomposes PSR reconstruction into independent sub-problems based on generalized alternating projection framework. An alternating projection operator and an enhancing neural network are derived to tackle the measurement fidelity and statistical prior regularization, respectively. PNP-PSR incorporates the advantages of individual operators, achieving both high efficiency and noise robustness. Extensive experiments show that PNP-PSR outperforms the existing techniques in both resolution enhancement and noise suppression.


## Usage

Please clone this repository by Git or download the zip file firstly. 

### main_PSR

Run `main_PSR.m` file to obtain the pixel super-resolution results of the compared algorithms (Conv-PSR, SR-SPAR, AS-PSR) and ours (PNPTV-PSR and PNPNet-PSR).
 
### Note
This demo code runs under gaussian noise

For different noise level, undersampling ratio and reslution, please adjustment parameters to acquire better results.


## Platform

All the experiments are implemented using MATLAB 2019b with an Intel i7-9700 processor at 3.0GHz and 16GB RAM. 


## Requirements and Dependencies

Notice that PNP-PSR algorithms require Matconvnet (1.0-beta25), CUDA (10.0) and cuDNN.
