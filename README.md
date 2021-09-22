# Plug-and-play optimization for complex-domain pixel super resolution

This repository contains the code for a complex-domain pixel super resolution method via Plug-and-play optimization. For more information, please contact Liheng Bian (bian at bit dot edu dot cn).

## About complex-domain pixel super resolution and reported method
In order to increase signal-to-noise ratio in measurement, most imaging detectors sacrifice resolution to increase pixel size in confined area. Although the pixel super-resolution technique (PSR) enables resolution enhancement in such as digital holographic imaging, it suffers from unsatisfied reconstruction quality. In this work, we report a high-fidelity plug-and-play optimization method for PSR phase retrieval, termed as PNP-PSR. It decomposes PSR reconstruction into independent sub-problems based on the generalized alternating projection framework. An alternating projection operator and an enhancing neural network are derived to tackle the measurement fidelity and statistical prior regularization, respectively. In this way, PNP-PSR incorporates the advantages of individual operators, achieving both high efficiency and noise robustness. We compare PNP-PSR with the existing PSR phase retrieval algorithms with a series of simulations and experiments, and PNP-PSR outperforms the existing algorithms with as much as 11dB on PSNR. The enhanced imaging fidelity enables one-order-of-magnitude higher cell counting precision.


## Usage

Please clone this repository by Git or download the zip file firstly. 

### main_PSR

Run `main_PSR.m` file to obtain the pixel super resolution results of compared algorithms (Conventional-PSR, SR-SPAR, AS-PSR) and ours (TV-PSR and PNP-PSR).
 
### Note
This demo code runs under gaussian noise

For different noise level, undersampling ratio and reslution, please adjustment parameters to acquire better results.


## Platform

All the experiments are implemented using MATLAB 2019b with an Intel i7-9700 processor at 3.0GHz and 16GB RAM. 


## Requirements and Dependencies

Notice that PNP-PSR algorithms require Matconvnet (1.0-beta25), CUDA (10.0) and cuDNN.
