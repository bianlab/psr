# Complex-domain super-resolution imaging with distributed optimization

This repository contains the code for a complex-domain pixel super resolution method via distributed optimization. For more information, please contact Liheng Bian (bian at bit dot edu dot cn).

## About complex-domain pixel super resolution and reported method
Complex-domain imaging has emerged as a valuable technique for investigating weak-scattered samples. However, due to the detector’s pursuit of large pixel size for high throughput, the resolution limitation impedes its further development. In this work, we report a lensless on-chip complex-domain imaging system, together with a distributed-optimization-based pixel super-resolution technique (DO-PSR). The system employs a diffuser shifting to realize phase modulation and increases observation diversity. The corresponding DO-PSR technique derives an alternating projection operator and an enhancing neural network to tackle the measurement fidelity and statistical prior regularization subproblems. Extensive experiments show that the system outperforms the existing techniques with as much as 11dB on PSNR, and one-order-of-magnitude higher cell counting precision.


## Usage

Please clone this repository by Git or download the zip file firstly. 

### main_PSR

Run `main_PSR.m` file to obtain the pixel super resolution results of compared algorithms (Conventional-PSR, SR-SPAR, AS-PSR) and ours (DOTV-PSR and DONet-PSR).
 
### Note
This demo code runs under gaussian noise

For different noise level, undersampling ratio and reslution, please adjustment parameters to acquire better results.


## Platform

All the experiments are implemented using MATLAB 2019b with an Intel i7-9700 processor at 3.0GHz and 16GB RAM. 


## Requirements and Dependencies

Notice that PNP-PSR algorithms require Matconvnet (1.0-beta25), CUDA (10.0) and cuDNN.
