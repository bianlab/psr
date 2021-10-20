% Comparison of complex-domain pixel super-resolution (CD-PSR) algorithms. 
% Algotithm lists: 1. Conventional PSR algorithm (Conv-PSR) [1]
%                  2. Super-resolution sparse phase-amplitude retrieval(SR-SPAR) [2]
%                  3. Adaptive smoothing PSR algorithm (AS-PSR) [3]
%                  4. Distributed-optimization Total-Variation PSR algorithm (DOTV-PSR) [4]
%                  5. Distributed optimization Denoising network PSR algorithm (DONet-PSR) [4]

% Reference
%   [1] Hu X, Li S, Wu Y. Resolution-enhanced subpixel phase retrieval 
%       method[J]. Applied Optics, 2008, 47(32): 6079-6087.
%   [2] Katkovnik V, Shevkunov I, Petrov N V, et al. Computational 
%       super-resolution phase retrieval from multiple phase-coded 
%       diffraction patterns: simulation study and experiments[J]. 
%       Optica, 2017, 4(7): 786-794.
%   [3] Gao Y, Cao L. High-fidelity pixel-super-resolved complex field 
%       reconstruction via adaptive smoothing[J]. Optics Letters, 
%       2020, 45(24): 6807-6810.
%   [4] Chang X, Bian L, Jiang S, et al. Complex-domain super-resolution 
%       imaging with distributed optimization[J]. arXiv preprint 
%       arXiv:2105.14746, 2021.

% Contact
%   Xuyang Chang, Beijing Institute of Technology (BIT), chang@bit.edu.cn
%   Liheng Bian, Beijing Institute of Technology (BIT), bian@bit.edu.cn

clc;clear;close all;
vl_setupnn;

%%  [1] Setup parameters
%   [1.1] Algorithm parameters
L = 60;                                           % number of masks
N = 2;                                            % undersampling factor              

%   [1.2] Optical wave field parameters
lambda = 532e-9;                                  % wavelength, in m
delta_detector = 1.4e-6;                          % the detector pixel size, in m
delta_computation = delta_detector / N;           % the computational pixel size, in m
ImageSize = 512;                                  % the computational pixel number 
d = 2.155e-2;                                     % propagation distance, in m

%%  [2] Simulate undersampling data (Rayleigh¨CSommerfeld model)
%   [2.1] Generate complex-domain singal
Amplitude = double(imread('./data/barbara.png')); 
Amplitude = imresize(Amplitude,[ImageSize,ImageSize]);
Phase = double(imread('./data/peppers.png'));
Phase = imresize(Phase,[ImageSize,ImageSize]);
Phase = (pi/2)*(Phase./max(Phase(:)));
x = Amplitude.*exp(1j*Phase);   

%   [2.2] Generate masks
rand('seed',101), randn('seed',11);
[yN,xN] = size(x);
Masks = zeros(yN,xN,L);        
Mask_scalled = zeros(yN,xN,L);
mask_STD = 0.3*pi;
for ll = 1:L
    for sy = 1:yN/N
        for sx = 1:xN/N
            Mask_scalled((sy-1)*N+(1:N),(sx-1)*N+(1:N),ll) = exp(1j*(randn(1,1)*(mask_STD)));
        end
    end
end
Masks = Mask_scalled;
Masks = exp(1j*(angle(Masks)));

%   [2.3] Generate measurements (Observation)
z = zeros(size(x));
for p = 1:L                                       % forward propagation    
    temp0 = AngularSpectrum(x.*(Masks(:,:,p)),d,lambda,delta_computation);  
    z(:,:,p) = abs(temp0).^2;                 
end
clear temp0
z = add_gaussion_noise(z,8);                      % add gaussion noise to measurements

if N~= 1                                          % simulate undersampling 
    for ll = 1:L                                  % average pixles intensities in NxN areas of the measurements
        for sy = 1:yN/N
            for sx = 1:xN/N
                temp1 = z(( sy-1)*N+(1:N),(sx-1)*N+(1:N),ll);
                z(( sy-1)*N+(1:N),(sx-1)*N+(1:N),ll) = ones(N,N)*mean(temp1(:));
            end
        end
    end
end
clear temp1

%%  [3] Start CD-PSR reconstruction
%   [3.1] Start Conv-PSR reconstruction
fprintf('------Start Conv-PSR reconstruction------ \n');
[rec_ConvPSR,PSNR_ConvPSR,SSIM_ConvPSR] = Conv_PSR(z,x,Masks,d,lambda,delta_computation);
%   [3.2] Start SR-SPAR reconstruction
fprintf('------Start SR-SPAR  reconstruction------ \n');
[rec_SRSPAR,PSNR_SRSPAR,SSIM_SRSPAR] = SR_SPAR(z,x,Masks,d,lambda,delta_computation);
%   [3.3] Start AS-PSR reconstruction
fprintf('------Start AS-PSR   reconstruction------ \n');
[rec_ASPSR,PSNR_ASPSR,SSIM_ASPSR] = AS_PSR(z,x,N,Masks,d,lambda,delta_computation);
%   [3.4] Start DOTV-PSR reconstruction
fprintf('------Start DOTV-PSR   reconstruction------ \n');
[rec_DOTVPSR,PSNR_DOTVPSR,SSIM_DOTVPSR] = DOTV_PSR(z,x,N,Masks,d,lambda,delta_computation);
%   [3.5] Start DONet-PSR reconstruction
fprintf('------Start DONet-PSR  reconstruction------ \n');
[rec_DONetPSR,PSNR_DONetPSR,SSIM_DONetPSR] = DONet_PSR(z,x,N,Masks,d,lambda,delta_computation);

%%  [4] Show results
figure;
subplot(2,5,1);imshow(abs(rec_ConvPSR),[]);title('Conv-PSR-Amp');
subplot(2,5,6);imshow(angle(rec_ConvPSR),[]);title('Conv-PSR-Pha');
subplot(2,5,2);imshow(abs(rec_SRSPAR),[]);title('SR-SPAR-Amp');
subplot(2,5,7);imshow(angle(rec_SRSPAR),[]);title('SR-SPAR-Pha');
subplot(2,5,3);imshow(abs(rec_ASPSR),[]);title('AS-PSR-Amp');
subplot(2,5,8);imshow(angle(rec_ASPSR),[]);title('AS-PSR-Pha');
subplot(2,5,4);imshow(abs(rec_DOTVPSR),[]);title('DOTV-PSR-Amp');
subplot(2,5,9);imshow(angle(rec_DOTVPSR),[]);title('DOTV-PSR-Pha');
subplot(2,5,5);imshow(abs(rec_DONetPSR),[]);title('DONet-PSR-Amp');
subplot(2,5,10);imshow(angle(rec_DONetPSR),[]);title('DONet-PSR-Pha');

fprintf('Comparison of reconstruction quality (Amplitude): \n');
fprintf('Algorithm: Conv-PSR ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_ConvPSR,SSIM_ConvPSR);
fprintf('Algorithm: SR-SPAR  ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_SRSPAR,SSIM_SRSPAR);
fprintf('Algorithm: AS-PSR   ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_ASPSR,SSIM_ASPSR);
fprintf('Algorithm: DOTV-PSR   ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_DOTVPSR,SSIM_DOTVPSR);
fprintf('Algorithm: DONet-PSR  ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_DONetPSR,SSIM_DONetPSR);