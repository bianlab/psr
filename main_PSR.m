% Comparison of pixel super-resolution (PSR) algorithms in digital holography. 
% Algotithm lists: 1. Conventional PSR algorithm (Conv-PSR) [1]
%                  2. Super-resolution sparse phase-amplitude retrieval(SR-SPAR) [2]
%                  3. Adaptive smoothing PSR algorithm (AS-PSR) [3]
%                  4. Plug-and-play Total-Variation PSR algorithm (PNPTV-PSR) [Ours]
%                  5. Plug-and-play (Denoising network, FFDNet) PSR algorithm (PNPNet-PSR) [Ours]

% Reference
%   [1] Fienup J R. Phase retrieval algorithms: a comparison[J]. 
%       Applied optics, 1982, 21(15): 2758-2769.
%   [2] Katkovnik V, Shevkunov I, Petrov N V, et al. Computational 
%       super-resolution phase retrieval from multiple phase-coded 
%       diffraction patterns: simulation study and experiments[J]. 
%       Optica, 2017, 4(7): 786-794.
%   [3] Gao Y, Cao L. High-fidelity pixel-super-resolved complex field 
%       reconstruction via adaptive smoothing[J]. Optics Letters, 
%       2020, 45(24): 6807-6810.


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

%%  [2] Simulate undersampling data (Rayleigh?CSommerfeld model)
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

%%  [3] Start PSR reconstruction
%   [3.1] Start Conv-PSR reconstruction
fprintf('------Start Conv-PSR reconstruction------ \n');
[rec_ConvPSR,PSNR_ConvPSR,SSIM_ConvPSR] = Conv_PSR(z,x,Masks,d,lambda,delta_computation);
%   [3.2] Start SR-SPAR reconstruction
fprintf('------Start SR-SPAR  reconstruction------ \n');
[rec_SRSPAR,PSNR_SRSPAR,SSIM_SRSPAR] = SR_SPAR(z,x,Masks,d,lambda,delta_computation);
%   [3.3] Start AS-PSR reconstruction
fprintf('------Start AS-PSR   reconstruction------ \n');
[rec_ASPSR,PSNR_ASPSR,SSIM_ASPSR] = AS_PSR(z,x,N,Masks,d,lambda,delta_computation);
%   [3.4] Start PNPTV-PSR reconstruction
fprintf('------Start PNPTV-PSR   reconstruction------ \n');
[rec_PNPTVPSR,PSNR_PNPTVPSR,SSIM_PNPTVPSR] = PNPTV_PSR(z,x,N,Masks,d,lambda,delta_computation);
%   [3.5] Start PNPNet-PSR reconstruction
fprintf('------Start PNPNet-PSR  reconstruction------ \n');
[rec_PNPNetPSR,PSNR_PNPNetPSR,SSIM_PNPNetPSR] = PNPNet_PSR(z,x,N,Masks,d,lambda,delta_computation);

%%  [4] Show results
figure;
subplot(2,5,1);imshow(abs(rec_ConvPSR),[]);title('Conv-PSR-Amp');
subplot(2,5,6);imshow(angle(rec_ConvPSR),[]);title('Conv-PSR-Pha');
subplot(2,5,2);imshow(abs(rec_SRSPAR),[]);title('SR-SPAR-Amp');
subplot(2,5,7);imshow(angle(rec_SRSPAR),[]);title('SR-SPAR-Pha');
subplot(2,5,3);imshow(abs(rec_ASPSR),[]);title('AS-PSR-Amp');
subplot(2,5,8);imshow(angle(rec_ASPSR),[]);title('AS-PSR-Pha');
subplot(2,5,4);imshow(abs(rec_PNPTVPSR),[]);title('PNPTV-PSR-Amp');
subplot(2,5,9);imshow(angle(rec_PNPTVPSR),[]);title('PNPTV-PSR-Pha');
subplot(2,5,5);imshow(abs(rec_PNPNetPSR),[]);title('PNPNet-PSR-Amp');
subplot(2,5,10);imshow(angle(rec_PNPNetPSR),[]);title('PNPNet-PSR-Pha');

fprintf('Comparison of reconstruction quality (Amplitude): \n');
fprintf('Algorithm: Conv-PSR ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_ConvPSR,SSIM_ConvPSR);
fprintf('Algorithm: SR-SPAR  ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_SRSPAR,SSIM_SRSPAR);
fprintf('Algorithm: AS-PSR   ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_ASPSR,SSIM_ASPSR);
fprintf('Algorithm: PNPTV-PSR   ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_PNPTVPSR,SSIM_PNPTVPSR);
fprintf('Algorithm: PNPNet-PSR  ---- PSNR: %2.2f dB ---- SSIM: %2.2f \n',PSNR_PNPNetPSR,SSIM_PNPNetPSR);