%% Conventional pixel super-resolution phase reteieval (Conv-PSR)
function [rec_ConvPSR,PSNR_ConvPSR,SSIM_ConvPSR] = Conv_PSR(z,x,Masks,d,lambda,delta_computation)
IterNum = 10;
[xN,yN,L] = size(z);
[yyy,xxx,zzz]=size(Masks);
% Initialization for the object wavefront
Us = sqrt(z).*exp(1j*zeros(size(z))); 
temp  =zeros(yN,xN,L);
for ww = 1:L          
    qqq = AngularSpectrum(Us(:,:,ww), -d,lambda,delta_computation);
    temp(:,:,ww) = (((Masks(:,:,ww))).^(-1)).*qqq;
end
clear qqq
V = mean(temp,3);    
      
for ss=1:IterNum   
    % STEP 1: Forward propagation    
    for ww=1:L          
        Us(:,:,ww)=AngularSpectrum((Masks(:,:,ww)).*V, d,lambda,delta_computation);
    end
        
    % STEP 2: Amplitude constraint 
    for s=1:L  
        Us(:,:,s)=(sqrt(z(:,:,s))).*exp(1j*angle(Us(:,:,s)));   
    end
    % STEP 3: Backward propagation
    for ww=1:L          
        qqq=AngularSpectrum(Us(:,:,ww), -d,lambda,delta_computation);
        temp(:,:,ww)=(((Masks(:,:,ww))).^(-1)).*qqq;
    end
    clear qqq 
    V = mean(temp,3);
    clear temp    
end
rec_ConvPSR = exp(-1i*(angle(trace(x((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:)'*V((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:))))) * V((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:);  %% Phase correction by an invariant shift according to x
PSNR_ConvPSR = psnr(abs(rec_ConvPSR)./max(max(abs(rec_ConvPSR))),abs(x)./max(max(abs(x))));
SSIM_ConvPSR = ssim(abs(rec_ConvPSR)./max(max(abs(rec_ConvPSR))),abs(x)./max(max(abs(x))));