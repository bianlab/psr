%%
function [Us]=NoiseFiltering_3D_propag(V,I,KAPPA,gamma_1,L)
%% [Us]=NoiseFiltering_3D_propag(Us,z,KAPPA,gamma_1,L).*Masks_Support+Us.*(1-Masks_Support);
Us=zeros(size(I));
threshold_ampl=10;N11=8; N22=8; Nstep=3; threshType='h';
for s=1:L
  
    temp1=V(:,:,s); %% for Us
    
    temp2=I(:,:,s);  %% for z
  %%  temp1=(temp1); temp2=(temp2); %% for proper zeroing the central part of the sensor
    W=abs(temp1); gamma_11=gamma_1*KAPPA;
    %%%%%%%%% BM3D filtering for W 
%     if 0
%      sigma_ampl=function_stdEst2D(W,2)
%      W_filtered=BM3D_SPAR_ABSp(W,threshType, sigma_ampl*threshold_ampl, N11, N22, Nstep,1,1,1);
%     else
%         W_filtered=W;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     W=(W/gamma_11+sqrt(W.^2/gamma_11^2+4*temp2*(1+1/gamma_11)/KAPPA))/2/(1+1/gamma_11);

 %%   b_{s}=((|v_{s}|+?(|v_{s}|?+4z_{s}?(1+???)))/(2(1+???)))
    
   %% W=(W/gamma_11+sqrt(W.^2/gamma_11^2+4*temp2*(1+1/gamma_11)/KAPPA))/2/(1+1/gamma_11);
    
   %% WW=abs(temp1).*(1-Masks_Support(:,:,s))+W.*Masks_Support(:,:,s);
    %WW=ifftshift(WW); temp1=ifftshift(temp1);    %% return fft2 location of frequencies
    
    Us(:,:,s)=W.*exp(1j*angle(temp1));
   
end

end
    
    