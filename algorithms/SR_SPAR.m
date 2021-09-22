%% DEMO_SUPER_RESOLUTION_SPARSE_PHASE_RETRIEVAL.m

%% SPAR (WRAPPED and UNWRAPPED) PHASE and AMPLITUDE RECONSTRUCTION for Phase Retrieval
%% Subpixel imaging
%% Random Phase Coded Aperture
%% Demo-Program for the paper "Computational super-resolution phase retrieval from multiple phase-coded diffractional patterns: experimental study", V.Katkovnik, I.Shevkunov, N.V.PEtrov, K.Egiazarian.
%% I.Shevkunov,  2017
%% 
function [rec_SRSPAR,PSNR_SRSPAR,SSIM_SRSPAR] = SR_SPAR(z,x,Masks,d,lambda,delta_computation)
%% PARAMETERS of BM3D-FRAME FILTER
Nstep = 1;  N11=2;  N22=2;
threshType='h';
threshold_ampl=2;
threshold_phi=2;
[yN,xN,L]=size(z);
[yyy,xxx,zzz]=size(Masks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS of ALGORITHMS
filtering=1;               %%  1 for BM3D filtering
unwrapping=0;              %%  1 -means unwrapping on before BM3D phase filtering, 0- means unwrapping off

                           
IterNum=15;                %% ITERATION NUMBER

zp=1;                     %% the factor of zero padding, (resulting one dimmension size)=(initial one dimmension size)*zp

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% SPAR algorithm %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic                 %% for SPAR
%% INITIALIZATION FOR SPAR
%%%%%%%%%%%% NEW INITIALIZATION from z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Us=sqrt(z).*exp(1j*zeros(size(z)));
temp=zeros(yN,xN,L);
for ww=1:L          %% backward propagation to the object plane
    qqq=AngularSpectrum(Us(:,:,ww), -d,lambda,delta_computation);
    temp(:,:,ww)=(((Masks(:,:,ww))).^(-1)).*qqq;
end
clear qqq
V=mean(temp,3); 
%%%%%%%%%%%%%%%%%%%%% INITIALIZATION from Z, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_init=V;           %%first estimate for the object wavefront

%% Phase unwrapping
if unwrapping==1
    potential.quantized = 'no';
    potential.threshold = pi;
    verbose='no';
    phi_Uo_init=angle(V_init);
    phi_Uo_init1=puma_ho(angle(V_init((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,1)),.5,'verbose',verbose); %% for demo unwrapping from raw data
    bb_init=mean(varphi(:) - phi_Uo_init1(:));
    phi_Uo_init((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)=phi_Uo_init1;
    RMSE_phi(1) = norm(varphi - phi_Uo_init((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,1)-bb_init, 'fro')/sqrt(yN/zp*xN/zp); %% RMSE Initial rel. phase error
    
else
    V_cor_init=exp(-1i*angle(trace(x'*V_init))) * V_init;   %% PHASE CORRECTED INITIAL GUESS
    bb1_init=-angle(trace(x'*V_init));
    phi_Uo_init=angle(V_init);
    Uo_corr_init=exp(-1i*angle(trace(x((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)'*V_init((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)))) * V_init((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2);   %% CORRECTED INITIAL GUESS
    RMSE_phi(1) = norm(wrap(angle(x((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)) - angle(Uo_corr_init)), 'fro')/sqrt(yN/zp*xN/zp); % EMSE Initial rel. phase error
end

%% RMSE of the phase and amplitude estimates for use in BM3D filtered estimates
sigma_phi(1)=function_stdEst2D(phi_Uo_init((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2),2);
sigma_ampl(1)=function_stdEst2D(abs(V((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)),2);

%% MAIN ITERATIONS of SPAR PHASE RETRIEVAL
KAPPA=1e-3;                                  %% noise parameter                  
gamma_1=1/KAPPA;
%%%%%%%%%%%%% TRUE INITIALIZATION
for ss=1:IterNum
    ss0 = ss;
    %% STEP 1:     %% Forward propagation
    for ww=1:L          %% forward propagation to the sensor plane
        Us(:,:,ww)=AngularSpectrum((Masks(:,:,ww)).*V, d,lambda,delta_computation);
    end
    
    
    %% STEP 2:     %% Filtering
       [Us]=NoiseFiltering_3D_propag(Us,z,KAPPA,gamma_1,L); %% Filtering at Sensor Plane
       
    %% STEP 3: Backward propagation
    for ww=1:L          %% backward propagation to the object plane
        qqq=AngularSpectrum(Us(:,:,ww), -d,lambda,delta_computation);
        temp(:,:,ww)=(((Masks(:,:,ww))).^(-1)).*qqq;
    end
    clear qqq 
    Uo=mean(temp,3);
    clear temp
    
    %% STEP 4: Phase unwrapping ;
    if unwrapping==1
        potential.quantized = 'no';
        potential.threshold = pi;
        verbose='no';
        phi_Uo=double(angle(Uo));
        phi_Uo((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,1)=double(puma_ho(angle(Uo((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,1)),.5,'verbose',verbose)); %% for demo unwrapping from raw data
        
    else
        
        phi_Uo=double(angle(Uo));
    end
    
    %% STEP 5: BM3D phase and amplitude filtering
    sigma_phi(ss+1)=function_stdEst2D(phi_Uo((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2),2); %% STD of phase
    sigma_ampl(ss+1)=function_stdEst2D(abs(Uo((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)),2);
    
    if filtering==1
            phi_u0_SPAR=phi_Uo;
            B_u0_SPAR=abs(Uo);
            phi_u0_SPAR((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)=BM3D_SPAR_UNWRAP_PHIp(phi_Uo((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2),threshType, sigma_phi(ss+1)*threshold_phi, N11, N22, Nstep,filtering,ss,ss0);
            B_u0_SPAR((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)=BM3D_SPAR_ABSp(abs(Uo((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)),threshType, sigma_ampl(ss+1)*threshold_ampl, N11, N22, Nstep,filtering,ss,ss0);
    else
        phi_u0_SPAR=phi_Uo;
        B_u0_SPAR=abs(Uo);
    end
    
    %% STEP 6: UPDATE of object estimate
    
    V=B_u0_SPAR.*exp(1j*phi_u0_SPAR);
    
    if 1  %%take only central part of the object plane; values outside object that was equal zero, zeros again
        Vtemp=V((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,:);
        V=zeros(size(V));
        V((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,:)=Vtemp;
        clear Vtemp
    end
    
    rec_SRSPAR = exp(-1i*(angle(trace(x((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,:)'*V((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,:))))) * V((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,:);  %% Phase correction by an invariant shift according to x
    
%     if unwrapping==1
%         phi_u0_SPARcenter=phi_u0_SPAR((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,:);
%         bb=mean(varphi(:) - phi_u0_SPARcenter(:));                                       %% absolute phase shift
%         RMSE_phi(ss+1) = norm(varphi - phi_u0_SPARcenter-bb, 'fro')/sqrt(yN/zp*xN/zp); %% RMSE absolute phase error
%         RMSE_phi_center(ss+1) = norm(wrap(angle(x((yyy*zp-yyy/4)/2+1:(yyy*zp+yyy/4)/2,(xxx*zp-xxx/4)/2+1:(xxx*zp+xxx/4)/2,:)) - (angle(V_cor(yyy/2-yyy/8+1:yyy/2+yyy/8,xxx/2-xxx/8+1:xxx/2+xxx/8)))), 'fro')/sqrt(yN/zp/4*xN/zp/4);
%     else
%         
%         RMSE_phi(ss+1) = norm(wrap(angle(x((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,:)) - (angle(V_cor))), 'fro')/sqrt(yN/zp*xN/zp); %% RMSE phase error
    end
%        RMSE_ampl(ss+1)=norm(abs(x((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2))-abs(V((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2,:)),'fro')/sqrt(yN/zp*xN/zp);                     %% RMSE amplitude error
    
    
%     if 1  %% ITERATION SHOW
%         if unwrapping==1
%            
%                 figure(1),
%                 subplot(2,2,1), imshow(phi_u0_SPAR((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)+bb,[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(ss+1),3)]),...
%                 subplot(2,2,2), imshow(varphi,[]), title('Original phase'),...
%                 subplot(2,2,3), imshow(abs(V((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(ss+1),3)])
%                 subplot(2,2,4), imshow(phi_u0_SPAR((yyy*zp-yyy)/2+1:(yyy*zp+yyy)/2,(xxx*zp-xxx)/2+1:(xxx*zp+xxx)/2)+bb-varphi,[]), title('Phase difference'), colorbar 
%                 pause(.1)
%         else
%                 figure(1), 
%                 subplot(2,2,1), imshow(angle(V_cor),[]), title(['Rec Phase, IterNum=', num2str(ss+1),'; RMSE_{phase}=', num2str(RMSE_phi(ss+1),3)]),...
%                 subplot(2,2,2), imshow(varphi,[]), title('Original Phase'),...    
%                 subplot(2,2,3), imshow(abs(V_cor),[]), title(['Ampl, IterNum=', num2str(ss+1),'; RMSE_{ampl}=', num2str(RMSE_ampl(ss+1),3)])
%                 subplot(2,2,4), imshow(angle(V_cor)-varphi,[]), title('Phase difference'), colorbar 
%                  pause(0.1)
%         end
%     end
PSNR_SRSPAR = psnr(abs(rec_SRSPAR)./max(max(abs(rec_SRSPAR))),abs(x)./max(max(abs(x))));
SSIM_SRSPAR = ssim(abs(rec_SRSPAR)./max(max(abs(rec_SRSPAR))),abs(x)./max(max(abs(x))));
end