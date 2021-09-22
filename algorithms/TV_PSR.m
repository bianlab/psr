%%  Total Variation optimization pixel super-resolution phase reteieval algorithm (TV-PSR)
function [rec_TVPSR,PSNR_TVPSR,SSIM_TVPSR] = TV_PSR(z,x,N,Masks,d,lambda,delta_computation)
[yN,xN,L]=size(z);
[yyy,xxx,zzz]=size(Masks);

IterNum_initial = 15;                             % initial iteration number
IterNum_GAP = 5;                                 % GAP iteration number (PNP optimization)
lambda_gap = 0.1;                                 % weight coefficient 

 
% Initialization for the object wavefront
Us = sqrt(z).*exp(1j*zeros(size(z)));
temp = zeros(yN,xN,L);
for ww=1:L          % backward propagation to the object plane
    qqq=AngularSpectrum(Us(:,:,ww), -d,lambda,delta_computation);
    temp(:,:,ww)=(((Masks(:,:,ww))).^(-1)).*qqq;
end
clear qqq
V=mean(temp,3); 

alpha = 1;
beta = 0.2;
for ss = 1:IterNum_initial

    % forward propagation to the sensor plane
    [V] =  PSR_solver (V,z,Masks,d,lambda,delta_computation,xN,yN,N,alpha,beta,L);    
    AAAA = 0;
    for ii1 = 1:L        
          for sy=1:yN/N
              for sx=1:xN/N
                  temp10 = Us(( sy-1)*N+(1:N),(sx-1)*N+(1:N),ii1);
                  U_sp = (abs(temp10)).^2;
                  
                  temp11 =  z(( sy-1)*N+(1:N),(sx-1)*N+(1:N),ii1);
                  i_p = temp11(1,1);
                  
                  AAAA = AAAA + ((  i_p  - sum(U_sp(:))   ).^2);
              end
          end
    end
    V_cor=exp(-1i*(angle(trace(x((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:)'*V((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:))))) * V((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:);  %% Phase correction by an invariant shift according to x
end

% GAP optimization (Total Variation denoising)
for isig = 1:IterNum_GAP 
	zb = A_nonliner (V_cor,Masks,d,lambda,delta_computation,N,L);
	z_new = (z-lambda_gap*zb);
       
	[V_cor] =  PSR_solver (V_cor,z_new,Masks,d,lambda,delta_computation,xN,yN,N,alpha,beta,L);
     
	V_cor_new = exp(-1i*(angle(trace(x((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:)'*V_cor((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:))))) * V_cor((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:);  %% Phase correction by an invariant shift according to x 
	V_cor =  V_cor_new+lambda_gap*V_cor; 
        
    phi_Uo=angle(V_cor);
    amp_Uo=abs(V_cor);

    phi_Uo = TV_denoising(phi_Uo,1,10);
    amp_Uo = TV_denoising(amp_Uo,4,10);
     
    V_cor = amp_Uo.*exp(1j*phi_Uo);
end
rec_TVPSR = V_cor;
PSNR_TVPSR = psnr(abs(rec_TVPSR)./max(max(abs(rec_TVPSR))),abs(x)./max(max(abs(x))));
SSIM_TVPSR = ssim(abs(rec_TVPSR)./max(max(abs(rec_TVPSR))),abs(x)./max(max(abs(x))));
end





