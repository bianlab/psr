%% Adaptive smoothing pixel super-resolution (AS-PSR)

function [rec_ASPSR,PSNR_ASPSR,SSIM_ASPSR] = AS_PSR(z,x,N,Masks,d,lambda,delta_computation)
% Parameters of algorithm
[yN,xN,L] = size(z);
[yyy,xxx,zzz]=size(Masks);

iters = 22;                % Iteration number
alpha = 0;  
beta = 0.3;  
err = 0;
threshold1 = 0.5*(-2.5e+03);
threshold2 = 0.7*(-2.5e+03);

% Initialization for the object wavefront
Us1=sqrt(z).*exp(1j*zeros(size(z)));
temp=zeros(yN,xN,L);
for ww=1:L                 % backward propagation to the object plane
    qqq=AngularSpectrum(Us1(:,:,ww), -d,lambda,delta_computation);
    temp(:,:,ww)=(((Masks(:,:,ww))).^(-1)).*qqq;
end
clear qqq
V=mean(temp,3); 

for ss = 1:iters
    %[1] forward propagation to the sensor plane
    for ww=1:L          
        Us(:,:,ww)=AngularSpectrum((Masks(:,:,ww)).*V, d,lambda,delta_computation);
        a_sp = abs(Us(:,:,ww));
        for sy=1:yN/N
            for sx=1:xN/N
                temp00 = a_sp(( sy-1)*N+(1:N),(sx-1)*N+(1:N));
                a_avg = sum(temp00(:))/(N.^2);

                a_sp(( sy-1)*N+(1:N),(sx-1)*N+(1:N))=a_avg+alpha*(a_sp(( sy-1)*N+(1:N),(sx-1)*N+(1:N))-a_avg);
                temp01 = a_sp(( sy-1)*N+(1:N),(sx-1)*N+(1:N));
                i_p1 = sum(temp01(:).^2);

                temp02 = z(( sy-1)*N+(1:N),(sx-1)*N+(1:N),ww);
                i_p = temp02(1,1);   

                a_sp(( sy-1)*N+(1:N),(sx-1)*N+(1:N)) = (sqrt(i_p/i_p1) ) .*  a_sp(( sy-1)*N+(1:N),(sx-1)*N+(1:N));
            end
        end
        Us(:,:,ww) = (  (1-beta).* a_sp + beta .* abs(Us(:,:,ww)) ) .* exp (1j.*(angle(Us(:,:,ww))));
    end     
    %[3] backward propagation to the object plane
    for ww=1:L         
        qqq=AngularSpectrum(Us(:,:,ww), -d,lambda,delta_computation);
        temp(:,:,ww)=(((Masks(:,:,ww))).^(-1)).*qqq;
    end
    clear qqq 
    V = mean(temp,3);
    clear temp
    
    % [4] Update alpha and beta
    M = xN*yN/N/N;
    
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
    delta_err  =   sqrt(AAAA/M/L) - err ;
    err =  sqrt(AAAA/M/L);
    
    % Update alpha
    if  delta_err< threshold1
        alpha = 1-(1-alpha)/2 ; 
        threshold1 = threshold1/2;
    end
    % Update beta
    if  delta_err< threshold2
        beta = 1-(1-beta)/2 ; 
        threshold2 = threshold2/2;
    end
   
end
rec_ASPSR = exp(-1i*(angle(trace(x((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:)'*V((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:))))) * V((yyy-yyy)/2+1:(yyy+yyy)/2,(xxx-xxx)/2+1:(xxx+xxx)/2,:);  %% Phase correction by an invariant shift according to x
PSNR_ASPSR = psnr(abs(rec_ASPSR)./max(max(abs(rec_ASPSR))),abs(x)./max(max(abs(x))));
SSIM_ASPSR = ssim(abs(rec_ASPSR)./max(max(abs(rec_ASPSR))),abs(x)./max(max(abs(x))));
end
