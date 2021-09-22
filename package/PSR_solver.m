function   [V]= PSR_solver (V,z,Masks,d,lambda,delta_x,xN,yN,N,alpha,beta,L)
parfor ww=1:L          
	Us(:,:,ww)=AngularSpectrum((Masks(:,:,ww)).*V, d,lambda,delta_x);
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
    % backward propagation to the object plane
for ww=1:L          
	qqq=AngularSpectrum(Us(:,:,ww), -d,lambda,delta_x);
	temp(:,:,ww)=(((Masks(:,:,ww))).^(-1)).*qqq;
end
clear qqq 
V = mean(temp,3);
clear temp
    
  