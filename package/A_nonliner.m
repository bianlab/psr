function z = A_nonliner (x,Masks,d,lambda,delta_x,N,L)

z=zeros(size(x));
for ww=1:L                             % Object wavefront propagation to sensor plane    
    temp=AngularSpectrum(x.*(Masks(:,:,ww)), d,lambda,delta_x);  %forward propagation
    z(:,:,ww)=abs(temp).^2;            % observations array
end

% average pixles intensities in NxN areas of the observations
[yN,xN]=size(x);
if N~=1
    for ll = 1:L                                 
        for sy=1:yN/N
            for sx=1:xN/N
                temp0=z(( sy-1)*N+(1:N),(sx-1)*N+(1:N),ll);
                z(( sy-1)*N+(1:N),(sx-1)*N+(1:N),ll)=ones(N,N)*mean(temp0(:));
            end
        end
    end
end

