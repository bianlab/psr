function [prop_wave, TF] = AngularSpectrum(wavefront,distance,lambda,dx)
% propagation using angular spectrum...  BY IGOR
% dx - the sample spacings in the x and y in meters 
% directions of the input object 'obj'.
% z is the propagation distance in m
% lambda is the wavelength in m


%% Setup matrices representing reciprocal space coordinates
[K,L] = size(wavefront);
k = -K/2:K/2 - 1;
l = -L/2:L/2 - 1;
[k,l] = meshgrid(k,l);
U=1 -lambda^2*((k/(dx*K)).^2+(l/(dx*L)).^2);
TF= exp(1i*2*pi/lambda*distance*sqrt(U));
TF(U<0)=0;%%%%TRANSFER FUNCTION

prop_wave =ifftshift(ifft2(ifftshift(TF.*fftshift(fft2(fftshift(wavefront))))));



