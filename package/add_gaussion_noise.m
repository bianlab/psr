function [Y] = add_gaussion_noise(X,SNR)

X = double(X);
signal_power = sum(X(:).^2);
noise_power = signal_power/(10^(SNR/10));
picture_size = size(X);
noise_var = noise_power/prod(picture_size(:));
Y = X+sqrt(noise_var)*randn(picture_size);
end

