function denoisedv = ffdnet_denoise(noisyv,para)
format compact;
global sigmas; % input noise level or input noise level map
useGPU = false;
ffdnetvnorm = true;
if isfield(para,'useGPU'),      useGPU = para.useGPU;             end
if isfield(para,'ffdnetvnorm'), ffdnetvnorm = para.ffdnetvnorm;   end

if ffdnetvnorm
    %%%%% [start] normalization, like VBM4D %%%%%
    maxz = max(noisyv(:));
    minz = min(noisyv(:));
    scale = 0.7;
    shift = (1-scale)/2;
    noisyv = (noisyv-minz)/(maxz-minz);
    noisyv = noisyv*scale+shift;

    sigmas = para.sigma/(maxz-minz)*scale;
    %%%%% [start] normalization, like VBM4D %%%%%
else
    % set noise level map
    sigmas = para.sigma;
end

[~,~,nframe] = size(noisyv); 

if isfield(para,'net') && ~isempty(para.net)
    net = para.net;
else
    load(fullfile('models','FFDNet_gray.mat'),'net');
    net = vl_simplenn_tidy(net);
    if useGPU
        net = vl_simplenn_move(net, 'gpu') ;
    end
end

for iframe = 1:nframe
    % input = double(noisyv(:,:,iframe));
    input = single(noisyv(:,:,iframe));
    
    % perform denoising
    % res    = vl_simplenn_de(net,input,[],[],'conserveMemory',true,'mode','test'); % matconvnet default
    if useGPU
        input = gpuArray(input);
        res    = vl_simplenn_de(net,input,[],[],'conserveMemory',true,'mode','test'); % matconvnet default
        output = res(end).x;
        output = gather(output);
    else
        res    = vl_simplenn_de(net,single(input),[],[],'conserveMemory',true,'mode','test'); % matconvnet default
        output = res(end).x;
    end
    % res    = vl_ffdnet_concise(net, input);    % concise version of vl_simplenn_de for testing FFDNet
    % res    = vl_ffdnet_matlab(net, input); % use this if you did  not install matconvnet; very slow
    
    denoisedv(:,:,iframe) = double(output);
end

if ffdnetvnorm
    %%%%% [start] de-normalization, like VBM4D %%%%%
    denoisedv = (denoisedv-shift)/scale;
    denoisedv = denoisedv*(maxz-minz)+minz;
    %%%%% [start] de-normalization, like VBM4D %%%%%
end

end

