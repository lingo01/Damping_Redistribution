% setting the network parameters
function settingNetwork(filename)
    % bus number settings
    global nNet;
    nNet = 39;
    
    global netY netG netB;
    load(['.\data&figure\', filename, '.mat']);

    if size(netY, 1) ~= size(netY, 2) || size(netY, 1) ~= nNet
        error('Initial error: settingNetwork error');
    end

    netG = real(netY);
    netB = imag(netY); 
    
    global bshunt;
    bshunt = zeros(nNet,1);
    for ii = 1:nNet
        bshunt(ii) = sum(netB(ii,:));
    end
end

