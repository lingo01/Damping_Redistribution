% compute the active and reactive power at each bus
% Parameters:
%   x = [delta_1 V_1 delta_2 V_2 ...]
function power = calcPower(x)
    global netG netB nNet;
    
    P = zeros(nNet,1);
    Q = zeros(nNet,1);
    
    % 功率平衡
    for ii = 1:nNet
        tmpP = 0;    tmpQ = 0;
        for jj = 1:nNet
            tmpP = tmpP + x(2*jj) * ( netG(ii,jj)*cos(x(2*ii-1)-x(2*jj-1)) + netB(ii,jj)*sin(x(2*ii-1)-x(2*jj-1)) );
            tmpQ = tmpQ + x(2*jj) * ( netB(ii,jj)*cos(x(2*ii-1)-x(2*jj-1)) - netG(ii,jj)*sin(x(2*ii-1)-x(2*jj-1)) );
        end
        P(ii) =  tmpP * x(2*ii);
        Q(ii) = -tmpQ * x(2*ii);

    end

    power = P + 1j*Q;
end

