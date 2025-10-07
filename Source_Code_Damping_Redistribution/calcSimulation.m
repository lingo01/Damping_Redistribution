% compute the derivatives of the power system dynamic state
% Parameters:
%   x = [delta_1 omega_1 V_1 delta_2 omega_2 V_2 ...]

function dx = calcSimulation(t,x)
    % dx = [delta, omega, V]
    
    global type Gen CDInv QDInv Load Cnct;
    global nNet;
    global LoadEquivTime LoadEquivDroop;

    % x(1:3:(3*nNet-2)) = x(1:3:(3*nNet-2)) - x(3*nNet-2);

    dx = zeros(3*nNet,1);
    
    
    % solve the power at each bus
    xReducted = funStateTrans(x);
    S = calcPower(xReducted);   
    Pout = real(S); Qout = imag(S);

    % compute the derivatives
    for ii = 1:nNet
        if strcmp(type{ii}, 'Gen')
            dx(3*ii-2) = x(3*ii-1);
            dx(3*ii-1) = 1/Gen{ii}.M  * (-Gen{ii}.D*x(3*ii-1) + Gen{ii}.Pm - Pout(ii));
            dx(3*ii)   = 1/Gen{ii}.Td * (-x(3*ii) - (Gen{ii}.xd-Gen{ii}.xdp)*Qout(ii)/x(3*ii) + Gen{ii}.Ef);
        elseif strcmp(type{ii}, 'CDInv')
            dx(3*ii-2) = 1/CDInv{ii}.t1 * ( -(x(3*ii-2)-CDInv{ii}.As) - CDInv{ii}.D1 * (Pout(ii) - CDInv{ii}.Ps) );
            dx(3*ii-1) = 0;
            dx(3*ii)   = 1/CDInv{ii}.t2 * ( -(x(3*ii)  -CDInv{ii}.Vs) - CDInv{ii}.D2 * (Qout(ii) - CDInv{ii}.Qs) );
        elseif strcmp(type{ii}, 'QDInv')
            dx(3*ii-2) = 1/QDInv{ii}.t1 * ( -(x(3*ii-2)-QDInv{ii}.As) - QDInv{ii}.D1 * (Pout(ii) - QDInv{ii}.Ps) );
            dx(3*ii-1) = 0;
            us = QDInv{ii}.Vs + QDInv{ii}.D2 * QDInv{ii}.Qs/QDInv{ii}.Vs;
            dx(3*ii)   = 1/QDInv{ii}.t2 * ( -QDInv{ii}.D2 * Qout(ii) - x(3*ii) * (x(3*ii)-us) );
        elseif strcmp(type{ii}, 'Load')
            dx(3*ii-2) = 1/LoadEquivTime * ( -(x(3*ii-2)-0) - LoadEquivDroop * (Pout(ii) - Load{ii}.Ps) );
            dx(3*ii-1) = 0;
            dx(3*ii)   = 1/LoadEquivTime * ( -(x(3*ii)  -1) - LoadEquivDroop * (Qout(ii) - Load{ii}.Qs) );
        elseif strcmp(type{ii}, 'Cnct')
            dx(3*ii-2) = 1/LoadEquivTime * ( -(x(3*ii-2)-0) - LoadEquivDroop * (Pout(ii) - 0) );
            dx(3*ii-1) = 0;
            dx(3*ii)   = 1/LoadEquivTime * ( -(x(3*ii)  -1) - LoadEquivDroop * (Qout(ii) - 0) );
        end
    end
    
end

