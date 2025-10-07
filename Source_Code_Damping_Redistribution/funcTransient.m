function [t, state, exit_flag] = funcTransient(state0, steady_exitflag, t0)
% return 0 if successful, otherwise return 1
    global dt te;
    global refPos;
    global netB netG;
    global shortPos shortRes shortCT shortType;

    if ~exist('t0', 'var')
        t0 = 0.1 * (steady_exitflag == 0) + 10 * (steady_exitflag == 1);
    end

    % before fault transient
    [t_1,state_1] = ode15s(@calcSimulation, [0:dt:t0], funStateTrans(state0, 'Expand') );  

    if strcmp(shortType, 'R')
        netG(shortPos,shortPos) = netG(shortPos,shortPos) - shortRes;   
    elseif strcmp(shortType, 'X')
        netB(shortPos,shortPos) = netB(shortPos,shortPos) - shortRes;   
    end

    % during fault transient
    [t_2,state_2] = ode15s(@calcSimulation, [t0:dt:(t0+shortCT)], state_1(end,:));

    if strcmp(shortType, 'R')
        netG(shortPos,shortPos) = netG(shortPos,shortPos) + shortRes;   
    elseif strcmp(shortType, 'X')
        netB(shortPos,shortPos) = netB(shortPos,shortPos) + shortRes;   
    end

    % post fault transient
    [t, state] = ode15s(@calcSimulation, [(t0+shortCT):dt:te],  state_2(end,:));

    % output results
    t = [t_1; t_2; t];
    state = [state_1; state_2; state];
    state(:,1:3:end) = state(:,1:3:end) - state(:,3*refPos-2);

    exit_flag = (t(end) ~= te);
end