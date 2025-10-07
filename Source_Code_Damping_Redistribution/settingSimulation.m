% transient simulation settings
function settingSimulation()
    global dt te;
    global LoadEquivTime LoadEquivDroop;
    global refPos;

    dt = 0.01;  % simulation time step
    te = 20;    % simulation duration

    LoadEquivTime  = 0.005;   % load dynamics approximation equivalent time constant, default=0.05
    LoadEquivDroop = 100; % load dynamics approximation equivalent droop coefficient, default=100

    refPos = 39; % phase angle reference node
end