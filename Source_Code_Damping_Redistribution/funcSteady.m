function [state0, steady_exitflag] = funcSteady()
% return 0 if successful, otherwise return 1
    options = optimoptions('fsolve','MaxIterations', 1200, 'MaxFunctionEvaluations', 4000, 'Display', 'off');

    [state0, steady_fval, steady_exitflag] = fsolve(@calcPowerflow, calcPowerflowInit(), options);

    if steady_exitflag <= 0
        warning('Steady value warning: funcSteady warning');
        steady_exitflag = 1;
    else
        steady_exitflag = 0;
    end
end