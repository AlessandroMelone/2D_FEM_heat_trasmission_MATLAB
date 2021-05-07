function x_timeseries = solve_explicitEuler_analyticSolution(fnc,q, x0,T,q0, Dt, time_initial, time_final)

    %solve the system: \dot x = fnc(x,u)

    time = (time_initial:Dt:time_final)';
    n_steps = size(time,1);
    n_dof = size(x0,1);
    x_evolution = zeros(n_dof, n_steps);
    x_evolution(:,1) = x0;
    
    u = [0; T; q0];                    %u = [mu; T; q0];
    
    controller_status = 0;
    for k = 1:(n_steps-1)
        mu = exp(q(time(k)));
        u(1) = mu;
        x_evolution(:,k+1) = x_evolution(:,k) + Dt*fnc(x_evolution(:,k), u);
    end
    
    x_timeseries = timeseries(x_evolution',time);
    
end

