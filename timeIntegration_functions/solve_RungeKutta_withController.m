function x_timeseries = solve_RungeKutta_withController(fnc, x0,T,freq_T,phase_T, q0,G_z, reference, idx_sensor, Dt, time_initial, time_final)
% Villani p.173: with x_evolution(k) = y_n and time(k) = x_n 

    %solve the system: \dot x = fnc(x,u)
    time = (time_initial:Dt:time_final)';
    n_steps = size(time,1);
    n_dof = size(x0,1);
    x_evolution = zeros(n_dof, n_steps);
    x_evolution(:,1) = x0;
    
    u = [0; T; q0];                    %u = [mu; T; q0];
    b = [1,2,2,1]/6;
    
    controller_status = 0;  z = 0;
    for k = 1:(n_steps-1)
        %         x_evolution(:,k+1) = x_evolution(:,k) + Dt*fnc(x_evolution(:,k), u);
        %evaluation k_i
        controller_out = controller([reference, x_evolution(idx_sensor,k), controller_status]);
        [mu, z] = filter(G_z.Numerator{1}, G_z.Denominator{1},controller_out(1),z);    
        controller_status = controller_out(2);
        u(1) = mu;
        u(2) = T + sin(freq_T * time(k) + phase_T);
        
        k1 = fnc(x_evolution(:,k), u);
        
        k2 = fnc(x_evolution(:,k)+Dt*k1/2, u);
        
        k3 = fnc(x_evolution(:,k)+Dt*k2/2, u);
        
        k4 = fnc(x_evolution(:,k)+Dt*k3, u);
       
        x_evolution(:,k+1) = x_evolution(:,k) + Dt*(b(1)*k1+b(2)*k2+b(3)*k3+b(4)*k4);
        
        if mod(k,500)==0
            disp([' integration step: ',num2str(k), ' of ',num2str(n_steps)]);
        end
    end
    
    x_timeseries = timeseries(x_evolution',time);
    
end

