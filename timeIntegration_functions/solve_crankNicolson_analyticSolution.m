function x_timeseries = solve_crankNicolson_analyticSolution(fnc,q, x0,T,q0, Dt, time_initial, time_final)
    %solve the system: \dot x = fnc(x)
    time = (time_initial:Dt:time_final)';
    n_steps = size(time,1);
    n_dof = size(x0,1);
    x_evolution = zeros(n_dof, n_steps);
    x_evolution(:,1) = x0;
    
    % Newton_Raphson parameters
    num_iter_max = 10; %max number of iteration
    toll = 0.001; %error tollerance
    
    u = [0; T; q0];
    for k = 1:(n_steps-1)
        mu = exp(q(time(k)));
        u(1) = mu;
       
        %define a function f(x)=0 s.t. for f(x*)=0 gives x*=x(k+1)
        %x_k = x_evolution(:,k); x_k+1 = x; 
        fnc_x = @(x) x-Dt/2*fnc(x,u)-x_evolution(:,k)-Dt/2*fnc(x_evolution(:,k),u);
%         %we need to use this when consider a tempo-variant disturbance T
%         fnc_x = @(x) x-Dt/2*fnc(x,time(k+1))-x_evolution(:,k)-Dt/2*fnc(x_evolution(:,k),time(k));
        
        
        %compute x(k+1)
        [x_evolution(:,k+1),~] = Newton_Raphson(fnc_x,x_evolution(:,k),num_iter_max,toll);
        
        if mod(k,20)==0
            disp([' integration step: ',num2str(k), ' of ',num2str(n_steps)]);
        end
    end
    
    x_timeseries = timeseries(x_evolution',time);

end

