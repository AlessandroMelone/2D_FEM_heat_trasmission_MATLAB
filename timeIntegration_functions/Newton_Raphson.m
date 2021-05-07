function [zero,niter]=Newton_Raphson(fnc_x,x0,num_iter_max,toll)
% Newton-Raphson method for system f(x)=0; 

    x_star = x0;%initialize x* s.t. f(x*)=0    
    niter = 0; 
    error = toll+1; %to start the loop
    while (error >= toll) && (niter <= num_iter_max)
        niter = niter+1;
        %compute x(k+1)-x(k)
        diff = -jacobiancsd(fnc_x,x_star) \ fnc_x(x_star);
        %diff = -jacobianest(fnc_x,x_star)\fnc_x(x_star);

        %compute x(k+1)
        x_star = x_star + diff;
        error = norm(diff);
        
        %Debug
%         disp([' iteration: ',num2str(niter),'; norm error: ',num2str(error), ...
%             '; f(x*): ',num2str(norm(fnc_x(x_star)))]);
    end
    zero = x_star;
    
    return

