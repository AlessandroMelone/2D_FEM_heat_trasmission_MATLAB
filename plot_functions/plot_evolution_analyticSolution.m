function [] = plot_evolution_analyticSolution(z_fem_evolution, z_analytic_evolution,difference_evolution,time, nodes_coordinates_full, Lx,Ly, NN)
    figure;
    set(gcf,'Position',[100 100 800 800])
    
    
    n_steps = length(time);
    

    
    x=nodes_coordinates_full(1,:);
    y=nodes_coordinates_full(2,:);
    
    
    for k = [1:NN:n_steps, n_steps]
        z_fem=z_fem_evolution(k,:);
        z_analytic = z_analytic_evolution(k,:);
        difference = difference_evolution(k,:);
        [X,Y] = meshgrid(linspace(min(x), max (x)),linspace(min(y), max (y)));
        z_fem_grid = griddata(x,y,z_fem, X, Y, 'cubic');
        z_analytic_grid = griddata(x,y,z_analytic,X,Y,'cubic');
        difference_grid = griddata(x,y,difference,X,Y,'cubic');


        %FEM solution
        subplot(3,2,1), surf(X,Y,z_fem_grid),
        axis([0 Lx 0 Ly min(min(z_fem_evolution)) .1+max(max(z_fem_evolution))]);shading interp; axis square;
        %contour3(X,Y,THETA_fem_grid,linspace(min(THETA_fem_grid(:)),max(THETA_fem_grid(:)),40),'k'); 
        xlabel('x'); ylabel('y'); zlabel('THETA approx (fem)');
        title(['FEM Solution,  step: ',num2str(k), ',  time = ', num2str(time(k)),'s'])
        subplot(3,2,2); contour(X,Y,z_fem_grid); axis square
        xlabel('x'); ylabel('y'); zlabel('THETA approx (fem)');
        
        
        %Analytic solution
        subplot(3,2,3), surf(X,Y,z_analytic_grid),
        axis([0 Lx 0 Ly min(min(z_analytic_evolution)) .1+max(max(z_analytic_evolution))]);shading interp; axis square;
        %contour3(X,Y,THETA_fem_grid,linspace(min(THETA_fem_grid(:)),max(THETA_fem_grid(:)),40),'k'); 
        xlabel('x'); ylabel('y'); zlabel('THETA analytic');
        title(['Analytic Solution,  step: ',num2str(k), ',  time = ', num2str(time(k)),'s'])
        subplot(3,2,4); contour(X,Y,z_analytic_grid); axis square
        xlabel('x'); ylabel('y'); zlabel('THETA analytic');
        
        
        %Difference
        subplot(3,2,5), surf(X,Y,difference_grid),
        axis([0 Lx 0 Ly min(min(difference_evolution)) .1+max(max(difference_evolution))]); shading interp; axis square;
        %contour3(X,Y,THETA_fem_grid,linspace(min(THETA_fem_grid(:)),max(THETA_fem_grid(:)),40),'k'); 
        xlabel('x'); ylabel('y'); zlabel('Difference');
        title(['Difference,  step: ',num2str(k), ',  time = ', num2str(time(k)),'s'])
        subplot(3,2,6); contour(X,Y,difference_grid); axis square
        xlabel('x'); ylabel('y'); zlabel('Difference');
       
        pause(0.0001);
    end
end

