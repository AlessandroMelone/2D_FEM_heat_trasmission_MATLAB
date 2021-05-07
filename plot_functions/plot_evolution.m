function [] = plot_evolution(z_fem_evolution,time, nodes_coordinates_full, Lx,Ly, NN)
    figure;
    
    n_steps = length(time);
    
    x=nodes_coordinates_full(1,:);
    y=nodes_coordinates_full(2,:);
    
    for k = [1:NN:n_steps, n_steps]
        z_fem=z_fem_evolution(k,:);
        [X,Y] = meshgrid(linspace(min(x), max (x)),linspace(min(y), max (y)));
        z_fem_grid = griddata(x,y,z_fem, X, Y, 'cubic');

        %FEM solution
        subplot(1,2,1), surf(X,Y,z_fem_grid),
        axis([0 Lx 0 Ly min(min(z_fem_evolution)) .1+max(max(z_fem_evolution))]);shading interp; axis square;
        %contour3(X,Y,THETA_fem_grid,linspace(min(THETA_fem_grid(:)),max(THETA_fem_grid(:)),40),'k'); 
        xlabel('x'); ylabel('y'); zlabel('THETA approx (fem)');
        title(['FEM Solution,  step: ',num2str(k), ',  hours = ', num2str(time(k)/3600),'h'])
        subplot(1,2,2); contour(X,Y,z_fem_grid); axis square
        xlabel('x'); ylabel('y'); zlabel('THETA approx (fem)');

        pause(0.01);
    end



end
