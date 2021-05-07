function [] = plot_stationaryResults_analyticSolution(THETA_fem_vector, THETA_analytic_vector, d_approx_anaylitic_vector, nodes_coordinates_full)
    x=nodes_coordinates_full(1,:); %1x28
    y=nodes_coordinates_full(2,:); %1x28; %1x28
    z=THETA_fem_vector'; %1x28
    xq = linspace(min(x), max(x));
    yq = linspace(min(y), max(y));
    [X,Y] = meshgrid(linspace(min(x), max (x)),linspace(min(y), max (y)));
    THETA_fem_grid = griddata(x,y,z, X, Y, 'cubic');
    z = THETA_analytic_vector';
    THETA_analytic_grid = griddata(x,y,z, X, Y, 'cubic');
    z = d_approx_anaylitic_vector';
    d_analytic_approx_grid = abs(THETA_analytic_grid - THETA_fem_grid);


    figure;

    %Soluzione FEM
    subplot(3,2,1), surf(X,Y,THETA_fem_grid), shading interp; hold on, axis square;
    contour3(X,Y,THETA_fem_grid,linspace(min(THETA_fem_grid(:)),max(THETA_fem_grid(:)),40),'k'); 
    xlabel('x'); ylabel('y'); zlabel('THETA approx (fem)');
    title('FEM Solution')
    subplot(3,2,2); contour(X,Y,THETA_fem_grid); axis square
    xlabel('x'); ylabel('y'); zlabel('THETA approx (fem)');

    %Soluzione analitica
    subplot(3,2,3), surf(X,Y,THETA_analytic_grid), shading interp; hold on, axis square;
    contour3(X,Y,THETA_analytic_grid,linspace(min(THETA_analytic_grid(:)),max(THETA_analytic_grid(:)),40),'k'); 
    xlabel('x'); ylabel('y'); zlabel('THETA analytic');
    title('Analytic solution')
    subplot(3,2,4); contour(X,Y,THETA_analytic_grid); axis square
    xlabel('x'); ylabel('y'); zlabel('THETA analyitc');

    %Difference
    subplot(3,2,5), surf(X,Y,d_analytic_approx_grid), shading interp; hold on, axis square;
    contour3(X,Y,d_analytic_approx_grid,linspace(min(d_analytic_approx_grid(:)),max(d_analytic_approx_grid(:)),40),'k'); 
    xlabel('x'); ylabel('y'); zlabel('Difference');
    title('Difference')
    subplot(3,2,6); contour(X,Y,d_analytic_approx_grid); axis square
    xlabel('x'); ylabel('y'); zlabel('Difference');


end

