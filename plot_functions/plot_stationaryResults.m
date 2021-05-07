function [] = plot_stationaryResults(THETA_fem_vector, nodes_coordinates_full)
    x=nodes_coordinates_full(1,:); %1x28
    y=nodes_coordinates_full(2,:); %1x28; %1x28
    z=THETA_fem_vector'; %1x28
    xq = linspace(min(x), max(x));
    yq = linspace(min(y), max(y));
    [X,Y] = meshgrid(linspace(min(x), max (x)),linspace(min(y), max (y)));
    THETA_fem_grid = griddata(x,y,z, X, Y, 'cubic');
    figure;

    %Soluzione FEM
    subplot(1,2,1), surf(X,Y,THETA_fem_grid), shading interp; hold on, axis square;
    contour3(X,Y,THETA_fem_grid,linspace(min(THETA_fem_grid(:)),max(THETA_fem_grid(:)),40),'k'); 
    xlabel('x'); ylabel('y'); zlabel('THETA approx (fem)');
    title('FEM Solution')
    subplot(1,2,2); contour(X,Y,THETA_fem_grid); axis square
    xlabel('x'); ylabel('y'); zlabel('THETA approx (fem)');

end

