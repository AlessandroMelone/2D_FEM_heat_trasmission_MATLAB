function e_L2 = compute_errorNormL2(THETA_fem,THETA_analyticFunction, idxs_elements_full,nodes_coordinates_full,Lx,Ly)


    %% Integration

    %Gaussian quadrature for computing surface integrals
    [xy_gaussian_f, w_gaussian_f] = get_gaussianPoints_triangle();  %gaussian point and weight for a triangular domain
    n_gaussianPoints_f = length(xy_gaussian_f(:,1));                %n. of point to consider for the gaussian integration
    
    %%
    e_L2 = 0;
    Ne = 3;                 %n. of nodes of each element (triangle)
    n_elements = size(idxs_elements_full,2);
    for idx_element = 1:n_elements
        idxs_nodesElement = idxs_elements_full(1:Ne,idx_element);
        element_coordinates = nodes_coordinates_full(:,idxs_nodesElement);

        for idx_gaussianPoint = 1:n_gaussianPoints_f
            %Create the parametric variables in the gaussian point
            csi1 = xy_gaussian_f(idx_gaussianPoint,1);
            csi2 = xy_gaussian_f(idx_gaussianPoint,2);

            %Value of the shape functions in the idx_gaussianPoint point for the idx_element element
            [N_3T, grad_N_3T] = compute_N_3T(csi1,csi2);    %eq. 6.7.4

            %Passing from parent(parametric) space to physical space (eq.6.4.3, figure 6.20)
            x = element_coordinates(1,:)*N_3T';
            y = element_coordinates(2,:)*N_3T';

            J = compute_jacobian(grad_N_3T, element_coordinates);  %eq. 6.4.6
            det_J = det(J);
            
            THETA_xy = N_3T * THETA_fem(idxs_nodesElement);
            
            e_L2 = e_L2 + det_J*(THETA_xy - THETA_analyticFunction(x,y))^2 * w_gaussian_f(idx_gaussianPoint);   %eq. 18.2.2
        end
        %-----------------------------------------------------------------------------------%

    end
    
    e_L2 = (e_L2/(Lx*Ly))^(1/2);


end

