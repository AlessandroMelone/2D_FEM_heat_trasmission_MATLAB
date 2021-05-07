function [L,K,MM,f_omega,f_gammaR1,f_gammaQ] = solve_diffusiveEquation_FEM(D,k,Rtot,s,nodes_coordinates_full, idxs_boundaries_full, idxs_elements_full, robin_segments_label, neumann_segments_label)

    %% Problem
    % \div ( [D] \grad \theta(x,y)) + s(x,y) * mu(t) = \der \theta
    % q * n = - k \grad \theta * n_ = q0    on  \gamma_Q                      (Neumann Condition)
    % q * n = - k \grad \theta * n_ = - 1/Rtot*( T - \theta)   on  \gamma_R      (Robin Condition)

    %% Integration

    %Gaussian quadrature for computing surface integrals
    [xy_gaussian_f, w_gaussian_f] = get_gaussianPoints_triangle();  %gaussian point and weight for a triangular domain
    n_gaussianPoints_f = length(xy_gaussian_f(:,1));                %n. of point to consider for the gaussian integration


    %Gaussian quadrature for computing the line integrals
    %and robin condition (eq. 6.2.19)
    [x_gaussian_line, w_gaussian_line] = get_gaussianPoints_line();      
    n_gaussianPoints_line = length(x_gaussian_line(:,1));

    %% FEM
    Ne = 3;                 %n. of nodes of each element (triangle)
    ndof = size(nodes_coordinates_full,2);  %n. of nodes
    n_elements = size(idxs_elements_full,2);

    K = sparse(ndof,ndof);
    C = sparse(ndof,ndof);
    f = zeros(ndof,1);
    f_gammaR1 = zeros(ndof,1);
    f_omega = zeros(ndof,1);
    f_gammaQ = zeros(ndof,1);
    MM = sparse(ndof,ndof);

    %for each element
    for idx_element = 1:n_elements
        idxs_nodesElement = idxs_elements_full(1:Ne,idx_element);
        element_coordinates = nodes_coordinates_full(:,idxs_nodesElement);

        %Defining di gather matrix L_e, eq. 6.1.9
        L_e = sparse([1 2 3]', idxs_nodesElement,[1 1 1]', Ne, ndof);
    %   Equivalent to:
    %     L_e = sparse(Ne, ndof);   
    %     L_e(1, idxs_nodesElement(1)) = 1;
    %     L_e(2, idxs_nodesElement(2)) = 1;
    %     L_e(3, idxs_nodesElement(3)) = 1;

        %% Compute stiffness matrix K        (term 1_v)
        %-----------------------------------------------------------------------------------%
        K_e = compute_K_e(element_coordinates, D);  %eq. 6.2.14
        K = K + L_e' * K_e * L_e;                   %eq. 6.1.27
        %-----------------------------------------------------------------------------------%



        %% Compute line integrals
        % Compute terms associated to Robin BC: -int(w*theta)     (term 1_R2) 
        %                                       +int(N*T/k)      (term 1_R1)
        %                    and to Neumann BC: -int(w*q)   (term 1_q)
        %-----------------------------------------------------------------------------------%
        [M_e, A] = compute_M_e(element_coordinates); % needed to compute N_e
        C_e = zeros(Ne,Ne);
        f_gammaR1_e = zeros(Ne,1);
        f_gammaQ_e = zeros(Ne,1);    
        [sideNodes_coordinate,idx_boundary_segment] = get_sidePoints(idxs_nodesElement,idxs_boundaries_full(1:2,:), element_coordinates);

        for i = 1:size(sideNodes_coordinate,1)  %for each boundary segment of the triangle
            x1 = sideNodes_coordinate(i,1); y1 = sideNodes_coordinate(i,2);     % (x1,y1) and (x3,y3) should detect 2 nodes which side
            x3 = sideNodes_coordinate(i,3); y3 = sideNodes_coordinate(i,4);     % is contained in the boundary domain       

            for idx_gaussianPoint = 1:n_gaussianPoints_line
                csi = x_gaussian_line(idx_gaussianPoint);       %csi is a gaussian point on the 1D domain
                [x,y] = map_gaussianPoint(csi,x1,y1,x3,y3);    %map csi into xy variable, eq. 6.2.15

                N_e = compute_N_e(x,y,M_e,A);   %eq. 6.2.10 or eq. 6.2.16

                l_31 = sqrt((x3-x1)^2 + (y3-y1)^2);                 %eq. 6.2.18

                %idxs_boundaries_full(5,idxs_boundary_segment(i)) contains the label of the boundary segment which the side of the triangle belong
                if(find(idxs_boundaries_full(5,idx_boundary_segment(i)) == robin_segments_label))
                    C_e = C_e + D/(Rtot*k) * N_e' * N_e * l_31/2 * w_gaussian_line(idx_gaussianPoint);    %for term 1_R2
                    %check constants and signs
                    f_gammaR1_e = f_gammaR1_e + D/(Rtot*k) * N_e'* 1 * l_31/2 * w_gaussian_line(idx_gaussianPoint); %for term 1_R1 (eq. 6.2.19)
                else
                    %segment belong to the neumann boundary
                    f_gammaQ_e = f_gammaQ_e - D/k * N_e' * 1 * l_31/2 * w_gaussian_line(idx_gaussianPoint); %for term 1_q
                end
            end
            %display(['Element: ', num2str(idx_element), ' . side: p1=[',num2str(x1),' ',num2str(y1),']',' p3=[',num2str(x3),' ',num2str(y3),']','  ROBIN']);

        end

        C = C + L_e' * C_e * L_e;
        %-----------------------------------------------------------------------------------%



        %% Compute surface integrals
        %%Compute term associated to <w,s> (forcement)      (term 2)
        %                     and to <?,?> (time derivative) (term 3)
        %-----------------------------------------------------------------------------------%
        f_omega_e = zeros(3,1);
        MM_e = zeros(Ne,Ne);
        for idx_gaussianPoint = 1:n_gaussianPoints_f
            %Create the parametric variables in the gaussian point
            csi1 = xy_gaussian_f(idx_gaussianPoint,1);
            csi2 = xy_gaussian_f(idx_gaussianPoint,2);

            %Value of the shape functions in the idx_gaussianPoint point for the idx_element element
            [N_3T, grad_N_3T] = compute_N_3T(csi1,csi2);    %eq. 6.7.4

            %grad_N_3T = [d N1 / d csi1,     d N1 / d csi2;
            %             d N2 / d csi1,     d N2 / d csi2;  
            %             d N3 / d csi1,     d N3 / d csi2;  

            %Passing from parent(parametric) space to physical space (eq.6.4.3, figure 6.20)
            x = element_coordinates(1,:)*N_3T';
            y = element_coordinates(2,:)*N_3T';

            J = compute_jacobian(grad_N_3T, element_coordinates);  %eq. 6.4.6
            det_J = det(J);
            
            f_omega_e = f_omega_e + det_J * N_3T' * s(x,y) * w_gaussian_f(idx_gaussianPoint); %eq. 6.7.5
            MM_e = MM_e + det_J * (N_3T'*N_3T) * w_gaussian_f(idx_gaussianPoint);         %mass matrix for the element 
        end

        MM = MM + L_e' * MM_e * L_e;
        %-----------------------------------------------------------------------------------%

        %% assembly f
        f_e = f_gammaR1_e + f_omega_e + f_gammaQ_e;      %eq. 6.1.21b
        f_gammaR1 = f_gammaR1 + L_e' * f_gammaR1_e;
        f_omega = f_omega + L_e' * f_omega_e;
        f_gammaQ = f_gammaQ + L_e' * f_gammaQ_e;
        f = f + L_e' * f_e;             %eq. 6.1.28

    end
    
L = K + C;
end
