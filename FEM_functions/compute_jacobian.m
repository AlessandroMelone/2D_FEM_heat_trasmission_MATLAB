%eq. 6.4.6
function J = compute_jacobian(grad_N, element_coordinates)

    J11 = 0;    J12 = 0;    J21 = 0;    J22 = 0;
    for i = 1:3
        J11 = J11 + grad_N(i,1) * element_coordinates(1,i);
        J21 = J21 + grad_N(i,2) * element_coordinates(1,i);
        J12 = J12 + grad_N(i,1) * element_coordinates(2,i);
        J22 = J22 + grad_N(i,2) * element_coordinates(2,i);
    end

    J = [J11, J12;
        J21, J22];
end


