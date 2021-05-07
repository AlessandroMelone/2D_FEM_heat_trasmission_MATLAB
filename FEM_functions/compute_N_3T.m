%eq. 6.7.4
function [N_3T,grad_N_3T] = compute_N_3T(csi1,csi2)
    N_3T(1) = csi1;
    N_3T(2) = csi2;
    N_3T(3) = 1-csi1-csi2;

    grad_N_3T(1,:) = [1, 0];
    grad_N_3T(2,:) = [0, 1];
    grad_N_3T(3,:) = [-1, -1];
    
end

