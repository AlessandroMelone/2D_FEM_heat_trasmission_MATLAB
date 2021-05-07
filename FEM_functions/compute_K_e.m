%eq. 6.2.14

function K_e = compute_K_e(element_coordinates,D)
    x1 = element_coordinates(1,1);
    x2 = element_coordinates(1,2);
    x3 = element_coordinates(1,3);
    
    y1 = element_coordinates(2,1);
    y2 = element_coordinates(2,2);
    y3 = element_coordinates(2,3);
    

    M = [1, x1 ,y1;         %eq. 6.2.4
        1, x2, y2; 
        1, x3, y3];     
    
    A = 1/2 * det(M);                       %eq. 6.2.11
        
    B_e = 1/(2*A) * [y2-y3, y3-y1, y1-y2;       %eq. 6.2.13
                    x3-x2, x1-x3, x2-x1];

    K_e = B_e' * D * B_e * A;      %eq. 6.2.14
end