%eq. 6.2.4, eq. 6.2.11
function [M_e, A_e] = compute_M_e(element_coordinates)
    %Compute the values of the 3 shape function in the specified (x,y)
    %coordinate for the specific element
    x1 = element_coordinates(1,1);
    x2 = element_coordinates(1,2);
    x3 = element_coordinates(1,3);
    
    y1 = element_coordinates(2,1);
    y2 = element_coordinates(2,2);
    y3 = element_coordinates(2,3);
    

    M_e = [1, x1 ,y1;         %eq. 6.2.4
        1, x2, y2; 
        1, x3, y3];     
    
    A_e = 1/2 * det(M_e);       %eq. 6.2.11



end

