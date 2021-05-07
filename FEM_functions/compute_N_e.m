%eq. 6.2.10
function N_e = compute_N_e(x,y,M_e,A_e)
    %Compute the values of the 3 shape function in the specified (x,y)
    %coordinate for the specific element
    x1 = M_e(1,2);    
    x2 = M_e(2,2);
    x3 = M_e(3,2);
    y1 = M_e(1,3);    
    y2 = M_e(2,3);
    y3 = M_e(3,3);

    %eq. 6.2.10
    N1_e = 1/(2*A_e)*(x2*y3 - x3*y2 + (y2-y3)*x + (x3-x2)*y); %value of N1_e in xy
    N2_e = 1/(2*A_e)*(x3*y1 - x1*y3 + (y3-y1)*x + (x1-x3)*y);
    N3_e = 1/(2*A_e)*(x1*y2 - x2*y1 + (y1-y2)*x + (x2-x1)*y);
    
    N_e = [N1_e, N2_e, N3_e];

end

