function [xy_gaussian, w_gaussian] = get_gaussianPoints_triangle()
    %order 3
    xy_gaussian = [1/3 1/3;
                    1/5 1/5;
                    3/5 1/5;
                    1/5 3/5];
                
    w_gaussian = [-27/96;
                    25/96;
                    25/96;
                    25/96];
    
end

