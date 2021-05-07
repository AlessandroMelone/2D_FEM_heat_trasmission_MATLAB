%eq. 6.2.15
function [x,y] = map_gaussianPoint(csi,x1,y1,x3,y3)

    x = (x1+x3)/2 + (x3-x1)/2*csi;
    y = (y1+y3)/2 + (y3-y1)/2*csi;

end

