function [sideNodes_coordinate,idxs_boundary_segment] = get_sidePoints(idxs_nodesElement,idxs_boundary_c,element_coordinates)
%     get_sidePoints: check if a point of the triangle belong to the
%     boundary
%     INPUT:
%         idxs_nodesElement: number of the index of the triangle
%         idxs_boundary_c: [2 x num_boundary_segments] each column contain 
%         the index of a boundary segment 
%         element_coordinates: coords of the tre vertices of the triangle
%     OUTPUT:
%         sideNodes_coordinate: [number of boundary segments x 4] at i-th
%         row contain sideNodes_coordinate(i,:)=[x1 y1 x2 y2] with (x1,y1)
%         coords of one point of the segment and (x2,y2) coords of the
%         other point of the segment
%
%         idxs_boundary_segment: array size=number of boundary segments,
%         each element is the index of the boundary segment in the 
%         idxs_boundary_c array
 

    %sideNodes_coordinate = [x11 y11 x13 y13;   <- first point
    %                        x21 y21 x23 y23];
    sideNodes_coordinate = [];
    idxs_boundary_segment = [];
    
    n_sides = 0;
    for i = 1:3 
        j = i+1; % at first iteration (i,j)=(1,2) at the second: (2,3); third: (3,1); 
        if j==4 
            j = 1;
        end
        
        %check if (i,j) or (j,i) belong to the boundary segment array (idxs_boundary_c)
        side = [idxs_nodesElement(i); idxs_nodesElement(j)]; %index segment (i,j)
        A1 = side == idxs_boundary_c; %A1 sparse?
        B1 = and(A1(1,:),A1(2,:)); %B1(i)!=0 if i-th segment of idxs_boundary_c(:,i) belong to the boundary 
        %B1 sparse?
        
        side = [idxs_nodesElement(j); idxs_nodesElement(i)]; %index segment (j,i)
        A2 = side == idxs_boundary_c; 
        B2 = and(A2(1,:),A2(2,:));
        
        B = or(B1,B2);
        idx = find(B);
        if idx
            n_sides = n_sides + 1;
            idxs_boundary_segment(n_sides) = idx; 
            sideNodes_coordinate(n_sides,:) = [element_coordinates(:,i)', element_coordinates(:,j)'];
        end

    end
    
    
end

