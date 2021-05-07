function index = find_nodeIndex(coordinates,nodes_coordinates_full)
    distance(1,:) = abs(nodes_coordinates_full(1,:) - coordinates(1));
    distance(2,:) = abs(nodes_coordinates_full(2,:) - coordinates(2));
    
    k = abs(distance(1,:) + distance(2,:));
    [value, index] = min(k);

end

