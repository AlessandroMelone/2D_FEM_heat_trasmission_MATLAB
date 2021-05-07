function out = radiator_forcement(x,y)
    
    global SOURCE_X_MIN SOURCE_X_MAX SOURCE_Y_MIN SOURCE_Y_MAX
    global SOURCE_X_MIN_2 SOURCE_X_MAX_2 SOURCE_Y_MIN_2 SOURCE_Y_MAX_2
    
    if (x >= SOURCE_X_MIN) && (x <= SOURCE_X_MAX) && (y >= SOURCE_Y_MIN) && (y <= SOURCE_Y_MAX)
        out = 1;
    else
        if (x >= SOURCE_X_MIN_2) && (x <= SOURCE_X_MAX_2) && (y >= SOURCE_Y_MIN_2) && (y <= SOURCE_Y_MAX_2)
            out = 1;
        else
            out = 0;
        end
    
    end


    
end

