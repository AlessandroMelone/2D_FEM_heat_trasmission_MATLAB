function out = controller(input)
    global RELE_THRESHOLD
    reference = input(1);   sensor = input(2);  status = input(3);
    
    
    if sensor >= reference + RELE_THRESHOLD
        status = 0;
    elseif  sensor <= reference - RELE_THRESHOLD
        status = 1;
    end
    
    if status
        u = 1;
    else
        u = 0;
    end
    
    out = [u; status];
    
    %%P controller
    %error = reference - sensor;
    %u = 100*error;
end

