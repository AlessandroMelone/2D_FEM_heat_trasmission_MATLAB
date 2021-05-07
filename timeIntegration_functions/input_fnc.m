function fnc = input_fnc(x,u)
    global A_lti B_lti
    fnc = A_lti * x + B_lti*u;
end

