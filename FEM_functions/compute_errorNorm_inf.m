function errorNorm_inf = compute_errorNorm_inf(THETA_fem, THETA_analytic)
    errorNorm_inf = max(THETA_fem - THETA_analytic);
end

