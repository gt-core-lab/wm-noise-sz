function k = vm_kappa_from_R(R)
    if R < 0.53
        k = 2*R + R^3 + (5*R^5)/6;
    elseif R < 0.85
        k = -0.4 + 1.39*R + 0.43/(1 - R);
    else
        k = 1/(R^3 - 4*R^2 + 3*R);
    end
    k = max(k, 1e-6);
end