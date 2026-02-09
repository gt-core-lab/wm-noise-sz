function th_out = add_circ_gauss(th_in, sigma)
% Add i.i.d. Gaussian noise with SD = sigma (radians), wrap to [0, 2π)
    if sigma <= 0
        th_out = mod(th_in, 2*pi);
    else
        th_out = mod(th_in + sigma .* randn(size(th_in)), 2*pi);
    end
end