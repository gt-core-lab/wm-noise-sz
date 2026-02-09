function beta1 = deming_slope(x, y, lambda)
    x = x(:); y = y(:);
    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    if numel(x) < 3, beta1 = NaN; return; end
    xbar = mean(x); ybar = mean(y);
    Sxx  = mean((x - xbar).^2);
    Syy  = mean((y - ybar).^2);
    Sxy  = mean((x - xbar).*(y - ybar));
    D = (Syy - lambda*Sxx)^2 + 4*lambda*Sxy^2;
    if ~isfinite(D) || Sxy == 0, beta1 = NaN; return; end
    beta1 = (Syy - lambda*Sxx + sqrt(D)) / (2*Sxy);
end