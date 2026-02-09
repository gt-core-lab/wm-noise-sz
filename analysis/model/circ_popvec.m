function th = circ_popvec(phi, r)
    x = sum(r(:).*cos(phi(:))); y = sum(r(:).*sin(phi(:)));
    th = atan2(y, x);
end