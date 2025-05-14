function prob = vonMises(x,mu,k)

prob = exp(k.*cos(x-mu))/(2*pi*besseli(0,k));

end