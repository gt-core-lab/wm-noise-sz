function log_pContrast= predictedcontrast(parameters,Ne,d)
r=parameters(1);
beta=parameters(2);
Na=parameters(3);
Nm=parameters(4);
Am=parameters(5);
Af=parameters(6);
Aa=parameters(7);

%Ne
%d

%pContrast=((((1+(Am*Nm)^2)*((Af*Ne)^(2*r))+(Aa*Na)^2)./(1/(d*d)-(Am*Nm)^2)).^(1/(2*r)))/beta;

%avoid complex number, so we do regularization here
x=1/(d^2)-(Am*Nm)^2;
if x<0
    x=0.01;
end

log_pContrast=1/(2*r)*log((1+(Am*Nm)^2)*((Af*Ne).^(2*r))+(Aa*Na)^2)-1/(2*r)*log(x)-log(beta);

%log_pContrast

if isnan(log_pContrast)
    save('data_temp','parameters','log_pContrast');
    error('Error: all contrast should be a value');
end


if ~isreal(log_pContrast)
    save('data_temp','parameters','log_pContrast');
    error('Error: all contrast should be a real value');
end


end
