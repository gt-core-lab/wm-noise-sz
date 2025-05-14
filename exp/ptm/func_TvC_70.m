function pContrast=func_TvC_70(parameters,Ne)
r=parameters{1};
beta=parameters{2};
Na=parameters{3};
Nm=parameters{4};
Am=1;
Af=1;
Aa=1;
%
d   =1.089;    %d' for 70%
%d    =1.3490;   %d' for 75%  
%d   =1.634;    %d' for 79%
 

pContrast=((((1+(Am*Nm).^2).*((Af*Ne).^(2*r))+(Aa*Na).^2)./(1/(d*d)-(Am.*Nm).^2)).^(1/(r*2)))./beta;
end