function loglikeli = get_likeli(params,data,whichmodel)

% degree to radian (xRange= -180 to 180; typically) 
% 2*pi/diff(p.xRange)

% Von Mises function 
% prob = exp(k*cos(x-mu))/(2*pi*besseli(0,k));

% Model (von mises + uniform) 
% prob = (1-lambda)*vonMises(x,mu,k) + lambda*(1/(2*pi)); 

mu = 0; 
k = params(1);
lambda = params(2); 
x = data; 

if whichmodel == 1
    likeli = (1-lambda).*vonMises(x,mu,k) + lambda*(1/(2*pi));
elseif whichmodel == 2
    likeli = vonMises(x,mu,k);
end

loglikeli = -sum(log(likeli)); 

end



