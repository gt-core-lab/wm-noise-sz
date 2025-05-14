clear all;
close all;

id = 'test';

whichmodel = 2; 
% 1: von mises + uniform 
% 2: von mises only 

xRange = [-90 90]; % will need to figure this out!! 


fname = strcat(id, '_VWMdat.mat'); 
load(fname); 

errInd = 6;
nsets = 3;
nparams = 2; 

set1 = matfile(matfile(:,1)==1, errInd); 
set2 = matfile(matfile(:,1)==2, errInd); 
set3 = matfile(matfile(:,1)==4, errInd); 

setR{1} = set1 .* 2*pi/diff(xRange); 
setR{2} = set2 .* 2*pi/diff(xRange); 
setR{3} = set3 .* 2*pi/diff(xRange);

for i = 1:nsets
    
    data = setR{i}; 

    params0 = [0.3 0.01]; 

    LB = [0.00001 0];
    UB = [1000 1];

    options = optimset('MaxIter',1000000, 'MaxFunEvals', 1000000); 

    [params_fit(i,:), fval(i)] = fmincon(@(params) get_likeli(params,data,whichmodel),params0,[],[],[],[],LB,UB,[],options); 
    % [params_fit fval] = fminsearch(@(params) get_likeli(params,data),params0,options);
    
    mu = 0; % mean of von mises (we assume that the mean of the error is 0)
    k = params_fit(i,1); % width of von mises: smaller the wider!!! (can be directly used as 'precision'); 
    lambda = params_fit(i,2); % guessing rate 
    
    plot_x = xRange(1):xRange(2);
    plot_xR = plot_x .* 2*pi/diff(xRange);
    
    if whichmodel == 1
        prob = (1-lambda)*vonMises(plot_xR,mu,k) + lambda*(1/(2*pi));
    elseif whichmodel == 2
        prob = vonMises(plot_xR,mu,k);
    end
    
    [f,x] = hist(data,10);
    
    subplot(1,nsets,i);
    hold on;
    bar(x,f/trapz(x,f));
    plot(plot_xR,prob,'r-'); 
    set(gca,'xtick',-pi:pi:pi);
    set(gca,'xticklabel',{'-pi','0','pi'})

end

 





