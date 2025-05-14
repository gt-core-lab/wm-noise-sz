%this model is to fit 2 curves to see the difference between 2 curves

%input actual data here
function [RSquare,paramsfit1,paramsfit2] = PTM(data,varargin) 

%% Perceptual template model for contrast threshold data.
% Instruction for input parameters

% data:  data should be a conditions * noise level matrix which contains
% all contrast threshold data between (0,1).If it's not contrast threshold
% data, please normlize data to range (0,1);
% Please range data maxtrix like [baseline_70%;expCondition_70%;baseline_79%;expCondition_79%]
% baseline_70% refers contrast thresholds at baseline condition at 70% performance level
% similarly, it can model perceptual learning as [pre_70%;post_70%;pre_79%;post_79%];

% varargin specify which parameters you want to change. You can input like
% PTM(data,'Aa'), means you specify the model that assumes only internal
% noise change between baseline and experiment condition.You can input
% like PTM(data,'Am','Ae','Aa') to derive full model which assumes all
% three parameters change. The default will be model with no parameters
% change, namely the least model.

%%
% you can use test Data here
testData=[0.15967725	0.1736934	0.15592305	0.18593935	0.22270735	0.2548816	0.40217735	0.49443375
0.1376227	0.1344836	0.1378777	0.15789365	0.18596365	0.21048745	0.30942795	0.393189
0.22438295	0.2210381	0.2178485	0.2640276	0.3300025	0.37167855	0.5111114	0.7147097
0.1736204	0.1707629	0.1804101	0.23623375	0.2704351	0.30184185	0.4042601	0.49448005];

%varargin={'Am','Ae','Aa'};

%if nargin<2 ||isempty(varargin), varargin=cell(3,1);end %if no reduced model was specified, We assume PTM fit least change model
if nargin<1 ||isempty(data), data=testData;end
        
close all;
p.modelOption=cell(3,1);
%sort input, we create p.modelOption if p.modelOption = {[],[],'Aa'}, this
%refers to Aa change model
if ~isempty(varargin)
    for i=1:length(varargin)
        switch varargin{i}
            case 'Am'
                p.modelOption{1} ='Am';
            case 'Ae'
                p.modelOption{2} ='Ae';
            case 'Aa'
                p.modelOption{3} ='Aa';
        end
    end
end
%
p.d                        =[1.089 1.634]';   %d' for 70% and 79%   
% p.Ne                       =[0 0.003 0.0061 0.0124 0.0251 0.051 0.103 0.21]';%you can change this based on your experiment 
p.Ne = [0 0.01 0.0166 0.0276 0.0458 0.0761 0.1264 0.21];



%initial parameter
params0=[1.3
    1.2
    0.1
    0.2
];


% 1.1 0.7 0.01 0.04 for data easy
% 1.3 0.7 0.01 0.1 for data easy no outlier

params0=[params0; ones(length(varargin),1)*2];

paramsLower=[0.0001; 0.00001; 0.00001;0.00001;0.00001;0.00001; 0.00001];
paramsUpper=[5; 5; 0.5;0.8;10;10; 10];

%Then do fitting
options=optimset('MaxFunEvals',100000,'MaxIter',100000);
%paramsfit=fminsearch(@(params) computeSSE(params,Actual_Data),params0);
% [params_temp sse]=fminsearch(@(params) costfun(params,data,p),params0,options);
[params_temp sse]=fmincon(@(params) costfun_fmincon(params,data,p),params0,[],[],[],[],paramsLower,paramsUpper,[],options);
paramsfit = ones(7,1);paramsfit(1:4)=params_temp(1:4);
count = 5;
%Determine parameters
for i=1:length(p.modelOption)
    if isempty(p.modelOption{i})
        paramsfit(count+i-1)=1;
    else
        switch p.modelOption{i}
            case 'Am'
                paramsfit(4+i)=params_temp(count);
                count=count+1;
            case 'Ae'
                paramsfit(4+i)=params_temp(count);
                count=count+1;
            case 'Aa'
                paramsfit(4+i)=params_temp(count);
                count=count+1;
        end
    end
end

%show best fiting parameter, the sequence is [r;beta;Na;Nm;Am;Ae;Aa];
fprintf('The fitted model option is: ');
p.modelOption
fprintf(['\n The parameters for fitted model are :' ...
    '\n r:   ', num2str(paramsfit(1))...
    '\n beta:', num2str(paramsfit(2))...
    '\n Na:  ', num2str(paramsfit(3))...
    '\n Nm:  ', num2str(paramsfit(4))...
    '\n Am:  ', num2str(paramsfit(5)),', refers to ',num2str(paramsfit(5)*100),'%% Multiplicative noise than baseline'...
    '\n Ae:  ', num2str(paramsfit(6)),', refers to ',num2str(paramsfit(6)*100),'%% External noise than baseline'...
    '\n Aa:  ', num2str(paramsfit(7)),', refers to ',num2str(paramsfit(7)*100),'%% Internal noise than baseline'...
    '\n The fitted model can account for ', num2str((1-sse)*100),'%% variance of data \n']);

%show RSqure for model fitting. Closer to 1, better the fitting.
RSquare=1-sse;



%visualize data

paramsfit1 = [paramsfit(1:4); 1; 1; 1]';
paramsfit2 = paramsfit';

ypredict11 = exp(predictedcontrast(paramsfit1,p.Ne, p.d(1)));
ypredict12 = exp(predictedcontrast(paramsfit1,p.Ne, p.d(2)));
ypredict21 = exp(predictedcontrast(paramsfit2,p.Ne, p.d(1)));
ypredict22 = exp(predictedcontrast(paramsfit2,p.Ne, p.d(2)));

Nlevels = p.Ne';
Nlevels(1) = 0.00001;

axismat=[0.00001 1.5 0.02 1.5];
subplot(1,2,1);
loglog(Nlevels,ypredict11,'r'); hold on;
loglog(Nlevels,ypredict21,'b'); hold on;
loglog(Nlevels,10.^(data(1,:))./100,'ro'); hold on;
loglog(Nlevels,10.^(data(2,:))./100,'bo'); hold on;
axis(axismat);
subplot(1,2,2);
loglog(Nlevels,ypredict12,'r'); hold on;
loglog(Nlevels,ypredict22,'b'); hold on;
loglog(Nlevels,10.^(data(3,:))./100,'ro'); hold on;
loglog(Nlevels,10.^(data(4,:))./100,'bo'); hold on;
axis(axismat);

sse

