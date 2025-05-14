function [stats, samples] = cluster_mPTM9_JAGS_2(condition)
    
    % fixed beta (=1), Nm (= 0.6) and r (=2.64)
    
    ngroups = 2;
    ntotalchains = 20;
    emat = expmat(1:ngroups, 1:ntotalchains);
    conds = [1:length(emat)]';
    emat = [conds emat];
    
    setnewseedmat = round(rand(1,length(emat))*1000000);
    setnewseed = setnewseedmat(condition);
    rng(setnewseed);

    idgroup = emat(condition,2);
    if idgroup == 1
        whichgroup = 'SCZ';
    elseif idgroup == 2
        whichgroup = 'CON';
    end
    whichchain = emat(condition,3);
    
    % MCMC parameters
    nchains = 1; 
    nburnin = 8000; % 15000
    nsamples = 1000; % 1000
    nthin = 100; % 200
    doparallel = 0; 
    
    loadname = strcat('stacked_PTM_',whichgroup,'.mat');
    savename = strcat('stats_PTM9_',whichgroup,'_',num2str(whichchain),'.mat');
    workingdir = strcat('PTM9_',whichgroup,'_',num2str(whichchain));
    
    % data
    load(loadname);
    data = data_all;
    nsubjs = length(data);
    ntrials = 40;
    nnoise = 8;
    naccuracy = 2;
    d = [1.089; 1.634];
    pthreshold = [0.7071 0.7937];
    Ne = unique(data{1}(:,1)); 
    x1 = zeros(nsubjs, nnoise, ntrials);
    resp1 = zeros(nsubjs, nnoise, ntrials);
    x2 = zeros(nsubjs, nnoise, ntrials);
    resp2 = zeros(nsubjs, nnoise, ntrials);
    for i = 1:nsubjs
        temp = data{i};
        for j = 1:naccuracy % accuracy level 
            temp2 = temp(temp(:,4)==j,:);
            for m = 1:nnoise
                if j == 1
                    x1(i,m,:) = temp2(temp2(:,1)==Ne(m),2); % contrast
                    resp1(i,m,:) = temp2(temp2(:,1)==Ne(m),3); % response
                elseif j == 2
                    x2(i,m,:) = temp2(temp2(:,1)==Ne(m),2); % contrast
                    resp2(i,m,:) = temp2(temp2(:,1)==Ne(m),3); % response
                end
            end
        end
    end
    
    x1(x1<0.01) = 0.01;
    x1 = log(x1*1000);
    x2(x2<0.01) = 0.01;
    x2 = log(x2*1000);
    
    % assign matlab variables to the observed nodes 
    datastruct = struct('nsubjs',nsubjs,'ntrials',ntrials, ...
        'd',d,'Ne',Ne,'x1',x1,'resp1',resp1,'x2',x2,'resp2',resp2,'pthreshold',pthreshold,'nnoise',nnoise);
    
    % assign initial values 
    for i = 1:nchains 
%         S.mur = randn(1)*0.1+2;
%         S.mubeta = randn(1)*0.001+0.5;
%         S.muNm = randn(1)*0.001+0.3;
        S.muNa = randn(1)*0.00001+0.005;
        S.muAf = randn(1)*0.001+1.2;
        S.muslope = randn(1)*0.1+5;
%         S.sigmar = 0.1;
%         S.sigmabeta = 0.001;
%         S.sigmaNm = 0.001;
        S.sigmaNa = 0.00001;
        S.sigmaAf = 0.001;
        S.sigmaslope = 0.1;
%         S.r = randn(1,nsubjs)*S.sigmar + S.mur;
%         S.beta = randn(1,nsubjs)*S.sigmabeta + S.mubeta;
%         S.Nm = randn(1,nsubjs)*S.sigmaNm + S.muNm;
        S.Na = randn(1,nsubjs)*S.sigmaNa + S.muNa;
        S.Af = randn(1,nsubjs)*S.sigmaAf + S.muAf;
        S.slope = randn(1,nsubjs)*S.sigmaslope + S.muslope;
        init0(i) = S;
    end
    
    % sampling
    tic
    fprintf(   'Running JAGS ...\n'   );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'anal_PTM9.txt'), ...
        init0, ...
        'doparallel', doparallel, ...
        'nchains', nchains, ...
        'nburnin', nburnin, ...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'Na','Af','slope','muNa','muAf','muslope','sigmaNa','sigmaAf','sigmaslope'}, ...
        'savejagsoutput', 1, ...
        'verbosity', 1, ... 
        'cleanup', 0, ...
        'workingdir', workingdir);
    toc

    save(savename,'stats','samples');
    
end

function mat = expmat(varargin)
	mat = [];
	for i = length(varargin):-1:1
		conds = length(varargin{i});
		mat = repmat(mat, conds, 1);
		trials = max([conds, size(mat, 1)]);
		mat = [reshape(repmat(varargin{i}, trials / conds, 1), trials, 1), mat]; %#ok<AGROW>
	end
end