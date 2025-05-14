function [] = cluster_PTM_weibull2(condition)

    ngroups = 3;
    ntotalchains = 5;
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
    elseif idgroup == 3
        whichgroup = 'YHC';
    end
    whichchain = emat(condition,3);
    
    % MCMC parameters
    nchains = 1;
    nburnin = 15000;
    nsamples = 2000;
    nthin = 200;
    doparallel = 0;
    
    % set directory and filenames
    loadname = strcat('stacked_PTM_',whichgroup,'.mat');
    savename = strcat('stats_PTM_Weibull2_',whichgroup,'_',num2str(whichchain),'.mat');
    workingdir = strcat('PTM_Weibull2_',whichgroup,'_',num2str(whichchain));
    
    % data
    load(loadname);
    data = data_all;
    nsubjs = length(data);
    
    for i = 1:nsubjs
       c(i,:) = unique(data{i}(:,1)); % noise levels
       tempdata = sortrows(data{i});
       for j = 1:length(c(i,:))
           r(i,:,j) = tempdata(tempdata(:,1)==c(i,j),3); % responses (correct/incorrect)
           x(i,:,j) = log10(tempdata(tempdata(:,1)==c(i,j),2)*1000); % contrasts 
       end
    end
    
    [nsubjs,ntrials,nc] = size(r);
    
    % assign matlab variables to the observed nodes 
    datastruct = struct('x',x,'r',r,'nsubjs',nsubjs,'ntrials',ntrials,'nc',nc);
    
    % assign initial values
    for i = 1:nchains
        for j = 1:nsubjs
            for k = 1:nc
                S.p2(j,k) = rand(1)+2; % 50% threshold
                S.p3(j,k) = rand(1)+3; % slope 
            end
        end
        init0(i) = S;
    end

%     for i = 1:nchains
%         for j = 1:nsubjs
%             S.p3(j) = rand(1)+3; % slope 
%             for k = 1:nc
%                 S.p2(j,k) = rand(1)+2; % 50% threshold
%             end
%         end
%         init0(i) = S;
%     end

%     for i = 1:nchains
%         S.mup3 = 3.5;
%         S.sigmap3 = 0.2;
%         for j = 1:nsubjs
%             S.p3(j) = rand(1)+3; % slope 
%             for k = 1:nc
%                 S.p2(j,k) = rand(1)+2; % 50% threshold
%             end
%         end
%         init0(i) = S;
%     end
    
    % sampling
    tic
    fprintf(   'Running JAGS ...\n'   );
    [samples, stats] = cluster_matjags( ...
        datastruct, ...
        fullfile(pwd, 'anal_PTM_weibull2.txt'), ...
        init0, ...
        'doparallel', doparallel, ...
        'nchains', nchains, ...
        'nburnin', nburnin, ...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'p2','p3'}, ...
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
