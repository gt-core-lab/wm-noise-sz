function PTM_out = getDataPTM(ptm_csv_path, trialSubIDs)
%
% Reads a PTM CSV and returns a table with columns:
%   SubID (string)   -- full IDs like '2001AB' to match Trials.SubID
%   sigma_add        -- PTM additive internal noise (Na-fix)
%   phi              -- PTM external-noise filtering (Af-fix)
%   k1var, k2var, k4var (optional) -- empirical 1/kappa for set sizes 1,2,4
%
% INPUTS
%   ptm_csv_path : path to PTM csv file
%   trialSubIDs  : string/cellstr array of SubIDs used in Trials (e.g., {'2001AB','3005XZ',...})
%
% The PTM CSV must contain:
%   studyID  (4-digit ID matching first 4 chars of SubID)
%   Na-fix   (or Na_fix)  --> additive internal noise
%   Af-fix   (or Af_fix)  --> external noise filtering
% It may also contain (optional, recommended):
%   k1var, k2var, k4var   --> empirical 1/kappa per set size
%
% OUTPUT
%   PTM_out : table aligned to the unique SubIDs in trialSubIDs
%             with columns: SubID, sigma_add, phi, [k1var,k2var,k4var if present]

% --- read PTM file
PTMraw = readtable(ptm_csv_path);

% normalize variable names (handle hyphens vs underscores)
vn = PTMraw.Properties.VariableNames;
vn_norm = strrep(vn, '-', '_');
PTMraw.Properties.VariableNames = vn_norm;

% required columns
req = {'studyID','Na_fix','Af_fix'};
has = ismember(req, PTMraw.Properties.VariableNames);
if ~all(has)
    missing = req(~has);
    error('PTM file is missing required column(s): %s', strjoin(missing, ', '));
end

% optional empirical columns
optCols = {'k1var','k2var','k4var'};
hasOpt  = ismember(optCols, PTMraw.Properties.VariableNames);

% cast studyID to string with zero padding if numeric
if isnumeric(PTMraw.studyID)
    % assume 4-digit numeric IDs
    PTMraw.studyID = arrayfun(@(x) sprintf('%04d', x), PTMraw.studyID, 'UniformOutput', false);
end
PTMraw.studyID = string(PTMraw.studyID);

% pull PTM parameters with standard names
PTMbase = table();
PTMbase.studyID   = PTMraw.studyID;
PTMbase.sigma_add = PTMraw.Na_fix;   % additive internal noise
PTMbase.phi       = PTMraw.Af_fix;   % external-noise filtering

% attach optional empirical 1/kappa columns if present
for i = 1:numel(optCols)
    if hasOpt(i)
        PTMbase.(optCols{i}) = PTMraw.(optCols{i});
    else
        PTMbase.(optCols{i}) = nan(height(PTMbase),1);
    end
end

% ensure trialSubIDs is string
if iscellstr(trialSubIDs), trialSubIDs = string(trialSubIDs); end
if ~isstring(trialSubIDs)
    trialSubIDs = string(trialSubIDs);
end
uSub = unique(trialSubIDs(:), 'stable');

% first 4 characters (studyID)
if exist('extractBefore','file') == 2
    studyFromSub = extractBefore(uSub, 5);   % '2000KN' -> '2000'
else
    % fallback
    studyFromSub = string(regexprep(cellstr(uSub), '^(.{4}).*$', '$1'));
end

% join: map each SubID (with letters) to its PTM row via studyID
[tf, loc] = ismember(studyFromSub, PTMbase.studyID);
if ~all(tf)
    miss = uSub(~tf);
    warning('PTM not found for %d SubID(s): %s', sum(~tf), strjoin(cellstr(miss), ', '));
end

% build output table aligned to *every* SubID present in Trials
PTM_out = table();
PTM_out.SubID     = uSub;
PTM_out.sigma_add = nan(size(uSub));
PTM_out.phi       = nan(size(uSub));
PTM_out.k1var     = nan(size(uSub));
PTM_out.k2var     = nan(size(uSub));
PTM_out.k4var     = nan(size(uSub));

ok = find(tf);
PTM_out.sigma_add(ok) = PTMbase.sigma_add(loc(ok));
PTM_out.phi(ok)       = PTMbase.phi(loc(ok));
PTM_out.k1var(ok)     = PTMbase.k1var(loc(ok));
PTM_out.k2var(ok)     = PTMbase.k2var(loc(ok));
PTM_out.k4var(ok)     = PTMbase.k4var(loc(ok));

% optional sanity messages
n_ok  = numel(ok);
n_all = numel(uSub);
fprintf('PTM mapping: %d/%d SubIDs matched via studyID (first 4 digits).\n', n_ok, n_all);
end
