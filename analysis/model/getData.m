function [errs, Ns] = getData(paths, subid)

% paths.remote = fullfile('/Users', 'wpark78', 'Dropbox', 'WoonJuPark', 'Projects', 'active');
% paths.project = fullfile(paths.remote, 'vandy-wm', 'wm');
% paths.data = fullfile(paths.project, '@Data');
% subid = '2000KN';

fname = fullfile(paths.data, [subid, '_VWMp.csv']);

T = readtable(fname);

% response error (degrees)
if ismember('RespError', T.Properties.VariableNames)
    err_deg = T.RespError;
else
    error('Column "RespError" not found in CSV.');
end

% set size
if ismember('SetSize', T.Properties.VariableNames)
    Ns = T.SetSize;
elseif ismember('N', T.Properties.VariableNames)
    Ns = T.N;
end

% convert to radians
errs = deg2rad(wrapTo180(err_deg));

