%This script generates the mex scripts to accelrate the computations

%check mvgc is on path
onPath = ~isempty(strfind(path, 'mvgc'));

if ~onPath
    disp('mvgc toolbox needs to be on path')
end

% where is adjmat_argminGV_mask
path_of_script = mfilename('fullpath');

% fine last slash
slashes = strfind(path_of_script,'/');
last_slash = slashes(end);
path_of_script = [path_of_script(1:last_slash) '/mex_stuff/'];

% path_of_script = path_of_script
eval(['codegen adjmat_argminGV_mask.m -d ' path_of_script])

% make sure noew mex stuff is on path
addpath(genpath(path_of_script))

