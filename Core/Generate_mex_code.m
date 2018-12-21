%This script generates the mex scripts to accelrate the computations
% **This requires the matlab MATLAB Coder toolbox**

%check mvgc is on path
onPath = ~isempty(strfind(path, 'mvgc'));

if ~onPath
    disp('mvgc toolbox needs to be on path')
end

% where is adjmat_argminGV_mask
path_of_script = mfilename('fullpath');

% find last slash
slashes = strfind(path_of_script,'/');
last_slash = slashes(end);
path_to_func = path_of_script(1:last_slash);
out_dir_script = [path_to_func '/mex_stuff/'];

%
eval(['codegen adjmat_argminGV_mask.m -d ' out_dir_script])

% make sure new mex stuff is on path
addpath(genpath(path_of_script))

