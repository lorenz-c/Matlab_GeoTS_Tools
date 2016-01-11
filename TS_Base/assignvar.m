function varargout = assignvar(filename, cmd_otpt)
%--------------------------------------------------------------------------
% This function allows to load a .mat-file with the name "filename"
% directly into the variables defined by varargout.
%--------------------------------------------------------------------------
% INPUT:
% - filename  Filename of the file which should be loaded
% - cmd_otpt  True: Write comments on the command-line, False: No comments
%--------------------------------------------------------------------------
% OUTPUT:
% - varargout Cell-array of output variables
%--------------------------------------------------------------------------
% Example:  Test_Var = assignvar('test_var.mat', true)
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         October 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% Turn command-line output off (by default)
if nargin < 2, cmd_otpt = false; end

% Load a file -> struct. with name "tmp"
tmp     = load(filename);

% Read the contents of the file (i.e., the different variables in "tmp")
fnmes   = fieldnames(tmp);

% Check the max. number output variables
mx_vars = max(nargout, length(fnmes));

if nargout > length(fnmes)
    % Show an error if there are more output arguments than variables.
    error('Too many output arguments!')
end

% Write some comments on the command line
if cmd_otpt == true
    fprintf('Reading from file %s .... \n', filename);
    fprintf('Assigning the following variables: \n');
    for i = 1:length(fnmes), fprintf('%i. %s \n', i, fnmes{i}); end
end

% Assign the variables in the input file to one (or multiple) output
% variable(s)
for i = 1:mx_vars
    varargout{i} = getfield(tmp, fnmes{i}); 
end
