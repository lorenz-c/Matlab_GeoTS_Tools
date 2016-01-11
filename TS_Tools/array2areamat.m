function [F, cindx] = array2areamat(array_in, mask, mval)
% The function rearranges the elements in an cell-array (e.g. an array of
% maps) to a big matrix, which has a number of rows equal to the number of
% time-steps in flds and a number of columns equal to the pixels of a
% single field. The function allows the consideration of a mask to reduce
% the size of a matrix, if e.g. large parts of the input fields contain
% missing values or if further computations are needed for a specific area
% only. 
%--------------------------------------------------------------------------
% Input:    flds    {m x 1}     Cell array (or single matrix) which contains 
%                               the input fields.                
%           mask    [r x c]     Binary mask for removing undesired pixels from 
%                               the flds-cells
%           mval    [1 x 1]     If mval is set, the function also searches for 
%                               missing values in flds and removes these 
%                               elements
%                               Default: -9999
%           arr     1 or 2      arr = 1: Longitude ordering
%                               arr = 2: Latitude ordering
%                               Detault: arr = 1;             
% Output:   F       [m x n]     Matrix which has a number of rows equal to
%                               the number of timesteps in flds (i.e. m)
%                               and a number of columns equal to the number
%                               of ones (and the number of missing values)
%                               in mask
%           cindx   [n x 1]     Column-vector which contains the positions
%                               of the valid elements in flds
%--------------------------------------------------------------------------
% Author: Christof Lorenz
% Date:   November 2015
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------
if nargin < 3
    mval = NaN;
end

nts  = size(array_in, 1);
mask = repmat(mask, 1, 1, nts);
mask = permute(mask, [3 1 2]);

if isnan(mval)
    mask(isnan(array_in))  = 0;
else
    mask(array_in == mval) = 0;
end

mask             = sum(mask, 1);
mask             = squeeze(mask);
mask(mask < nts) = 0;
mask(mask > 0)   = 1;

mask_vec = mask(:);

% Find the positions of the elements which are ~= 0 
cindx = find(mask_vec == 1);

% Create the matrix F...
F = zeros(nts, sum(sum(mask)));
    
for i = 1:nts
	tmp = squeeze(array_in(i, :, :));
    tmp = tmp(:);
    
    F(i,:) = tmp(cindx);  
end
    
 
 
 

    












