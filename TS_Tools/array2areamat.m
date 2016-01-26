function [F, cindx, mask] = array2areamat(array_in, mask, area_ids, mval)
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
%           mask    [r x c]     Binary mask for removing undesired pixels 
%                               from 
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
if nargin < 4, mval = NaN; end

% Length of the time-series
nts  = size(array_in, 1);
% Copy the mask nts-times such that the dimensions of array_in and mask
% have the same length
mask = repmat(mask, 1, 1, nts);
mask = permute(mask, [3 1 2]);

% Set all grid cells of mask where array_in contains missing values to zero
if isnan(mval)
    mask(isnan(array_in))  = 0;
else
    mask(array_in == mval) = 0;
end

% Compute the sum of mask along the first (time) dimension
mask              = sum(mask, 1);
% Remove the first dimension
mask              = squeeze(mask);
% Set all elements with non-continuous time-series to zero
mask(mask < nts)  = 0;
% Set the remaining elements to 1
mask(mask == nts) = 1;
% Transform the map into a vector
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
    
 
 
 

    












