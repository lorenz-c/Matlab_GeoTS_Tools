function R = nanrms(fld, dim)
%--------------------------------------------------------------------------
% Simple function which computes the root mean square from a given input
% field by taking into account NaNs.
%--------------------------------------------------------------------------
% INPUT:
% - fld       Input vector, matrix, or array
% - dim       Dimension from which the rms should be derived (default: 1)
%--------------------------------------------------------------------------
% OUTPUT:
% - rm        RMS of fld
%--------------------------------------------------------------------------
% EXAMPLE:
% >> rm = nanrms(rand(100, 10, 10), 3);
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         October 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% Create a mask with ones
mask = ones(size(fld));
% Remove all missing elements in fld from mask
mask(isnan(fld)) = 0;
% Sum over fld
nr_obs = sum(mask, dim);

R = sqrt((1./nr_obs).*nansum(fld.^2, dim));







