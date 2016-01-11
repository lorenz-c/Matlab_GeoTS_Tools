function map = signal_map(inpt, indx_vec, id_map, miss)
%--------------------------------------------------------------------------
% The function assigns a value (an element in the inpt-vector) to an area
% of a map refferenced by the id_map. 
%--------------------------------------------------------------------------
% INPUT:
% - inpt      Vector which contains the values to be assigned to a map
% - indx_vec  Vector which contains the ids of the appropriate areas#
% - id_map    Map where each "continuous" region has an ID which 
%             corresponds to the indx_vec
% - miss      Defines the missing elements in the output map (i.e. the
%             regions which are not considered by either the id_map or the
%             indx_vec; default: NaN).
%--------------------------------------------------------------------------
% OUTPUT:
% - map       Map with the same dimension as the id_map, where each pixel 
%             within a continuous area has the same value, as defined by 
%             inpt
%--------------------------------------------------------------------------
% EXAMPLE:
% >> map = signal_map([1 2 3], [111 222 333], id_map);
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         August 2012
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

if nargin < 4, miss = NaN; end

if length(inpt) ~= length(indx_vec)
    error('Inpt and indx_vec must have the same length!')
end

A = zeros(size(id_map)).*miss;

for i = 1:length(indx_vec)
    A(id_map == indx_vec(i)) = inpt(i);
end



