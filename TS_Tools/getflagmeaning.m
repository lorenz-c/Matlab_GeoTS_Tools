function flags_out = getflagmeaning(flag_values, flag_meanings, base, pos)

if nargin < 4, pos = 1:length(flag_values); end
if nargin < 3, base = length(flag_meanings); end

% Split the flag-meanings string into a cell array
if ~iscell(flag_meanings)
    flag_meanings = strsplit(flag_meanings, ' ');
end

flag_meanings = fliplr(flag_meanings);
if length(flag_meanings) < base
   flag_meanings = [cell(1, base - length(flag_meanings)) flag_meanings];
end

% Convert the decimal flag values to binary strings
bin_strings = dec2bin(flag_values(pos), base);
bin_strings = cellstr(bin_strings);

% Search for '1's in the bin_strings
true_pos = strfind(bin_strings, '1');

for i = 1:length(pos)
    flags_out{i, 1} = pos(i);
    flags_out{i, 1} = strjoin(flag_meanings(true_pos{i}));
end

