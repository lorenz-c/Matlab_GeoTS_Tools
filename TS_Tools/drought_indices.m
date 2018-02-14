function out = drought_indices(inpt, scale, type)

nts = size(inpt, 1);
out = NaN(size(inpt));

if strcmp(type, 'SPI_daily')
    
elseif strcmp(type, 'SPI_param_monthly')
    % Fit a gamma distribution to the input data
    %pd = fitdist(inpt, 'Gamma');
elseif strcmp(type, 'SPI_emp_monthly')
    
    accum = zeros(nts - scale + 1, size(inpt, 2), size(inpt, 3));
    
    % Create the accumulated precipitation
    for i = 1:scale
        accum = accum + inpt(i : nts - scale + i, :, :);
    end

    p_trafo = zeros(size(accum));
    
    
    % Now, let's fit an empirical distribution
    % TBA: Right now, it is assumed that the first month is January!
    for i = 1:12
        % Get the precipitation of the current month
        mnthly_p = accum(i:12:end, :, :);
        dst      = zeros(size(mnthly_p));
        % Loop over each included timestep
        for j = 1:size(mnthly_p, 1)
            % Substract the actual month from the remaining months
            delta = bsxfun(@minus, mnthly_p, mnthly_p(j, :, :));
            % ...and count the number of elements which are <= 0 --> eCDF
            dst(j, :, :) = sum(delta <= 0);
        end
        % Compute the marginal probability of precipitation using the
        % empirical Gringorten plotting position
        dst = (dst - 0.44)/(size(mnthly_p, 1) + 0.12);
        p_trafo(i:12:end, :, :) = dst;
    end
    % Compute the SPIs using the inverse standard normal distribution
    out(scale : end, :, :) = norminv(p_trafo);
  
end