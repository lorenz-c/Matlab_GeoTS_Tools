function [] = single_error_output(errors, filename)


fid = fopen(filename, 'wt');

if ismember('lt_1d', fieldnames(errors))
    err_nms = fieldnames(errors.lt_1d);
    fprintf(fid, '%s \n', 'Long term mean errors from 2D fields')
    for i = 1:length(err_nms)
        fprintf(fid, '%-8s %g \n', err_nms{i}, errors.lt_1d.(err_nms{i}));
    end
    fprintf(fid, ' \n');
end

if ismember('ssnl_1d', fieldnames(errors))
    err_nms = fieldnames(errors.ssnl_1d);
    fprintf(fid, '%s \n', 'Seasonal mean errors from 2D fields')
    for i = 1:length(err_nms)
        fprintf(fid, '%-8s %g %g %g %g \n', err_nms{i}, ...
            errors.ssnl_1d.(err_nms{i})(1), errors.ssnl_1d.(err_nms{i})(2), ...
               errors.ssnl_1d.(err_nms{i})(3), errors.ssnl_1d.(err_nms{i})(4));
    end
    fprintf(fid, ' \n');
end

if ismember('spataverage', fieldnames(errors))
    err_nms = fieldnames(errors.spataverage);
    fprintf(fid, '%s \n', 'Erros from spatial averages')
    for i = 1:length(err_nms)
        fprintf(fid, '%-8s %g \n', err_nms{i}, errors.spataverage.(err_nms{i}));
    end
    fprintf(fid, ' \n');
end

if ismember('spataverage_anom', fieldnames(errors))
    err_nms = fieldnames(errors.spataverage_anom);
    fprintf(fid, '%s \n', 'Erros from spatial averages (anomalies)')
    for i = 1:length(err_nms)
        fprintf(fid, '%-8s %g \n', err_nms{i}, errors.spataverage_anom.(err_nms{i}));
    end
    fprintf(fid, ' \n');
end