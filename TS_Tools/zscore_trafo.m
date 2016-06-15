function z = zscore_trafo(inpt, vars)


% First, compute monthly averages from the data
inpt_mn  = ts_average(inpt, 'monthly_lt');

% Compute the mean monthly standard devition
inpt_std = ts_average(inpt, 'monthly_lt', 'nanstd');

inpt_mn  = anncycle2fullts(inpt_mn,  inpt.Data.time);
inpt_std = anncycle2fullts(inpt_std, inpt.Data.time);


z = inpt;

dta = inpt.Data.(vars);
dta = (dta - inpt_mn.Data.(vars))./inpt_std.Data.(vars);


z.Data.(vars) = dta;