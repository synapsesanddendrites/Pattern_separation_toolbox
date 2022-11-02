% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [input_spiketimes] = temporally_correlated(max_time,n_trial,rate,regularity)
%TEMPORALLY_CORRELATED Generates a set of spike times with correlations within trains.

input_spiketimes=cell(n_trial,1);
for trial_ind=1:n_trial
    these_spiketimes=zeros(1);
    running_time=exprnd(1/rate);
    n_spike=0;
    while running_time<max_time
        n_spike=n_spike+1;
        these_spiketimes(n_spike)=running_time;
        running_time=running_time+gamrnd(regularity,1/(rate*regularity));
    end
    input_spiketimes{trial_ind}=these_spiketimes;
end
end