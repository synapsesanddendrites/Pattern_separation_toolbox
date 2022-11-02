% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [input_spiketimes] = spatially_correlated(max_time,n_trial,rate,corr)
%SPATIALLY_CORRELATED Generates a set of spike times with correlations within trains.

master_rate=rate/corr;
input_spiketimes=cell(n_trial,1);
master_spiketimes=zeros(1);
running_time=exprnd(1/master_rate);
n_spike=0;
while running_time<max_time
    n_spike=n_spike+1;
    master_spiketimes(n_spike)=running_time;
    running_time=running_time+exprnd(1/master_rate);
end

for trial_ind=1:n_trial
    these_spiketimes=zeros(1);
    spikes_here=0;
    for spike_ind=1:n_spike
       tst=rand(1);
       if tst<corr
           spikes_here=spikes_here+1;
           these_spiketimes(spikes_here)=master_spiketimes(spike_ind)+0.005*randn(1);
       end
    end    
    input_spiketimes{trial_ind}=these_spiketimes;
end
end