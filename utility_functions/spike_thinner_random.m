% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [output_spiketimes] = spike_thinner_random(input_spiketimes,n_rep,p_del)
%SPIKE_THINNER_RANDOM Randomly removes spikes from input trains
if nargin<2 || isempty(n_rep)
    n_rep=1;
end
if nargin<3 || isempty(p_del)
    p_del=0.5;
end

n_trial=length(input_spiketimes);
output_spiketimes=cell(n_trial,n_rep);
for trial_ind=1:n_trial
    n_spike=length(input_spiketimes{trial_ind});
    for rep_ind=1:n_rep
        n_thinned=0;
        spike_vec=[];
        for spike_ind=1:n_spike
            p_tst=rand();
            if p_tst>p_del
                n_thinned=n_thinned+1;
                spike_vec(n_thinned)=input_spiketimes{trial_ind}(spike_ind);
            end
        end
        output_spiketimes{trial_ind,rep_ind}=spike_vec;
    end                
end

