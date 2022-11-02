% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [output_spiketimes] = spike_thinner_temporal(input_spiketimes,n_rep,filt_len)
%SPIKE_THINNER_TEMPORAL Removes spikes from input trains if they are within
%'filt_len' of the previous spike.
if nargin<2 || isempty(n_rep)
    n_rep=1;
end
if nargin<3 || isempty(filt_len)
    filt_len=0.05; % 50 ms
end

n_trial=length(input_spiketimes);
output_spiketimes=cell(n_trial,n_rep);
for trial_ind=1:n_trial
    n_spike=length(input_spiketimes{trial_ind});
    for rep_ind=1:n_rep
        n_thinned=0;
        spike_vec=[];
        run_ISI=filt_len*rand;
        for spike_ind=1:n_spike
            if spike_ind>=2
                run_ISI=run_ISI+(input_spiketimes{trial_ind}(spike_ind)-input_spiketimes{trial_ind}(spike_ind-1));
            else
                run_ISI=run_ISI+(input_spiketimes{trial_ind}(spike_ind));
            end
           if run_ISI>filt_len
               n_thinned=n_thinned+1;
               spike_vec(n_thinned)=input_spiketimes{trial_ind}(spike_ind);
               run_ISI=0;
           end
        end
        output_spiketimes{trial_ind,rep_ind}=spike_vec;
    end                
end

