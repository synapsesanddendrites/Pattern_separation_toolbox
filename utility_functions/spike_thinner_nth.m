% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [output_spiketimes] = spike_thinner_nth(input_spiketimes,n_rep,n)
%SPIKE_THINNER_NTH Removes spikes from input trains if they are within 'n'
%of the previous allowed spike.
if nargin<2 || isempty(n_rep)
    n_rep=1;
end
if nargin<3 || isempty(n)
    n=2;
end

n_trial=length(input_spiketimes);
output_spiketimes=cell(n_trial,n_rep);
for trial_ind=1:n_trial
    n_spike=length(input_spiketimes{trial_ind});
    for rep_ind=1:n_rep
        n_thinned=0;
        spike_vec=[];
        run_count=floor(n*rand);
        for spike_ind=1:n_spike
           run_count=run_count+1;
           if run_count==n
               n_thinned=n_thinned+1;
               spike_vec(n_thinned)=input_spiketimes{trial_ind}(spike_ind);
               run_count=0;
            end
        end
        output_spiketimes{trial_ind,rep_ind}=spike_vec;
    end                
end

