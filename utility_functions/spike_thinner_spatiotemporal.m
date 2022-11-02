% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [output_spiketimes] = spike_thinner_spatiotemporal(input_spiketimes,n_rep,filt_len)
%SPIKE_THINNER_SPATIOTEMPORAL Removes spikes from input trains if they are within
%'filt_len' of any other spike in the ensemble.
if nargin<2 || isempty(n_rep)
    n_rep=1;
end
if nargin<3 || isempty(filt_len)
    filt_len=0.005; % 5 ms
end

n_trial=length(input_spiketimes);
output_spiketimes=cell(n_trial,n_rep);
% Do spike filtering
for rep_ind=1:n_rep
    f_inputs=input_spiketimes;
    n_spike_vec=zeros(n_trial,1);
    for trial_ind=1:n_trial
        n_spike=length(f_inputs{trial_ind});
        n_spike_vec(trial_ind)=n_spike;
    end
    n_total=sum(n_spike_vec(:)); % Total number of spikes
    n_thinned=0;
    unsorted_spikes=cell(0);
    while n_total>0
        n_cum_total=cumsum(n_spike_vec(:)); % Cumulative total number of spikes
        % Establish first spike
        st_spike_ind=ceil(n_total*rand);
        found=0;
        trial_ind=1;
        while found==0
            running_count=n_cum_total(trial_ind);
            if st_spike_ind<=running_count
                st_trial_ind=trial_ind;
                if trial_ind>=2
                    st_spike_time=f_inputs{trial_ind}(st_spike_ind-n_cum_total(trial_ind-1));
                else
                    st_spike_time=f_inputs{trial_ind}(st_spike_ind);
                end
                n_thinned=n_thinned+1;
                unsorted_spikes{n_thinned}=[st_trial_ind,st_spike_time];
                found=1;
            else
                trial_ind=trial_ind+1;
            end
        end
        % Now remove all other close spikes
        for trial_ind=1:n_trial
            for spike_ind=n_spike_vec(trial_ind):-1:1
                if abs(f_inputs{trial_ind}(spike_ind)-st_spike_time)<(filt_len/2)
                    f_inputs{trial_ind}(spike_ind)=[];
                end
            end
        end
        for trial_ind=1:n_trial
            n_spike=length(f_inputs{trial_ind});
            n_spike_vec(trial_ind)=n_spike;
        end
        n_total=sum(n_spike_vec(:)); % Total number of spikes
    end
    % Recapitulate existing spikes
    n_filt_tot=length(unsorted_spikes);
    for trial_ind=1:n_trial
        n_thinned=0;
        spike_vec=[];
        for spike_ind=1:n_filt_tot
            if unsorted_spikes{spike_ind}(1)==trial_ind
                n_thinned=n_thinned+1;
                spike_vec(n_thinned)=unsorted_spikes{spike_ind}(2);
            end
        end
        output_spiketimes{trial_ind,rep_ind}=sort(spike_vec);
    end
end

end

