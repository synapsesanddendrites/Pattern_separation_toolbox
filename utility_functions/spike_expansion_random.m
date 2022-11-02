function [shuffled_output_spiketimes] = spike_expansion_random(input_spiketimes,n_exp,n_rep,p_del)
%SPIKE_EXPANSION_RANDOM Expands input spike trains into random numbers of output units
if nargin<2 || isempty(n_exp)
    n_exp=5;
end
if nargin<3 || isempty(n_rep)
    n_rep=8;
end
if nargin<4 || isempty(p_del)
    p_del=0.8;
end

n_trial=length(input_spiketimes);
output_spiketimes=cell(0,n_rep);
output_index=1;
for trial_ind=1:n_trial
    n_spike=length(input_spiketimes{trial_ind});
    real_n_exp=ceil(n_exp*rand(1)); % Number of units to expand to
    for exp_ind=1:real_n_exp
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
            output_spiketimes{output_index,rep_ind}=spike_vec;
        end
        output_index=output_index+1;
    end
end
shuffled_output_spiketimes=cell(output_index-1,n_rep);

for rep_ind=1:n_rep
    shuffled_order=randperm(output_index-1);
    for ward=1:(output_index-1)
        shuffled_output_spiketimes{shuffled_order(ward),rep_ind}=output_spiketimes{ward,rep_ind};
    end
end
end
