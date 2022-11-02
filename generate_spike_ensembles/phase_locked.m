% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [input_spiketimes] = phase_locked(max_time,n_trial,rate,amplitude,num_phase)
%PHASE_LOCKED Generates a set of spike times with sinusoidally varying rates.
dt=0.001;
input_spiketimes=cell(n_trial,1);
for trial_ind=1:n_trial
    these_spiketimes=zeros(1);
    n_spike=0;
    for time_ind=1:ceil(max_time/dt)
        t_here=(time_ind-0.5)*dt;
        rate_here=rate*(1+amplitude*sin(2*pi*num_phase*t_here/max_time));
        p=rand(1);
        if p<(dt*rate_here)
            n_spike=n_spike+1;
            these_spiketimes(n_spike)=t_here;           
        end
    end
    input_spiketimes{trial_ind}=these_spiketimes;
end
end