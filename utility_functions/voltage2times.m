% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [spiketimes]=voltage2times(voltage_trace)
% Converts a voltage trace to a vector of spike times using the findpeaks function of the signal processing toolbox.
warning('off','all');
ts=voltage_trace(1,:);
vs=voltage_trace(2,:);
dts=diff(ts);
nums=find(dts<=0);
for sooth=1:length(nums)
    ts((nums(sooth)+1):end)=(ts((nums(sooth)+1):end)-2*dts(nums(sooth))+realmin); % Correct time resolution (for example for a 'Neuron' simulation with variable timesteps)
end
[~,spiketimes] = findpeaks(vs,ts,'MinPeakHeight',-10,'MinPeakDistance',2);
warning('on','all');
end
