% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [binned_spikes]=bin_outputs_delay(output,output_code,binsize_ms,delay,max_time,wordsize)
n_trial=size(output,1);
n_rep=size(output,2);
if nargin<4 || isempty(max_time)
    max_time=0;
    for ward=1:n_trial
        fspike=output{ward,n_rep}(end);
        if fspike>max_time
            max_time=fspike;
        end
    end
end
if nargin<5 || wordsize
    wordsize=5;
end

% Establish grid
if ~strcmp(output_code,'temporal')
    binsize=binsize_ms/1000;
    grid_len=ceil((max_time+binsize/2)/binsize);
    binned_spikes=zeros(n_trial,n_rep,grid_len);
else
    binsize=binsize_ms/(1000*wordsize); % Correct for different codes
    grid_len=wordsize*ceil((max_time+wordsize*binsize/2)/(wordsize*binsize));
    binned_spikes=zeros(n_trial,n_rep,grid_len);
end

delayed_output=cell(n_trial,n_rep);
for trial_ind=1:n_trial
    for rep_ind=1:n_rep
        delayed_spikes=output{trial_ind,rep_ind}-delay/1000;
        delayed_output{trial_ind,rep_ind}=delayed_spikes(delayed_spikes>0);
    end
end

% Do binning
if strcmp(output_code,'spatial')
    for ward=1:n_trial
        for tune=1:n_rep
            nspike=length(delayed_output{ward,tune});
            for sooth=1:nspike
                binned_spikes(ward,tune,ceil(delayed_output{ward,tune}(sooth)/binsize))=1;
            end
        end
    end
elseif strcmp(output_code,'rate_ensemble')
    for ward=1:n_trial
        for tune=1:n_rep
            nspike=length(delayed_output{ward,tune});
            for sooth=1:nspike
                binned_spikes(ward,tune,ceil(delayed_output{ward,tune}(sooth)/binsize))=binned_spikes(ward,tune,ceil(delayed_output{ward,tune}(sooth)/binsize))+1;
            end
        end
    end
elseif strcmp(output_code,'rate_local')
    for ward=1:n_trial
        for tune=1:n_rep
            nspike=length(delayed_output{ward,tune});
            for sooth=1:nspike
                binned_spikes(ward,tune,ceil(delayed_output{ward,tune}(sooth)/binsize))=binned_spikes(ward,tune,ceil(delayed_output{ward,tune}(sooth)/binsize))+1;
            end
        end
    end
elseif strcmp(output_code,'temporal')
    for ward=1:n_trial
        for tune=1:n_rep
            nspike=length(delayed_output{ward,tune});
            for sooth=1:nspike
                binned_spikes(ward,tune,ceil(delayed_output{ward,tune}(sooth)/binsize))=1;
            end
        end
    end
end
end