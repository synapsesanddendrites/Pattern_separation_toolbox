% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [binned_spikes,max_time]=bin_inputs(input,input_code,binsize_ms,wordsize)
if nargin<4 || wordsize
    wordsize=5;
end
n_trial=size(input,1);

% Determine last spike
max_time=0;
for ward=1:n_trial
    fspike=input{ward}(end);
    if fspike>max_time
        max_time=fspike;
    end
end

% Establish grid
if ~strcmp(input_code,'temporal')
    binsize=binsize_ms/1000;
    grid_len=ceil((max_time+binsize/2)/binsize);
    binned_spikes=zeros(n_trial,grid_len);
else
    binsize=binsize_ms/(1000*wordsize); % Correct for different codes
    grid_len=wordsize*ceil((max_time+wordsize*binsize/2)/(wordsize*binsize));
    binned_spikes=zeros(n_trial,grid_len);
end

% Do binning
if strcmp(input_code,'spatial')
    for ward=1:n_trial
        nspike=length(input{ward});
        for sooth=1:nspike
            binned_spikes(ward,ceil(input{ward}(sooth)/binsize))=1;
        end
    end
elseif strcmp(input_code,'rate_ensemble')
    for ward=1:n_trial
        nspike=length(input{ward});
        for sooth=1:nspike
            binned_spikes(ward,ceil(input{ward}(sooth)/binsize))= binned_spikes(ward,ceil(input{ward}(sooth)/binsize))+1;
        end
    end
elseif strcmp(input_code,'rate_local')
    for ward=1:n_trial
        nspike=length(input{ward});
        for sooth=1:nspike
            binned_spikes(ward,ceil(input{ward}(sooth)/binsize))= binned_spikes(ward,ceil(input{ward}(sooth)/binsize))+1;
        end
    end
elseif strcmp(input_code,'temporal')  
    for ward=1:n_trial
        nspike=length(input{ward});
        for sooth=1:nspike
            binned_spikes(ward,ceil(input{ward}(sooth)/binsize))=1;
        end
    end
end
end