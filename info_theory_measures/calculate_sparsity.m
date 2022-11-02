% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [sparsity_obj ] = calculate_sparsity(input,output)
% CALCULATE_SPARSITY Calculates relative sparsity between input and output
% spike ensembles
n_spike_in=0;
for sooth=1:size(input,1)
    n_spike_in=n_spike_in+length(input{sooth});
end

n_spike_out=0;
for sooth=1:size(output,1)
    for rep=1:size(output,2)
        n_spike_out=n_spike_out+length(output{sooth,rep});
    end
end
sparsity=(n_spike_in-n_spike_out/size(output,2))/n_spike_in;

sparsity_obj.sparsity=sparsity;
sparsity_obj.n_spike_in=n_spike_in;
sparsity_obj.n_spike_out=n_spike_out;
end

