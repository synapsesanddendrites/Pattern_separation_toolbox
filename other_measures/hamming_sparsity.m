% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [hamming_obj] = hamming_sparsity(input,output,num_params)
%HAMMING_SPARSITY  Records Hamming distance (corrected for sparsity) between spike trains

if nargin<3 || isempty(num_params)
   num_params=struct;
end
if ~isfield(num_params,'binsize')
    num_params.binsize=10;
end
if ~isfield(num_params,'n_trial')
    num_params.n_trial=size(output,1);
end    

[binned_inputs,max_time]=bin_inputs(input,'spatial',num_params.binsize);
[binned_outputs]=bin_outputs(output,'spatial',num_params.binsize,max_time); 

n_trial_in=size(binned_inputs,1);
n_trial_out=size(binned_outputs,1);
n_rep=size(binned_outputs,2);

HDis=zeros(n_trial_in);
HDi_mn=0;
for ii=1:(n_trial_in-1)
    for jj=(ii+1):n_trial_in
        HDi_vec=abs(binned_inputs(ii,:)-binned_inputs(jj,:));
        HDis(ii,jj)=sum(HDi_vec(:));
        HDi_mn=HDi_mn+HDis(ii,jj)/(n_trial_in*(n_trial_in-1)/2);
    end
end
fr_i=sum(binned_inputs(:))/numel(binned_inputs);

HDos=zeros(n_rep,n_trial_out,n_trial_out);
fr_o=zeros(n_rep,1);
HDo_mn=zeros(n_rep,1);

for ward=1:n_rep
    for ii=1:(n_trial_out-1)
        for jj=(ii+1):n_trial_out
            HDo_vec=abs(binned_outputs(ii,ward,:)-binned_outputs(jj,ward,:));
            HDos(ward,ii,jj)=sum(HDo_vec(:));
            HDo_mn(ward)=HDo_mn(ward)+HDos(ward,ii,jj)/(n_trial_out*(n_trial_out-1)/2);
        end
    end
    fr_o(ward)=sum(binned_outputs(:,2,:),'all')/numel(binned_outputs(:,2,:));
end

f1_vec=zeros(n_rep,1);
for ward=1:n_rep
   f1_vec(ward)=HDi_mn/(2*(1-fr_i)*n_trial_in)-HDo_mn(ward)/(2*(1-fr_o(ward))*n_trial_out);
end

hamming_obj.in_Hd=HDi_mn/(2*(1-fr_i));
hamming_obj.out_Hd=HDo_mn./(2*(1-fr_o));
hamming_obj.hammings=f1_vec;
hamming_obj.hamming_est=mean(f1_vec);

end



