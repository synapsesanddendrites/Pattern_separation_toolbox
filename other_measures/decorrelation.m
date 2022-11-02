% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [decorr_obj] = decorrelation(input,output,num_params)
%DECORRELATION Records normalised decorrelation of the spike trains

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

input_corrs=corrcoef(binned_inputs');

n_trial_in=size(binned_inputs,1);
n_trial_out=size(binned_outputs,1);
n_rep=size(binned_outputs,2);
t_len=size(binned_outputs,3);

all_corrs=zeros(n_rep,n_trial_out,n_trial_out);
for tune=1:n_rep
    i_vec=zeros(n_trial_out,t_len);
    for bear=1:n_trial_out
        i_vec(bear,:)=binned_outputs(bear,tune,:);
    end
    icorr=corrcoef(i_vec');
    all_corrs(tune,:,:)=icorr;
end

decorr_vec=zeros(n_rep,1);
try % Standard system
    for ward=1:n_rep
        decorr_here=0;
        for ii=2:n_trial_in
            for jj=1:(ii-1)
                i_decorr=input_corrs(jj,ii)-all_corrs(ward,jj,ii);
                decorr_here=decorr_here+i_decorr/((n_trial_in)*(n_trial_in-1)/2);
            end
        end
        decorr_vec(ward)=decorr_here;
    end
catch % Expansion system with different numbers of units
    input_corr=0;
    for ii=2:n_trial_in
        for jj=1:(ii-1)
            input_corr=input_corr+input_corrs(jj,ii)/((n_trial_in)*(n_trial_in-1)/2);
        end
    end

    for ward=1:n_rep
        outcorr_here=0;
        for ii=2:n_trial_out
            for jj=1:(ii-1)
                outcorr_here=outcorr_here+all_corrs(ward,jj,ii)/((n_trial_out)*(n_trial_out-1)/2);
            end
        end
        decorr_vec(ward)=input_corr-outcorr_here;
    end
end

decorr_obj.out_corr_mat=all_corrs;
decorr_obj.input_corrs=input_corrs;
decorr_obj.decorrs=decorr_vec;
decorr_obj.decorr_est=mean(decorr_vec);

end

