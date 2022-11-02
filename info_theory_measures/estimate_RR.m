% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [estRR_obj] = estimate_RR(input,output,options)
%ESTIMATE_RR Estimates the redundancy of the input and output ensembles.

if nargin<3 || isempty(options)
    options='-par';
end

n_trial_in=size(input,1);
n_trial_out=size(output,1);
max_time=0;
for ward=1:n_trial_in
    fspike=input{ward}(end);
    if fspike>max_time
        max_time=fspike;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input redundancy
all_in_RR=cell(n_trial_in,1);
RRs_in=zeros(n_trial_in,1);
for ward=1:n_trial_in
    RR_in=estimate_MI(input(ward,:),input(setdiff(1:n_trial_in,ward),:),options);
    all_in_RR{ward}=RR_in;
    RRs_in(ward)=RR_in.MI;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output redundancy
all_out_RR=cell(n_trial_out,1);
RRs_out=zeros(n_trial_out,1);
for ward=1:n_trial_out
    RR_out=estimate_MI(output(ward,:),output(setdiff(1:n_trial_out,ward),:),options,max_time);
    all_out_RR{ward}=RR_out;
    RRs_out(ward)=RR_out.MI;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare everything

estRR_obj=struct;
estRR_obj.RR=min(RRs_in)-min(RRs_out);
estRR_obj.RRin=min(RRs_in);
estRR_obj.RRout=min(RRs_out);
[~,amin_in]=min(RRs_in);
[~,amin_out]=min(RRs_out);
estRR_obj.RRin_bounds=all_in_RR{amin_in,1}.MI_bounds;
estRR_obj.RRout_bounds=all_out_RR{amin_out,1}.MI_bounds;
estRR_obj.RRin_vector=all_in_RR{amin_in,1}.MI_vector;
estRR_obj.RRout_vector=all_out_RR{amin_out,1,1}.MI_vector;

end
