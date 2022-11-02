% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [orth_obj] = orthogonalisation(input,output,num_params)
% ORTHOGONALISATION Records normalised dot product of the spike trains

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

n_rep=size(binned_outputs,2);
t_len=size(binned_outputs,3);

[in_NDP,in_SF]=dot_product(binned_inputs);

all_NDPs=zeros(n_rep,1);
all_SFs=zeros(n_rep,1);
    
for tune=1:n_rep
    i_vec=zeros(num_params.n_trial,t_len);
    for bear=1:num_params.n_trial
        i_vec(bear,:)=binned_outputs(bear,tune,:);
    end
    [o_NDP,o_SF] = dot_product(i_vec);

    all_NDPs(tune)=in_NDP-o_NDP;
    all_SFs(tune)=in_SF-o_SF;
end

orth_obj.in_NDP=in_NDP;
orth_obj.in_SF=in_SF;

orth_obj.NDPs=all_NDPs(~isnan(all_NDPs));
orth_obj.SFs=all_SFs(~isnan(all_SFs));
orth_obj.NDP_est=mean(orth_obj.NDPs);
orth_obj.SF_est=mean(orth_obj.SFs);

end

function [NDP_av,SF_av] = dot_product(input)

n_trial=size(input,1);
NDP_av=0;
SF_av=0;

for ii=2:n_trial
    for jj=1:(ii-1)
        X=input(ii,:);
        Y=input(jj,:)';

        Xbar=norm(X);
        Ybar=norm(Y);
        if Xbar>0 && Ybar>0
            NDP=(X*Y)/(Xbar*Ybar);
            if Xbar>Ybar
                SF=Ybar/Xbar;
            else
                SF=Xbar/Ybar;
            end
        NDP_av=NDP_av+NDP/((n_trial)*(n_trial-1)/2);
        SF_av=SF_av+SF/((n_trial)*(n_trial-1)/2);    
        else
            NDP_av=NaN;
            SF_av=NaN;
        end
        
    end
    
end

end