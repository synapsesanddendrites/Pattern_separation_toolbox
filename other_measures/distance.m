% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [dist_obj] = distance(input,output,num_params)
%DISTANCE Records Wasserstein distance between spike trains

if nargin<3 || isempty(num_params)
   num_params=struct;
end
if ~isfield(num_params,'n_trial_in')
    num_params.n_trial_in=size(input,1);
end    
if ~isfield(num_params,'n_trial_out')
    num_params.n_trial_out=size(output,1);
end    
n_trial_in=num_params.n_trial_in;
n_trial_out=num_params.n_trial_out;
n_rep=size(output,2);

in_WS=zeros(n_trial_in);
in_mean_dist=0;
for ii=1:(n_trial_in-1)
    for jj=(ii+1):n_trial_in
        d_ii_jj=spike_wasserstein(input{ii},input{jj});
        in_mean_dist=in_mean_dist+d_ii_jj*2/(n_trial_in*(n_trial_in-1));
        in_WS(ii,jj)=d_ii_jj;
    end
end

out_WS=zeros(n_rep,n_trial_out,n_trial_out);
out_mean_dists=zeros(n_rep,1);
for tune=1:n_rep
    out_mean_dist=0;
    for ii=1:(n_trial_out-1)
        for jj=(ii+1):n_trial_out
            d_ii_jj=spike_wasserstein(output{ii,tune},output{jj,tune});
            out_mean_dist=out_mean_dist+d_ii_jj*2/(n_trial_out*(n_trial_out-1));
            out_WS(tune,ii,jj)=d_ii_jj;
        end
    end
    out_mean_dists(tune)=out_mean_dist;
end

dist_obj.in_WS_mat=in_WS;
dist_obj.out_WS_mat=out_WS;
dist_obj.in_WS_mu=in_mean_dist;
dist_obj.out_WS_mu_vec=out_mean_dists;
dist_obj.out_WS_mu=mean(out_mean_dists(:));
dist_obj.dist_est=1-min(in_mean_dist/mean(out_mean_dists(:)),mean(out_mean_dists(:))/in_mean_dist);

end

function [distance] = spike_wasserstein(st1,st2)
n_sp1=length(st1);
n_sp2=length(st2);

cdf_1=linspace(0,1,n_sp1+1);
cdf_2=linspace(0,1,n_sp2+1);

st_vec1_full=zeros(2*n_sp1,2);
for ii=1:n_sp1
    st_vec1_full(2*(ii-1)+1,1)=st1(ii)-realmin;
    st_vec1_full(2*(ii-1)+2,1)=st1(ii);
    st_vec1_full(2*(ii-1)+1,2)=cdf_1(ii);
    st_vec1_full(2*(ii-1)+2,2)=cdf_1(ii+1);
end

st_vec2_full=zeros(2*n_sp2,2);
for ii=1:n_sp2
    st_vec2_full(2*(ii-1)+1,1)=st2(ii)-realmin;
    st_vec2_full(2*(ii-1)+2,1)=st2(ii);
    st_vec2_full(2*(ii-1)+1,2)=cdf_2(ii);
    st_vec2_full(2*(ii-1)+2,2)=cdf_2(ii+1);
end

all_spikes=sort([st1,st2]);
n_all=length(all_spikes);
st_vec_diff=zeros(2*(n_all),2);
for ii=1:n_all
    st_vec_diff(2*(ii-1)+1,1)=all_spikes(ii)-2*realmin;
    st_vec_diff(2*(ii-1)+2,1)=all_spikes(ii);
    st_vec_diff(2*(ii-1)+1,2)=abs(step_function_interp(st_vec1_full,all_spikes(ii)-2*realmin)-step_function_interp(st_vec2_full,all_spikes(ii)-2*realmin));
    st_vec_diff(2*(ii-1)+2,2)=abs(step_function_interp(st_vec1_full,all_spikes(ii))-step_function_interp(st_vec2_full,all_spikes(ii)));
end
distance=trapz(st_vec_diff(:,1),st_vec_diff(:,2));
end

function [val] = step_function_interp(svec,q)
t=svec(:,1);
t_ind=find(t>=q,1);
if length(t_ind)==1
    val=svec(t_ind,2);
else
    val=1;
end
end