function [estMI_obj] = estimate_MI(input,output,options,max_time)
%ESTIMATE_RR Estimates the mutual information between the input and output ensembles.

if nargin<3 || isempty(options)
    options='-par';
end

if contains(options,'-par') % Establish number of estimates to use in each step
    p=gcp();
    num_chunk=p.NumWorkers;
else
    num_chunk=8;
end

if nargin<4 || isempty(max_time)
    % Determine last spike
    n_trial=size(input,1);
    max_time=0;
    for ward=1:n_trial
        fspike=input{ward}(end);
        if fspike>max_time
            max_time=fspike;
        end
    end
end
n_chunk_vec=zeros(0,3);
max_chunks=ceil(max_time*10); % Start with 100ms chunks
min_chunks=floor(max_time*4); % 250 ms chunks;

fit_not_found=true;
while fit_not_found
    try % Only works with Statistics toolbox
        n_chunk=randsample(min_chunks:max_chunks,num_chunk,false);
    catch % Otherwise use a simpler method
        A=min_chunks:max_chunks;
        n_chunk=A(randi(numel(A),num_chunk,1));
    end
    ests_here=cell(length(n_chunk),1);
    if contains(options,'-par')
        parfor ward=1:length(n_chunk)
            [Imax,hmax] = KL_MI(n_chunk(ward),max_time,input,output);
            ests_here{ward}=[n_chunk(ward) , Imax , hmax];
        end
    else
        for ward=1:length(n_chunk)
            [Imax,hmax] = KL_MI(n_chunk(ward),max_time,input,output);
            ests_here{ward}=[n_chunk(ward) , Imax , hmax];
        end
    end
    nS=size(n_chunk_vec,1);
    for ward=1:length(n_chunk)
        n_chunk_vec(nS+ward,:)=ests_here{ward};
    end
    i=n_chunk_vec(:,3)<=(0.5*n_chunk_vec(:,1)); % Values with a consistent smoothing parameter
    if nnz(i)>=5
        Xi=n_chunk_vec(i,1);
        Yi=n_chunk_vec(i,2);
        [xData, yData] = prepareCurveData( Xi, Yi );
        ft = fittype( 'poly1' );
        [fitresult, gof] = fit( xData, yData, ft );
        if gof.rsquare>0.95
            estMI_obj=struct;
            estMI_obj.MI=fitresult.p2;
            estMI_obj.MI_vector=n_chunk_vec;
            ci=confint(fitresult);
            estMI_obj.MI_bounds=[ci(1,2) , ci(2,2)];
            fit_not_found=false;
        end
    else
        min_chunks=max_chunks;
        max_chunks=ceil(1.25*max_chunks); % Increase resolution
    end
end

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
t1=st_vec1_full(:,1);
t2=st_vec2_full(:,1);

all_spikes=sort([st1,st2]);
n_all=length(all_spikes);
st_vec_diff=zeros(2*(n_all),2);
for ii=1:n_all
    st_vec_diff(2*(ii-1)+1,1)=all_spikes(ii)-2*realmin;
    st_vec_diff(2*(ii-1)+2,1)=all_spikes(ii);
    st_vec_diff(2*(ii-1)+1,2)=abs(step_function_interp(t1,st_vec1_full,all_spikes(ii)-2*realmin)-step_function_interp(t2,st_vec2_full,all_spikes(ii)-2*realmin));
    st_vec_diff(2*(ii-1)+2,2)=abs(step_function_interp(t1,st_vec1_full,all_spikes(ii))-step_function_interp(t2,st_vec2_full,all_spikes(ii)));
end
distance=trapz(st_vec_diff(:,1),st_vec_diff(:,2));
end

function [val] = step_function_interp(t,svec,q)
t_ind=find(t>=q,1);
if length(t_ind)==1
    val=svec(t_ind,2);
else
    val=1;
end
end

function[Imax,hmax] = KL_MI(n_chunk,max_time,input,output)

% Measure all pairwise distances
grid_times=linspace(0,max_time,n_chunk+1);
grid_period=grid_times(2);
n_trial_in=size(input,1);
n_trial_out=size(output,1);
n_rep=size(output,2);

% Divide spike trains into chunks
chunked_trains_in=cell(n_trial_in,n_chunk);
chunked_trains_out=cell(n_trial_out,n_rep,n_chunk);
for trial_ind=1:n_trial_in
    sps_in=input{trial_ind};
    for chunk_ind=1:n_chunk
        shfted_sps_in=sps_in-grid_times(chunk_ind);
        chunked_trains_in{trial_ind,chunk_ind}=shfted_sps_in(shfted_sps_in>=0 & shfted_sps_in<grid_period);
    end
end
for trial_ind=1:n_trial_out
    for rep_ind=1:n_rep
        sps_out=output{trial_ind,rep_ind};
        for chunk_ind=1:n_chunk
            shfted_sps_out=sps_out-grid_times(chunk_ind);
            chunked_trains_out{trial_ind,rep_ind,chunk_ind}=shfted_sps_out(shfted_sps_out>=0 & shfted_sps_out<grid_period);
        end
    end
end

% Compute pairwise distances
p_dists_in=zeros(n_trial_in,n_chunk,n_chunk);
p_dists_out=zeros(n_trial_out,n_rep,n_chunk,n_chunk);
for trial_ind=1:n_trial_in
    for ii_chunk_ind=1:n_chunk
        for jj_chunk_ind=ii_chunk_ind:n_chunk
            p_dists_in(trial_ind,ii_chunk_ind,jj_chunk_ind)=spike_wasserstein(chunked_trains_in{trial_ind,ii_chunk_ind},chunked_trains_in{trial_ind,jj_chunk_ind});
        end
    end
end
for trial_ind=1:n_trial_out
    for rep_ind=1:n_rep
        for ii_chunk_ind=1:n_chunk
            for jj_chunk_ind=ii_chunk_ind:n_chunk
                p_dists_out(trial_ind,rep_ind,ii_chunk_ind,jj_chunk_ind)=spike_wasserstein(chunked_trains_out{trial_ind,rep_ind,ii_chunk_ind},chunked_trains_out{trial_ind,rep_ind,jj_chunk_ind});
            end
        end
    end
end

% Collate over ensembles
ens_dists_in=zeros(n_chunk);
ens_dists_out=zeros(n_rep,n_chunk,n_chunk);
for ii_chunk_ind=1:n_chunk
    for jj_chunk_ind=ii_chunk_ind:n_chunk
        ens_here_in=0;
        for trial_ind=1:n_trial_in
            ens_here_in=ens_here_in+p_dists_in(trial_ind,ii_chunk_ind,jj_chunk_ind)^2;
        end
        ens_dists_in(ii_chunk_ind,jj_chunk_ind)=sqrt(ens_here_in);
        for rep_ind=1:n_rep
            ens_here_out=0;
            for trial_ind=1:n_trial_out
                ens_here_out=ens_here_out+p_dists_out(trial_ind,rep_ind,ii_chunk_ind,jj_chunk_ind)^2;
            end
            ens_dists_out(rep_ind,ii_chunk_ind,jj_chunk_ind)=sqrt(ens_here_out);
        end
    end
end
ens_dists_in=ens_dists_in+ens_dists_in';
ens_dists_out=ens_dists_out+permute(ens_dists_out,[1 3 2]);

% Compute MI using KL method
full_grid=zeros(n_chunk-2,2);
full_grid(:,1)=2:(n_chunk-1); % Keep track of everything already
h_space=2:(n_chunk-2);
I_space=zeros(length(h_space),1);
for ward=1:length(h_space)
    if full_grid(h_space(ward),2)==0
        I_space(ward)=h_fun(h_space(ward),ens_dists_in,ens_dists_out,n_chunk,n_rep);
    else
        I_space(ward)=full_grid(h_space(ward),2);
    end
end
[~,minH]=min(I_space);
Imax=-I_space(minH);
hmax=h_space(minH);

end

function[fact_out]=log_factorial(n)

if n<=10
    fact_out=log(factorial(n));
else
    fact_out=log(sqrt(2*pi*n))+n*log(n/exp(1))+log(1+1/(12*n)+1/(288*n^2)-139/(51840*n^3)-571/(2488320*n^4)); % Stirling's approximation
end

end

function h_score = h_fun(h_in,ens_dists_in,ens_dists_out,n_chunk,n_rep)

h_score_here=zeros(n_rep,1);
for i=1:n_chunk
    C_in=sort(ens_dists_in(:,i));
    C_in_val=C_in(h_in);
    for rept_ind=1:n_rep
        C_out=sort(ens_dists_out(rept_ind,:,i));
        C_out_val=C_out(h_in);

        hash_here=nnz((ens_dists_in(:,i)<=C_in_val).*(ens_dists_out(rept_ind,:,i)<=C_out_val));
        h_score_here(rept_ind)=h_score_here(rept_ind)+(1/n_chunk)*log2(n_chunk*hash_here/(h_in^2));
    end
end
I0_h=0;
for i=1:h_in
    if (h_in-i)<=(n_chunk-h_in)
        log_prob=2*log_factorial(h_in-1)+2*log_factorial(n_chunk-h_in)-log_factorial(n_chunk-1)-log_factorial(i-1)-2*log_factorial(h_in-i)-log_factorial(n_chunk-2*h_in+i);
        I0_h=I0_h+exp(log_prob)*log2(n_chunk*i/(h_in^2));
    end
end
h_score=-(mean(h_score_here(:))-I0_h);
end



