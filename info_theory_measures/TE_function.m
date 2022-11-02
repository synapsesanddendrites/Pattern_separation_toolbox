% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [TE_obj] = TE_function(input,output,code,num_params,options)
%TE_FUNCTION Calculates the transfer entropy from the input to the output ensembles
if nargin<3 || isempty(code)
    code={'spatial','spatial'};
elseif ischar(code)
    code={code,code};
end
if size(input,1)~=size(output,1)
    if strcmp(code{1},'rate_local') || strcmp(code{2},'rate_local') || strcmp(code{1},'temporal') || strcmp(code{2},'temporal') 
        error('Inputs and outputs are encoded by different numbers of units and an invalid local code is specified.')
    end
end
if nargin<4 || isempty(num_params)
  num_params=struct;
end
if ~isfield(num_params,'binsize')
    num_params.binsize=50;
end
if ~isfield(num_params,'wordsize')
    num_params.wordsize=5;
end
if ~isfield(num_params,'delay')
    num_params.delay=50;
end
if ~isfield(num_params,'n_trial')
    num_params.n_trial=size(output,1);
end
if nargin<5 || isempty(options)
    options='-max_bins -par';
end


function[temp_TE] = eval_TE(bsz,delay,input,output,code,wordsize)
        [binned_input,max_time]=bin_inputs(input,code{1},bsz,wordsize);
        [binned_output]=bin_outputs(output,code{2},bsz,max_time,wordsize);
        [binned_delayed_output]=bin_outputs_delay(output,code{2},bsz,delay,max_time,wordsize);
        [joint_distribution,delayed_distribution]=joint_distribution_TE(binned_input,code{1},binned_output,binned_delayed_output,code{2},wordsize);
        tTE = compute_TE(joint_distribution,delayed_distribution);
        temp_TE = -tTE;
    end


if contains(options,'-par')
    go_options = optimoptions(@patternsearch,'UseParallel',true,'Display','none');
else
    go_options = optimoptions(@patternsearch,'Display','none');
end


if contains(options,'-max_code') % Maximise over spike codes, bin size, and delay
    if size(input,1)==size(output,1)
        possible_codes={'spatial','rate_ensemble','rate_local','temporal'};
    else
        possible_codes={'spatial','rate_ensemble'};
    end
    TE_max=0;
    for pre_code_ind=1:length(possible_codes)
        for post_code_ind=1:length(possible_codes)
            anon_eval_TE=@(X)eval_TE(X(1),X(2),input,output,{possible_codes{pre_code_ind},possible_codes{post_code_ind}},num_params.wordsize);
            x = patternsearch(anon_eval_TE,[num_params.binsize,num_params.delay],[],[],[],[],[1,1],[1000,1000],[],go_options);
            TE_here = eval_TE(x(1),x(2),input,output,{possible_codes{pre_code_ind},possible_codes{post_code_ind}},num_params.wordsize);

            TE_obj.all_codes{pre_code_ind,post_code_ind}={possible_codes{pre_code_ind},possible_codes{post_code_ind}};
            TE_obj.all_binsizes(pre_code_ind,post_code_ind)=x(1);
            TE_obj.all_delays(pre_code_ind,post_code_ind)=x(2);
            TE_obj.all_TEs(pre_code_ind,post_code_ind)=abs(TE_here);
            if abs(TE_here)>TE_max
                best_code={possible_codes{pre_code_ind},possible_codes{post_code_ind}};
                best_binsize=x(1);
                best_delay=x(2);
                TE_max=abs(TE_here);
            end
        end
    end
    TE_obj.best_code=best_code;
    TE_obj.best_binsize=best_binsize;
    TE_obj.best_delay=best_delay;
    TE_obj.TE=TE_max;
elseif contains(options,'-max_bins') % Maximise over bin size and delay
    anon_eval_TE=@(X)eval_TE(X(1),X(2),input,output,code,num_params.wordsize);
    x = patternsearch(anon_eval_TE,[num_params.binsize,num_params.delay],[],[],[],[],[1,1],[1000,1000],[],go_options);
    TE_here = eval_TE(x(1),x(2),input,output,code,num_params.wordsize);
    TE_obj.best_code=code;
    TE_obj.best_binsize=x(1);
    TE_obj.best_delay=x(2);
    TE_obj.MI=abs(TE_here);
else
    TE_here = eval_TE(num_params.binsize,num_params.delay,input,output,code,num_params.wordsize);
    TE_obj.best_code=code;
    TE_obj.best_binsize=num_params.binsize;
    TE_obj.best_delay=num_params.delay;
    TE_obj.MI=abs(TE_here);
end
end