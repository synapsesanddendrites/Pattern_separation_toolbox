% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [MI_obj] = MI_function(input,output,code,num_params,options)
%MI_FUNCTION Calculates the mutual information between the input and output ensembles

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
if ~isfield(num_params,'n_trial')
    num_params.n_trial=size(output,1);
end
if nargin<5 || isempty(options)
    options='-max_bins -par';
end

    function[temp_MI] = eval_MI(bsz,input,output,code,wordsize)
        [binned_input,max_time]=bin_inputs(input,code{1},bsz,wordsize);
        [binned_output]=bin_outputs(output,code{2},bsz,max_time,wordsize);
        [joint_distribution,input_distribution,output_distribution]=joint_distribution_MI(binned_input,code{1},binned_output,code{2},wordsize);
        tMI = compute_MI(joint_distribution,input_distribution,output_distribution);
        temp_MI = -tMI;
    end

if contains(options,'-par')
    go_options = optimoptions(@patternsearch,'UseParallel',true,'Display','none');
else
    go_options = optimoptions(@patternsearch,'Display','none');
end

if contains(options,'-max_code') % Maximise over spike codes and bin size
    if size(input,1)==size(output,1)
        possible_codes={'spatial','rate_ensemble','rate_local','temporal'};
    else
        possible_codes={'spatial','rate_ensemble'};
    end
    MI_max=0;
    for pre_code_ind=1:length(possible_codes)
        for post_code_ind=1:length(possible_codes)
            anon_eval_MI=@(bsz)eval_MI(bsz,input,output,{possible_codes{pre_code_ind},possible_codes{post_code_ind}},num_params.wordsize);
            x = patternsearch(anon_eval_MI,num_params.binsize,[],[],[],[],1,1000,[],go_options);
            MI_here = eval_MI(x,input,output,{possible_codes{pre_code_ind},possible_codes{post_code_ind}},num_params.wordsize);

            MI_obj.all_codes{pre_code_ind,post_code_ind}={possible_codes{pre_code_ind},possible_codes{post_code_ind}};
            MI_obj.all_binsizes(pre_code_ind,post_code_ind)=x;
            MI_obj.all_MIs(pre_code_ind,post_code_ind)=abs(MI_here);
            if abs(MI_here)>MI_max
                best_code={possible_codes{pre_code_ind},possible_codes{post_code_ind}};
                best_binsize=x;
                MI_max=abs(MI_here);
            end
        end
    end
    MI_obj.best_code=best_code;
    MI_obj.best_binsize=best_binsize;
    MI_obj.MI=MI_max;
elseif contains(options,'-max_bins') % Maximise over bin size
    anon_eval_MI=@(bsz)eval_MI(bsz,input,output,code,num_params.wordsize);
    x = patternsearch(anon_eval_MI,num_params.binsize,[],[],[],[],1,1000,[],go_options);
    MI_here = eval_MI(x,input,output,code,num_params.wordsize);
    MI_obj.best_code=code;
    MI_obj.best_binsize=x;
    MI_obj.MI=abs(MI_here);
else
    MI_here = eval_MI(num_params.binsize,input,output,code,num_params.wordsize);
    MI_obj.best_code=code;
    MI_obj.best_binsize=num_params.binsize;
    MI_obj.MI=abs(MI_here);
end
end



