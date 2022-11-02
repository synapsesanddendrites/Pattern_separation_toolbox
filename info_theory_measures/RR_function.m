% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [RR_obj] = RR_function(input,output,code,num_params,options)
%RR_FUNCTION Calculates the redundancy of the input and output ensembles
if nargin<3 || isempty(code)
    code={'spatial','spatial'};
elseif ischar(code)
    code={code,code};
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
if size(input,1)~=size(output,1)
    if strcmp(code{1},'rate_local') || strcmp(code{2},'rate_local') || strcmp(code{1},'temporal') || strcmp(code{2},'temporal') 
        error('Inputs and outputs are encoded by different numbers of units and an invalid local code is specified.')
    end
end

    function [temp_RR] = eval_RR(bsz,input,output,code,wordsize)
        [binned_input,max_time]=bin_inputs(input,code{1},bsz,wordsize);
        [binned_output]=bin_outputs(output,code{2},bsz,max_time,wordsize);
        [redundancies_input]=joint_distribution_RR(binned_input,code{1},wordsize);
        [redundancies_output]=joint_distribution_RR(binned_output,code{2},wordsize);
        temp_RR=min(redundancies_output(redundancies_output>0))-min(redundancies_input(redundancies_input>0));
    end

    function [RR_in,RR_out] = eval_RR_final(bsz,input,output,code,wordsize)
        [binned_input,max_time]=bin_inputs(input,code{1},bsz,wordsize);
        [binned_output]=bin_outputs(output,code{2},bsz,max_time,wordsize);
        [redundancies_input]=joint_distribution_RR(binned_input,code{1},wordsize);
        [redundancies_output]=joint_distribution_RR(binned_output,code{2},wordsize);
        RR_out=min(redundancies_output(redundancies_output>0));
        RR_in=min(redundancies_input(redundancies_input>0));
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
    RRs_by_code=zeros(length(possible_codes));
    RR_max=0;
    for code_ind_in=1:length(possible_codes)
        for code_ind_out=1:length(possible_codes)
            anon_eval_RR=@(bsz)eval_RR(bsz,input,output,{possible_codes{code_ind_in},possible_codes{code_ind_out}},num_params.wordsize);
            x = patternsearch(anon_eval_RR,num_params.binsize,[],[],[],[],1,1000,[],go_options);
            RR_here = eval_RR(x,input,output,{possible_codes{code_ind_in},possible_codes{code_ind_out}},num_params.wordsize);
            RRs_by_code(code_ind_in,code_ind_out)=abs(RR_here);
            RR_obj.all_binsizes(code_ind_in,code_ind_out)=x;
            if abs(RR_here)>RR_max
                best_code={possible_codes{code_ind_in},possible_codes{code_ind_out}};
                best_binsize=x;
                RR_max=abs(RR_here);
                [RR_in,RR_out] = eval_RR_final(x,input,output,{possible_codes{code_ind_in},possible_codes{code_ind_out}},num_params.wordsize);
            end
        end
    end
    RR_obj.best_code=best_code;
    RR_obj.best_binsize=best_binsize;
    RR_obj.RR=RR_max;
    RR_obj.RRin=RR_in;
    RR_obj.RRout=RR_out;
    RR_obj.RRbycode=RRs_by_code;
elseif contains(options,'-max_bins') % Maximise over bin size
    anon_eval_RR=@(bsz)eval_RR(bsz,input,output,code,num_params.wordsize);
    x = patternsearch(anon_eval_RR,num_params.binsize,[],[],[],[],1,1000,[],go_options);
    RR_here = eval_RR(x,input,output,code,num_params.wordsize);
    [RR_in,RR_out] = eval_RR_final(x,input,output,code,num_params.wordsize);
    RR_obj.best_code=code;
    RR_obj.best_binsize=x;
    RR_obj.RR=abs(RR_here);
    RR_obj.RRin=RR_in;
    RR_obj.RRout=RR_out;
else
    RR_here = eval_RR(num_params.binsize,input,output,code,num_params.wordsize);
    [RR_in,RR_out] = eval_RR_final(num_params.binsize,input,output,code,num_params.wordsize);
    RR_obj.best_code=code;
    RR_obj.best_binsize=num_params.binsize;
    RR_obj.RR=abs(RR_here);
    RR_obj.RRin=RR_in;
    RR_obj.RRout=RR_out;
end
end

