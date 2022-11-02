% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [joint_distribution,input_distribution,output_distribution]=joint_distribution_MI(binned_inputs,input_code,binned_outputs,output_code,wordsize)
    n_trial_in=size(binned_inputs,1);
    n_trial_out=size(binned_outputs,1);
    n_rep=size(binned_outputs,2);
    if ~strcmp(output_code,'temporal')
        t_len=size(binned_outputs,3);
    else
        t_len=size(binned_outputs,3)/wordsize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Map inputs to codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(input_code,'spatial')
        all_words=double(dec2bin(0:2^n_trial_in-1)' - '0'); % All possible input spiking patterns
        [~,i_inputs]=ismember(binned_inputs',all_words','rows'); % Match input patterns
        edge_vec=linspace(0.5,2^n_trial_in+0.5,2^n_trial_in+1);
        input_distribution=histcounts(i_inputs,edge_vec)';
        cd_inputs=reshape(i_inputs',[1,1,t_len]);
    elseif strcmp(input_code,'rate_ensemble')
        binned_inputs=sum(binned_inputs); % Take ensemble rates
        max_input_rate=max(binned_inputs(:));
        edge_vec=linspace(-0.5,max_input_rate+0.5,max_input_rate+2);
        input_distribution=histcounts(binned_inputs(:),edge_vec)';
        i_inputs=binned_inputs+1;
        cd_inputs=reshape(i_inputs,[1,1,t_len]);
    elseif strcmp(input_code,'rate_local')
        max_input_rate=max(binned_inputs(:));
        edge_vec=linspace(-0.5,max_input_rate+0.5,max_input_rate+2);
        i_inputs=zeros(n_trial_in,t_len);
        for trial_ind=1:n_trial_in
            i_inputs(trial_ind,:)=binned_inputs(trial_ind)+1;
        end
        input_distribution=histcounts(binned_inputs(:),edge_vec)';
        cd_inputs=reshape(i_inputs,[n_trial_in,1,t_len]);
    elseif strcmp(input_code,'temporal')
        all_words=double(dec2bin(0:2^wordsize-1)' - '0'); % All possible input spiking patterns
        i_inputs=zeros(n_trial_in,t_len);
        for trial_ind=1:n_trial_in
            ordered_inputs=reshape(binned_inputs(trial_ind,:),[wordsize,t_len]);
            [~,j_inputs]=ismember(ordered_inputs',all_words','rows'); % Match input patterns
            i_inputs(trial_ind,:)=j_inputs;
        end
        edge_vec=linspace(0.5,2^wordsize+0.5,2^wordsize+1);
        input_distribution=histcounts(i_inputs(:),edge_vec)';
        cd_inputs=reshape(i_inputs,[n_trial_in,1,t_len]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% Map outputs to codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(output_code,'spatial')
        all_words=double(dec2bin(0:2^n_trial_out-1)' - '0'); % All possible input spiking patterns
        i_outputs=zeros(n_rep,t_len);
        for rep_ind=1:n_rep
            [~,j_outputs]=ismember(squeeze(binned_outputs(:,rep_ind,:))',all_words','rows'); % Match output patterns
            i_outputs(rep_ind,:)=j_outputs;
        end
        edge_vec=linspace(0.5,2^n_trial_out+0.5,2^n_trial_out+1);
        output_distribution=histcounts(i_outputs(:),edge_vec)';
        cd_outputs=reshape(i_outputs,[1,n_rep,t_len]);
    elseif strcmp(output_code,'rate_ensemble')
        binned_outputs=sum(binned_outputs,1); % Take ensemble rates
        max_output_rate=max(binned_outputs(:));
        i_outputs=zeros(n_rep,t_len);
        for rep_ind=1:n_rep
            i_outputs(rep_ind,:)=binned_outputs(1,rep_ind,:)+1;
        end
        edge_vec=linspace(-0.5,max_output_rate+0.5,max_output_rate+2);
        output_distribution=histcounts(binned_outputs(:),edge_vec)';
        cd_outputs=reshape(i_outputs,[1,n_rep,t_len]);
    elseif strcmp(output_code,'rate_local')
        max_output_rate=max(binned_outputs(:));
        i_outputs=zeros(n_trial_out,n_rep,t_len);
        for trial_ind=1:n_trial_out
            for rep_ind=1:n_rep
                i_outputs(trial_ind,rep_ind,:)=binned_outputs(trial_ind,rep_ind,:)+1;
            end
        end
        edge_vec=linspace(-0.5,max_output_rate+0.5,max_output_rate+2);
        output_distribution=histcounts(binned_outputs(:),edge_vec)';
        cd_outputs=reshape(i_outputs,[n_trial_out,n_rep,t_len]);
    elseif strcmp(output_code,'temporal')
        all_words=double(dec2bin(0:2^wordsize-1)' - '0'); % All possible input spiking patterns
        i_outputs=zeros(n_trial_out,n_rep,t_len);
        for trial_ind=1:n_trial_out
            for rep_ind=1:n_rep
                ordered_outputs=reshape(binned_outputs(trial_ind,rep_ind,:),[wordsize,t_len]);
                [~,j_outputs]=ismember(ordered_outputs',all_words','rows'); % Match input patterns
                i_outputs(trial_ind,rep_ind,:)=j_outputs;
            end
        end
        edge_vec=linspace(0.5,2^wordsize+0.5,2^wordsize+1);
        output_distribution=histcounts(i_outputs(:),edge_vec)';
        cd_outputs=reshape(i_outputs,[n_trial_out,n_rep,t_len]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% Count joint occurrences %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ex_in=size(cd_inputs,1);
    ex_out=size(cd_outputs,1);   

    joint_distribution=zeros(length(input_distribution),length(output_distribution));
    if ex_in==ex_out % Matched inputs and outputs
        for trial_ind=1:ex_in
            for rep_ind=1:n_rep
               joint_distribution=joint_distribution+accumarray([squeeze(cd_inputs(trial_ind,1,:)) , squeeze(cd_outputs(trial_ind,rep_ind,:))],1,[length(input_distribution),length(output_distribution)]);
            end
        end
    elseif ex_in==1 % Input ensemble code
        for trial_ind=1:ex_out
            for rep_ind=1:n_rep
                joint_distribution=joint_distribution+accumarray([squeeze(cd_inputs(1,1,:)) , squeeze(cd_outputs(trial_ind,rep_ind,:))],1,[length(input_distribution),length(output_distribution)]);
            end
        end
    elseif ex_out==1 % Output ensemble code
        for trial_ind=1:ex_in
            for rep_ind=1:n_rep
                joint_distribution=joint_distribution+accumarray([squeeze(cd_inputs(trial_ind,1,:)) , squeeze(cd_outputs(1,rep_ind,:))],1,[length(input_distribution),length(output_distribution)]);
            end
        end
    else % Other mismatch in sizes; all-to-all comparison
        for trial_ind_in=1:ex_in
            for trial_ind_out=1:ex_out
                for rep_ind=1:n_rep
                    joint_distribution=joint_distribution+accumarray([squeeze(cd_inputs(trial_ind_in,1,:)) , squeeze(cd_outputs(trial_ind_out,rep_ind,:))],1,[length(input_distribution),length(output_distribution)]);
                end
            end
        end
    end

end
