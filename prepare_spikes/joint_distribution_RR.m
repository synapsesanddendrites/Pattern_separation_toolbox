function [redundancies]=joint_distribution_RR(binned_spikes,code,wordsize)
    if ismatrix(binned_spikes)
        sz=size(binned_spikes);
        binned_spikes=reshape(binned_spikes,[sz(1),1,sz(2)]);
    end
    n_trial=size(binned_spikes,1);
    n_rep=size(binned_spikes,2);
    if ~strcmp(code,'temporal')
        t_len=size(binned_spikes,3);
    else
        t_len=size(binned_spikes,3)/wordsize;
    end

    redundancies=zeros(n_trial,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Map spikes to redundancies %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(code,'spatial') % 'Leave-one-out' ensemble code
        all_words=double(dec2bin(0:2^(n_trial-1)-1)' - '0'); % All possible input spiking patterns
        edge_vec=linspace(0.5,2^(n_trial-1)+0.5,2^(n_trial-1)+1);
        L1=2;
        L2=2^(n_trial-1);
        for out_ind=1:n_trial
            joint_distribution=zeros(L1,L2);
            out_distribution=zeros(L1,1);
            rem_distribution=zeros(L2,1);
            for rep_ind=1:n_rep            
                left_out_vector=squeeze(binned_spikes(out_ind,rep_ind,:)); % Take out vector
                rem_matrix=squeeze(binned_spikes(:,rep_ind,:));
                rem_matrix(out_ind,:)=[]; % Delete left-out vector
                
                %%%%%%%%%%% Get codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [~,i_rems]=ismember(rem_matrix',all_words','rows'); % Match input patterns  
                cd_out=reshape(left_out_vector+1,[1,1,t_len]);
                cd_rem=reshape(i_rems',[1,1,t_len]);
                
                out_distribution=out_distribution+histcounts(left_out_vector,[-0.5,0.5,1.5])';
                rem_distribution=rem_distribution+histcounts(i_rems,edge_vec)';
                joint_distribution=joint_distribution+accumarray([squeeze(cd_out) , squeeze(cd_rem)],1,[L1,L2]);
            end
            redundancies(out_ind) = compute_RR(joint_distribution,out_distribution,rem_distribution);
        end
    elseif strcmp(code,'rate_ensemble') % 'Leave-one-out' ensemble code
        for out_ind=1:n_trial
            left_out_matrix=binned_spikes(out_ind,:,:); % Take out vector
            rem_tensor=binned_spikes;
            rem_tensor(out_ind,:,:)=[];
            rem_ensemble=sum(rem_tensor,1);

            max_rem_rate=max(rem_ensemble(:));
            max_out_rate=max(left_out_matrix(:));
            edge_vec_rem=linspace(0.5,max_rem_rate+1.5,max_rem_rate+2);
            edge_vec_out=linspace(0.5,max_out_rate+1.5,max_out_rate+2);

            joint_distribution=zeros(max_out_rate+1,max_rem_rate+1);
            out_distribution=zeros(max_out_rate+1,1);
            rem_distribution=zeros(max_rem_rate+1,1);
            for rep_ind=1:n_rep            
                left_out_vector=squeeze(left_out_matrix(:,rep_ind,:)); % Take out vector
                rem_matrix=squeeze(rem_tensor(:,rep_ind,:));
                
                %%%%%%%%%%% Get codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                i_rems=sum(rem_matrix,1)+1;
                i_outs=left_out_vector+1;
                cd_out=reshape(i_outs,[1,1,t_len]);
                cd_rem=reshape(i_rems,[1,1,t_len]);

                out_distribution=out_distribution+histcounts(i_outs,edge_vec_out)';
                rem_distribution=rem_distribution+histcounts(i_rems,edge_vec_rem)';
                joint_distribution=joint_distribution+accumarray([squeeze(cd_out) , squeeze(cd_rem)],1,[max_out_rate+1,max_rem_rate+1]);
            end
            redundancies(out_ind) = compute_RR(joint_distribution,out_distribution,rem_distribution);
        end    
    elseif strcmp(code,'rate_local') % Pairwise redundancy code
        redundancies=zeros(n_trial);
        for out_ind=1:(n_trial-1)
            left_out_matrix=binned_spikes(out_ind,:,:); % Take out vector
            max_out_rate=max(left_out_matrix(:));
            edge_vec_out=linspace(0.5,max_out_rate+1.5,max_out_rate+2);
            out_distribution=zeros(max_out_rate+1,1);
            for comp_ind=(out_ind+1):n_trial
                comp_matrix=binned_spikes(comp_ind,:,:); % Comparison vector                
                max_comp_rate=max(comp_matrix(:));
                comp_distribution=zeros(max_comp_rate+1,1);

                edge_vec_comp=linspace(0.5,max_comp_rate+1.5,max_comp_rate+2);
                joint_distribution=zeros(max_out_rate+1,max_comp_rate+1);            
                for rep_ind=1:n_rep            
                    left_out_vector=squeeze(left_out_matrix(1,rep_ind,:)); % Take out vector
                    comp_vector=squeeze(comp_matrix(1,rep_ind,:));
                    %%%%%%%%%%% Get codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    i_outs=left_out_vector+1;
                    i_comps=comp_vector+1;

                    cd_out=reshape(i_outs,[1,1,t_len]);
                    cd_comp=reshape(i_comps,[1,1,t_len]);

                    out_distribution=out_distribution+histcounts(i_outs,edge_vec_out)';
                    comp_distribution=comp_distribution+histcounts(i_comps,edge_vec_comp)';
                    joint_distribution=joint_distribution+accumarray([squeeze(cd_out) , squeeze(cd_comp)],1,[max_out_rate+1,max_comp_rate+1]);
                end
                redundancies(out_ind,comp_ind) = compute_RR(joint_distribution,out_distribution,comp_distribution);
            end
        end
    elseif strcmp(code,'temporal')
        redundancies=zeros(n_trial);
        all_words=double(dec2bin(0:2^wordsize-1)' - '0'); % All possible input spiking patterns
        edge_vec=linspace(0.5,2^wordsize+0.5,2^wordsize+1);
        for out_ind=1:(n_trial-1)
            left_out_matrix=binned_spikes(out_ind,:,:); % Take out vector
            comp_distribution=zeros(2^wordsize,1);
            out_distribution=zeros(2^wordsize,1);
            joint_distribution=zeros(2^wordsize);
            for comp_ind=(out_ind+1):n_trial
                comp_matrix=binned_spikes(comp_ind,:,:); % Comparison vector                       
                for rep_ind=1:n_rep            
                    left_out_vector=squeeze(left_out_matrix(1,rep_ind,:)); % Take out vector
                    comp_vector=squeeze(comp_matrix(1,rep_ind,:));
                
                    %%%%%%%%%%% Get codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ordered_leftouts=reshape(left_out_vector,[wordsize,t_len]);
                    ordered_comps=reshape(comp_vector,[wordsize,t_len]);

                    [~,i_outs]=ismember(ordered_leftouts',all_words','rows'); % Match input patterns
                    [~,i_comps]=ismember(ordered_comps',all_words','rows'); % Match input patterns
                    
                    cd_out=reshape(i_outs,[1,1,t_len]);
                    cd_comp=reshape(i_comps,[1,1,t_len]);

                    out_distribution=out_distribution+histcounts(i_outs,edge_vec)';
                    comp_distribution=comp_distribution+histcounts(i_comps,edge_vec)';
                    joint_distribution=joint_distribution+accumarray([squeeze(cd_out) , squeeze(cd_comp)],1,[2^n_trial,2^n_trial]);
                end
                redundancies(out_ind,comp_ind) = compute_RR(joint_distribution,out_distribution,comp_distribution);
            end
        end    
    end
end
