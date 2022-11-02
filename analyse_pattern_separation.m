% 'analyse_pattern_separation' Master function to analyse pattern separation
% between input and output spiking ensembles using all classical and
% information-theoretic measures.
% 
% (pattern separation toolbox)
%
%[patsep] = analyse_pattern_separation(input,output,options)
% --------------------------------------------------------------
%
% Input
% -----
% - input: Nx1 cell array of input spike times or voltage vectors, where N
% is the number of spike trains in the input ensemble. If inputs are spike
% times, each cell should be a 1xn vector (where n is the number of spikes in 
% that train). If inputs are voltage traces, each cell should be a 2xn
% matrix, where the first row is time, the second row is voltage, and n is the
% number of timesteps recorded.
%
% - output: MxR cell array of output spike times or voltage vectors, where M
% is the number of spike trains in the output ensemble and R is the number of
% repeats. If M is also the number of spike trains in the input, it is assumed
% that each output train corresponds to the single input train with the same index
% (so the first trains are matched and so on); in this case R can also be the number 
% of output units for each input. Matching numbers of input and output trains allow
% the 'local' spike codes to be applied to exact mutual information and transfer
% entropy calculations. If outputs are spike times, each cell should be a 1xn vector
% (where n is the number of spikes in that train). If outputs are voltage traces,
% each cell should be a 2xn matrix, where the first row is time and the second row is voltage.
%
% - options: Name-value pairs specifying desired options.
%          - 'transfer_entropy': boolean specifying if transfer entropy should be
%          calculated instead of mutual information. Overwritten if 'estimate' is
%          set to 'true'. Default: false
%          - 'estimate': boolean specifying if Kozachenko-Leonenko
%          estimators should be used for mutual information and redundancy
%          instead of exact calculation. Overwrites 'transfer_entropy' and 'wordsize'
%          if set to 'true'. Default: false
%          - 'parallel': boolean specifying if some operations should be parallelised
%          to speed up calculations. Default: true
%          - 'classical_binsize': double specifying the binsize (in ms) used to
%          discretise spike times for some classical measures. Default: 10
%          - 'wordsize': integer specifying the length of binary words
%          used in the temporal spiking codes. Overwritten if 'estimate' is
%          set to 'true'. Default: 5
%
% Output
% ------
% - patsep: Structure with fields and subfields quantifying the pattern separation
% between inputs and outputs. See 'Information theoretic measures of pattern 
% separation in the dentate gyrus' (2022) by Bird, Cuntz, and Jedlicka for a
% full description of output fields. Fields and subfields are:
%     'info'   : Information theoretic measures with subfields:
%           'spMI'   : Sparsity-weighted mutual information between inputs and
%           outputs
%           'spTE'   : Sparsity-weighted transfer entropy between inputs and
%           outputs (computed instead of spMI if the optional argument 'transfer_entropy'
%           is set to 'true' and 'estimate' is set to 'false').
%           'RRR'   : Relative redundancy reduction between inputs and
%           outputs
%     'classical' : Classical pattern separation measures with subfields:
%           'orthogonalisation': Orthogonalisation measure (increased
%           cosine distances)
%           'scaling': Scaling measure (increased differences in norms)
%           'decorrelation': Decorrelation measure
%           'hamm_dist': Hamming distance measure
%           'wass_dist': Wasserstein distance measure
%      'sparsity' : Numbers of spikes in input and output ensembles with
%      subfields:
%           'sparsity': Sparsity between input and output
%           'n_spike_in': Total number of spikes in the input ensemble
%           'n_spike_out': Mean number of spikes in each repetion of the output
%           ensemble
%       'info_details' : Details of the information theoretic pattern
%       separation calculations. Subfields are specific to the choice of
%       calculation method (for example, only transfer entropy has a delay/lag parameter,
%       and the Kozachenko-Leonenko estimators use smoothing parameters instead of
%       codes and bin sizes):
%            'MI': Raw mutual information
%            'input_code': Most informative spiking code for the input
%            ensemble
%            'output_code': Most informative spiking code for the output
%            ensemble
%            'best_binsize': Most informative bin size for the mutual
%            information or transfer entropy calculation
%            'best_lag': Most informative lag for the transfer entropy
%            calculation
%            'RR': Raw redundancy reduction
%            'input_redundancy': Redundancy of the input ensemble
%            'input_redundancy_code': Code for the redundancy of the input ensemble
%            'redundancy_binsize': Bin size for the redundancies
%            'output_redundancy': Redundancy of the output ensemble
%            'output_redundancy_code': Code for the redundancy of the output ensemble
% Example
% -------
% input_spiketimes = phase_locked(60,5,5,0.75,10); % Generate input ensemble
% output_spiketimes = spike_thinner_random(input_spiketimes,8,0.5); % Filter input ensemble
% patsep = analyse_pattern_separation(input_spiketimes,output_spiketimes);
% [patsep.info.spMI , patsep.info.RRR] % Print key information theoretic measures
%
% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird



function [patsep] = analyse_pattern_separation(input,output,varargin)
% Initialise outputs
patsep=struct; 
patsep.info = struct; % Initialise information-theoretic measures
patsep.classical = struct; % Initialise classical measures
patsep.sparsity = struct; % Initialise sparsity measures 
patsep.info_details = struct; % Initialise supporting details for informational measures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%% Validate and parse inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Default values of option name-value pairs
    defaultTE=false;
    defaultEST=false;
    defaultPAR=true;
    defaultCBS=10;
    defaultWS=5;
    
    % Check major input arguments are of the corrrect form
    if ~iscell(input)
        error('Input is not a cell array of spike times or voltage vectors.')
    end
    if size(input,2)~=1
        error('The second dimension of the input cell array should be one. Different components (spike trains) of the input pattern should use the first dimension.')
    end
    for train_ind=1:size(input,1)
        if ~ismatrix(input{train_ind}) || size(input{train_ind},1)>2 
            error('One or more elements of the input are not of the correct form (a 1xn vector of spiketimes or a 2xn matrix with time and voltage).')
        end    
    end
    if ~iscell(output)
        error('Output is not a cell array of spike times or voltage vectors.')
    end
    if size(output,1)~=size(input,1)
        warning('The first dimension of the output cell array does not match the first dimension of the input cell array, local rate and local temporal codes cannot be used.')
    end
    for train_ind=1:size(output,1)
        for rep_ind=1:size(output,2)
            if ~ismatrix(output{train_ind,rep_ind}) || size(output{train_ind,rep_ind},1)>2
                error('One or more elements of the output are not of the correct form (a 1xn vector of spiketimes or a 2xn matrix with time and voltage).')
            end
        end
    end

    % Construct and parse options
    p=inputParser;
    validEnsemble = @(x) iscell(x);
    addRequired(p,'input',validEnsemble);
    addRequired(p,'output',validEnsemble);
    validateTE = @(x) assert(islogical(x),'transfer_entropy should be set to true or false.');
    addParameter(p,'transfer_entropy',defaultTE,validateTE);
    validateEST = @(x) assert(islogical(x),'estimate should be set to true or false.');
    addParameter(p,'estimate',defaultEST,validateEST);
    validatePAR = @(x) assert(islogical(x),'parallel should be set to true or false.');
    addParameter(p,'parallel',defaultPAR,validatePAR);
    validateCBS = @(x) assert(isnumeric(x),'classical_binsize should be set to a numeric value in milliseconds.');
    addParameter(p,'classical_binsize',defaultCBS,validateCBS);
    validateMWS = @(x) assert(floor(x)==x && x>=1,'wordsize should be set to a positive integer.');
    addParameter(p,'wordsize',defaultWS,validateMWS);
    parse(p,input,output,varargin{:});
    
    % Map voltage traces to spike times if necessary
    input_spiketimes=cell(size(p.Results.input,1),1);
    for train_ind=1:size(p.Results.input,1)
        if size(p.Results.input{train_ind},1)==1 % Inputs are already spike times
            input_spiketimes{train_ind}=sort(p.Results.input{train_ind});
        elseif size(p.Results.input{train_ind},1)==2 % Inputs are voltage vectors and need to be converted to times
            input_spiketimes{train_ind}=voltage2times(p.Results.input{train_ind}); 
        end
    end
    output_spiketimes=cell(size(p.Results.output,1),size(p.Results.output,2));
    for train_ind=1:size(p.Results.output,1)
        for rep_ind=1:size(p.Results.output,2)
            if size(p.Results.output{train_ind,rep_ind},1)==1 % Outputs are already spike times
                output_spiketimes{train_ind,rep_ind}=sort(p.Results.output{train_ind,rep_ind});
            elseif size(p.Results.output{train_ind,rep_ind},1)==2 % Outputs are voltage vectors and need to be converted to times
                output_spiketimes{train_ind,rep_ind}=voltage2times(p.Results.output{train_ind,rep_ind});
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%% Information calculations %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mutual information or transfer entropy
if ~p.Results.transfer_entropy && ~p.Results.estimate % Calculate explicit mutual information
    num_params=struct;
    num_params.wordsize=p.Results.wordsize;
    if p.Results.parallel
        MIoptions='-max_bins -max_code -par'; % Use parallelisation
    else
        MIoptions='-max_bins -max_code'; % Don't use parallelisation
    end
    MI_obj=MI_function(input_spiketimes,output_spiketimes,[],num_params,MIoptions); % Calculate mutual information
    patsep.info_details.MI=MI_obj.MI;
    patsep.info_details.best_binsize=MI_obj.best_binsize;
    patsep.info_details.input_code=MI_obj.best_code{1};
    patsep.info_details.output_code=MI_obj.best_code{2};
elseif p.Results.transfer_entropy && ~p.Results.estimate % Calculate explicit transfer entropy
    num_params=struct;
    num_params.wordsize=p.Results.wordsize;
    if p.Results.parallel
        TEoptions='-max_bins -max_code -par'; % Use parallelisation
    else
        TEoptions='-max_bins -max_code'; % Don't use parallelisation
    end
    TE_obj=TE_function(input_spiketimes,output_spiketimes,[],num_params,TEoptions); % Calculate transfer entropy
    patsep.info_details.TE=TE_obj.TE;
    patsep.info_details.best_binsize=TE_obj.best_binsize;
    patsep.info_details.input_code=TE_obj.best_code{1};
    patsep.info_details.output_code=TE_obj.best_code{2};
    patsep.info_details.best_lag=TE_obj.best_delay;
elseif p.Results.estimate % Estimate mutual information
    if p.Results.parallel
        MIESToptions='-par'; % Use parallelisation
    else
        MIESToptions=''; % Don't use parallelisation
    end
    estMI_obj = estimate_MI(input_spiketimes,output_spiketimes,MIESToptions); % Estimate mutual information
    patsep.info_details.MI=estMI_obj.MI;
    patsep.info_details.MI_bounds=estMI_obj.MI_bounds;
end
% Redundancies
if ~p.Results.estimate % Calculate explicit redundancies
    num_params=struct;
    num_params.wordsize=p.Results.wordsize;
    if p.Results.parallel
        RRoptions='-max_bins -max_code -par'; % Use parallelisation
    else
        RRoptions='-max_bins -max_code'; % Don't use parallelisation
    end
    RR_obj=RR_function(input_spiketimes,output_spiketimes,[],num_params,RRoptions); % Calculate redundancies
    patsep.info_details.RR=RR_obj.RR;
    patsep.info_details.redundancy_binsize=RR_obj.best_binsize;
    patsep.info_details.input_redundancy_code=RR_obj.best_code{1};
    patsep.info_details.output_redundancy_code=RR_obj.best_code{2};
    patsep.info_details.input_redundancy=RR_obj.RRin;
    patsep.info_details.output_redundancy=RR_obj.RRout;
    if ~isfield(patsep.info_details,'MI') % Calculate mutual information for throughput if not already done
        MI_obj=MI_function(input_spiketimes,output_spiketimes,[],num_params,RRoptions); % Calculate mutual information
        patsep.info_details.MI=MI_obj.MI;
    end
elseif p.Results.estimate % Estimate redundancies
    if p.Results.parallel
        RRESToptions='-par'; % Use parallelisation
    else
        RRESToptions=''; % Don't use parallelisation
    end
    estRR_obj = estimate_RR(input_spiketimes,output_spiketimes,RRESToptions); % Estimate mutual information
    patsep.info_details.RR=estRR_obj.RR;
    patsep.info_details.RRin_bounds=estRR_obj.RRin_bounds;
    patsep.info_details.RRout_bounds=estRR_obj.RRout_bounds;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%% Sparsity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sparsity_obj=calculate_sparsity(input_spiketimes,output_spiketimes);
patsep.sparsity.sparsity=sparsity_obj.sparsity;
patsep.sparsity.n_spike_in=sparsity_obj.n_spike_in;
patsep.sparsity.n_spike_out=sparsity_obj.n_spike_out;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%% Classical measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classical_params=struct;
classical_params.binsize=p.Results.classical_binsize;
% Orthogonalsiation and scaling
orth_obj = orthogonalisation(input_spiketimes,output_spiketimes,classical_params);
patsep.classical.orthogonalisation=orth_obj.NDP_est;
patsep.classical.scaling=orth_obj.SF_est;
% Decorrelation
decorr_obj = decorrelation(input_spiketimes,output_spiketimes,classical_params);
patsep.classical.decorrelation=decorr_obj.decorr_est;
% Hamming distance
hamming_obj = hamming_sparsity(input_spiketimes,output_spiketimes,classical_params);
patsep.classical.hamm_dist=hamming_obj.hamming_est;
% Wasserstein distance
dist_obj = distance(input_spiketimes,output_spiketimes);
patsep.classical.wass_dist=dist_obj.dist_est;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%% Information-theoretic measures %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p.Results.transfer_entropy && ~p.Results.estimate
    patsep.info.spTE=patsep.sparsity.sparsity*patsep.info_details.TE;
else
    patsep.info.spMI=patsep.sparsity.sparsity*patsep.info_details.MI;
end
patsep.info.RRR=patsep.info_details.MI*patsep.info_details.RR;
end