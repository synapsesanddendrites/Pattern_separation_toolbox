% 'start_pattern_separation' Script to check if pattern separation
% dependencies are installed on this system and add necessary folders to
% the path.
% 
% (pattern separation toolbox)
%
% Example
% -------
% start_pattern_separation % Check dependencies are installed and add necessary subfolders to the path
%
% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Optimization toolbox %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_Matlab=ver; % Check Matlab installation
if ~license('test','Optimization_Toolbox')
    warning('The license for the MATLAB (local) optimization toolbox is missing. Please add this toolbox to your Matlab installation to use the pattern separation toolbox.')
end
opt_found=false;
for ward=1:length(full_Matlab)
    if strcmp(full_Matlab(ward).Name,'Optimization Toolbox')
        opt_found=true;
    end
end
if ~opt_found
    warning('The MATLAB (local) optimization toolbox is not installed on your system, but you may have a license for it. Please add this toolbox to your Matlab installation to use the pattern separation toolbox.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Global optimization toolbox %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~license('test','GADS_Toolbox')
    warning('The license for the MATLAB global optimization toolbox is missing. Please add this toolbox to your Matlab installation to use the pattern separation toolbox.')
end
opt_found=false;
for ward=1:length(full_Matlab)
    if strcmp(full_Matlab(ward).Name,'Global Optimization Toolbox')
        opt_found=true;
    end
end
if ~opt_found
    warning('The MATLAB global optimization toolbox is not installed on your system, but you may have a license for it. Please add this toolbox to your Matlab installation to use the pattern separation toolbox.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Parallel computing toolbox %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~license('test','Distrib_Computing_Toolbox')
    warning('The license for the MATLAB parallel toolbox is missing. Without this toolbox, the parallelised features of the pattern separation toolbox will not work. Please set the option "parallel" in the function "analyse_pattern_separation" to "false", or add the parallel computing toolbox to your Matlab installation.')
end
opt_found=false;
for ward=1:length(full_Matlab)
    if strcmp(full_Matlab(ward).Name,'Parallel Computing Toolbox')
        opt_found=true;
    end
end
if ~opt_found
    warning('The MATLAB parallel toolbox is not installed on your system, but you may have a license for it. Without this toolbox, the parallelised features of the pattern separation toolbox will not work. Please set the option "parallel" in the function "analyse_pattern_separation" to "false", or add the parallel computing toolbox to your Matlab installation.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Curve fitting toolbox %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~license('test','Curve_Fitting_Toolbox')
    warning('The license for the MATLAB curve fitting toolbox is missing. Without this toolbox, the estimators for the mutual information and redundancy in the pattern separation toolbox cannot be used (but exact measures can still be computed). If you would like to use these estimators, please add the curve fitting toolbox to your Matlab installation.')
end
opt_found=false;
for ward=1:length(full_Matlab)
    if strcmp(full_Matlab(ward).Name,'Curve Fitting Toolbox')
        opt_found=true;
    end
end
if ~opt_found
    warning('The MATLAB curve fitting toolbox is not installed on your system, but you may have a license for it. Without this toolbox, the estimators for the mutual information and redundancy in the pattern separation toolbox cannot be used. If you would like to use these estimators, please add the curve fitting toolbox to your Matlab installation.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Signal processing toolbox %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~license('test','Signal_Toolbox')
    warning('The license for the MATLAB signal processing toolbox is missing. Without this toolbox, the pattern separation toolbow will not be able to process voltage traces. If you need to analyse voltage traces (instead of spike times), please add the signal processing toolbox to your Matlab installation.')
end
opt_found=false;
for ward=1:length(full_Matlab)
    if strcmp(full_Matlab(ward).Name,'Signal Processing Toolbox')
        opt_found=true;
    end
end
if ~opt_found
    warning('The MATLAB signal processing toolbox is not installed on your system, but you may have a license for it. Without this toolbox, the pattern separation toolbow will not be able to process voltage traces. If you need to analyse voltage traces (instead of spike times), please add the signal processing toolbox to your Matlab installation.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Add necessary folders to path %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('generate_spike_ensembles','info_theory_measures','other_measures','prepare_spikes','utility_functions','tutorials','example_data')
