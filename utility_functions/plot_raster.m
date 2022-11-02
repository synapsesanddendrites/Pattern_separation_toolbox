% 'plot_raster' Plotting function to visualise spike train ensembles.
% 
% (pattern separation toolbox)
%
% plot_raster(spike_ensemble)
% --------------------------------------------------------------
%
% Input
% -----
% - spike_ensemble: Nx1 cell array of input spike times or voltage vectors, where N
% is the number of spike trains in the input ensemble. If inputs are spike
% times, each cell should be a 1xn vector (where n is the number of spikes in 
% that train). If inputs are voltage traces, each cell should be a 2xn
% matrix, where the first row is time, the second row is voltage, and n is the
% number of timesteps recorded.
%
% - margin (optional): Relative distance between trains in the plot (default=0.025).
% 
% - filename (optional): Where to save the figure if required (default
% 'None', no save). If a filepath is specified, the raster is saved as both
% a .png and .eps.
% 
% - fig_size (optional): Size (in cm) of saved figure (default [4,4]).
%
% Example
% -------
% input_spiketimes = phase_locked(60,5,5,0.75,10); % Generate input ensemble
% plot_raster(input_spiketimes) % Plot raster
%
% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird


function [] = plot_raster(spike_ensemble,margin,filename,fig_size)
% PLOT_RASTER Plots a raster of spikes in an ensemble
if nargin<2 || isempty(margin)
    margin=0.025;
end
if nargin<3 || isempty(filename)
    filename='None';
end
if nargin<4 || isempty(fig_size)
    fig_size=[4,4];
end

spikes=cell(size(spike_ensemble,1),1);
for train_ind=1:size(spike_ensemble,1)
    if size(spike_ensemble{train_ind},1)==1 % Inputs are already spike times
        spikes{train_ind}=sort(spike_ensemble{train_ind});
    elseif size(spike_ensemble{train_ind},1)==2 % Inputs are voltage vectors and need to be converted to times
        spikes{train_ind}=voltage2times(spike_ensemble{train_ind});
    end
end

t_max=0;
max_spike=0;
n_train=size(spikes,1);
for train_ind=1:n_train
    try
        if spikes{train_ind}(end)>t_max
            t_max=spikes{train_ind}(end);
        end
        if length(spikes{train_ind})>max_spike
            max_spike=length(spikes{train_ind});
        end
    catch
    end
end


hstep=(1-margin)/n_train;
figure
hold on
for ward=1:n_train
    n_spike=length(spikes{ward});
    for sooth=1:n_spike
        tspike=spikes{ward}(sooth)/((1+margin)*t_max);
        plot((margin/2+tspike)*ones(10,1),linspace(margin/2+(ward-1)*hstep+margin*hstep,margin/2+ward*hstep-margin*hstep,10),'Color','black','LineWidth',0.5)
    end
end
set(gca,'xlim',[0,1])
set(gca,'ylim',[0,1])
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'box','on')
set(gca,'DataAspectRatio',[1 1 1])

if ~strcmp(filename,'None')
   set(gcf,'PaperUnits','centimeters');
   set(gcf,'PaperPosition',[0 0 fig_size(1) fig_size(2)])
   exportgraphics(gcf,strcat(filename,'.png'),'Resolution',300)
   saveas(gcf,filename,'eps')
end
end


