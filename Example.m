addpath('Core')
addpath('Display')
addpath('Utils')
%% First compute array-wise artifact hyperparameters (equation 7), load model parameters and store them in params structure

%Art0 Contains an artifact 'pseudosample' constructed by translating mean of traces
% to stimulation at different stimulating electrodes (as described in the 'Artifact Learning' section')
% "Artifact Learning" section
ArtStr='Art0.mat';

% arrayObj MEA-specific parameters (e.g. positions of electrodes, number of
% electrodes, etc). Array(512) returns those parameters for the
% 512-electrode stimulation system.
arrayObj=Array(512);

%InitializeArrayRobust for details
[params]=InitializeArrayRobust(arrayObj,ArtStr);


%% The params are now stored in the params structure. They can be modified
%here if needed.
% e.g. 
params.global.tarray=[0 [7:33]]; %expanded the set to look for spikes a bit

%see below for a comprehensive list of values in the params structure)


%% Now consider the Data (actual traces, stimulating electrode information, breakpoint list and list of stimulus is stored here).
DataStr='Data.mat';
load(DataStr)
%% Do spike Sorting and Store it in Output structure. Notice the optional fields (e.g. 'useNaiveExtrapolate') see the code in DoSpikeSortingLargeScale for details
[Output,params]=DoSpikeSortingLargeScale(params,arrayObj,DataStr,'useNaiveExtrapolate',0);
%Here I updated the params structure for pedagogical reasons. This may not be necessary in practice.
%% Output.Log indicates the number of coordinate ascent iterations per stimulus amplitude

%% Output.neuronInfo
%templates: templates
%ActiveElectrodes: array containing the 'useful' electrodes for each neuron (i.e. the ones that have a signal >thresEI)
%ActiveElectrodesAll = all useful electrodes collapsed
%spikes: spike information (cell array, each element is a J*max(nTrials) matrix containing spike times (0 if no spike) for each all repetitions at each of the J amplitude of stimulation

%% Output.stimInfo
% ActiveElectrodes: active electrodes (for the Array(501) case it will  always be the vector [1:512];
% breakpoints: breakpoint amplitudes (the 0-th amplitude and last one are considered breakpoints) 
% Kers: cell array with the non-stimulating electrode Kernels. Kers{1}: temporal kernel, Kers{2}: spatial Kernel (numbered by electrode number, does not consider stimulating electrode), Kers{3}=amplitude Kernels {[55x55 double]  [511x511 double]  [30x30 double]}
% KersStim: cell array with the stimulating electrode Kernels (one for each inter-breakpoint range). Each of the cells in turn another cell array containing the time (first component) and spatial (second component) kernels, respectively
% xStim: [3x5 double]
% stimElecs: stimulating electrode
% patternNo: pattern number. Here, it is the same as stimulating electrode.
% listAmps: list of amplitudes 
% var0: estimate of variance sigma^2
% Arts: Final stimate of the artifact (After subtracting spikes and filtering)
% nTrials: vector containing the number of trials (repetitions) for each stimulation amplitude

%% params.global
%useStimElectrode: use stim electrode? if zero only non-stimulating electrodes are considered
%filterStimElectrode: if 1 apply the filtering operation (see validation) to the stimulating electrode. if not estimates of the artifact are the simple average over traces.
%extraStimElectrode: if 0 uses the naive extrapolation (see validation section) if 1 uses kernel-based extrapolation. For the stimualting electrode.
%filterNoStimElectrodes:  if 1 apply the filtering operation (see validation) to the non-stimulating electrode. if not estimates of the artifact are the simple average over traces.
%extraNoStimElectrodes: if 0 uses the naive extrapolation (see validation section) if 1 uses kernel-based extrapolation. For the non-stimualting electrodes
%Tmax: number of time-steps (sampling rate=20KhZ)
%tarray: set of time-steps to look for spikes
%options: optimization options
%positions: Load the electrode positions
%maxIter: maximum number of coordinate ascent iterations per amplitude of stimulation
%thresEI: EI strenght treshold: for spike inference, only electrodes such that for some neuron the EI at that electrode exceeds the value thresEI are considered for spike sorting.

%% params.arrayInfo : see initializeArrayRobust

%% params.patternInfo: contains more structure but are non-relevant (see the code for their explanation)



%% Now plot the artifact at the stimulation amplitude n25 (2.01 mA)

center=404; %relative center of the array
patternNo=404; %stimulating electrode (plotted on a different scale)
size=[0 0 0.8 0.4]; %size of the window 
factor=0.2;
colors{1}=[0 0 1]; %color for the stimulating electrode
colors{2}=[0 0 0.5]; %color for the non-stimulating electrode
axiss{2}=[1 20 -50 400]; %axis for the non-stimulating electrode
axiss{1}=[1 20 -700 200]; %axis for the stimualting electrode
nh=10; nv=10; % neihborhoods in the x and y axis)
figuren=1; %figure 1
PlotArtifact(Output.stimInfo.Arts,25,center,patternNo,nh,nv,figuren,size,axiss,factor,colors)

%% Finally plot the activation curves of each neuron (notice only neuron 17 was activated). Insets show the electrical image (EI) of each neuron in relation to the stimulating electrode (yellow dot). The soma of the activated neuron is the closest to the stimulating electrode.
hf=figure(2);
set(hf,'units','normalized');
set(hf,'position',[0 0 1 0.6])
for i=1:24
    
    
    h2=subplot(3,8,i);
    plot(listAmps,nanmean(Output.neuronInfo.spikes{i}'>1),'linewidth',2,'color','black')
    hold on
    scatter(listAmps,nanmean(Output.neuronInfo.spikes{i}'>1),20,[0 0 1],'filled')
    
    hold on
    axis([0 max(listAmps)+0.0001 0 1])
    title(['Neuron ' num2str(i)])
    set(gca,'fontsize',13)
    set(gca,'TickLength',[0.08 0.08])
    
    pos=get(h2,'position');
    yminnew=pos(2)+pos(4)/2+pos(4)/32;
    h3=axes('position',[pos(1)+pos(3)/64 yminnew pos(3)/2 pos(4)/2-pos(4)/8]);
    mat=nanmax(abs(templates{i})');
    mat(Output.stimInfo.stimElecs)=100;
    [~,EIm_view]   = ei2matrix(log(mat)');
    imagesc(EIm_view,[0 log(100)])
    box on
    axis 'off'
    axis square
end