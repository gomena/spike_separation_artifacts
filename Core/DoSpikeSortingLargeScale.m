function [Output,params]=DoSpikeSortingLargeScale(params,arrayObj,DataStr,varargin)
%DoSpikeSortingLargeScale is a wrapper of SpikeSortingAllCases
%  SpikeSortingAllCases; it performs spike identification to responses to
%    electrical stimulation at a single electrical stimulus
%  params: structure with details on parameters.
%  arrayObj: structure containing details for array
%  DataStr: direction to a .mat file that contains all relevant data
%    (traces, list of stimulating electrodes, etc).
% Gonzalo Mena, 09/2017

disp('Now loading relevant data structures for the stimulating electrode')
%% First load optional parameters:
load(DataStr)
useNaiveExtrapolate=0;
useNaiveExtrapolateStimElec=0;
useStimElectrode=1;
filterStimElectrode=1;
filterNoStimElectrode=1;
% Read in optional inputs.
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'usenaiveextrapolate'
            useNaiveExtrapolate= varargin{j*2};
        case 'usenaiveextrapolatestimelec'
            useNaiveExtrapolateStimElec= varargin{j*2};
        case 'usestimelectrode'
            useStimElectrode= varargin{j*2};
        case 'filternostimelectrode'
            filterStimElectrode=varargin{j*2};
        case 'filterstimelectrode'
            filterNoStimElectrode=varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

%positions=params.global.positions;
x=params.arrayInfo.x;
params.global.useStimElectrode = useStimElectrode; %use stim electrode? if zero only non-stimulating electrodes are considered
params.global.filterStimElectrode= filterStimElectrode; %if 1 apply the filtering operation (see validation) to the stimulating electrode. if not estimates of the artifact are the simple average over traces.
params.global.extraStimElectrode=1-useNaiveExtrapolateStimElec;%if 0 uses the naive extrapolation (see validation section) if 1 uses kernel-based extrapolation. For the stimualting electrode.
params.global.filterNoStimElectrodes= filterNoStimElectrode; %if 1 apply the filtering operation (see validation) to the non-stimulating electrode. if not estimates of the artifact are the simple average over traces.
params.global.extraNoStimElectrodes = 1 -useNaiveExtrapolate; %if 0 uses the naive extrapolation (see validation section) if 1 uses kernel-based extrapolation. For the non-stimualting electrodes

Tmax=params.global.Tmax;
var0=params.arrayInfo.var0;


% Get active electrode list from the array object (see initializeArrayRobust) and distance between electrodes and stimulating electrodes)
activeElecs = arrayObj.getElectrodes;
numActiveElecs = length(activeElecs);
positions=arrayObj.getPositions;
[~,rho] = cart2pol(positions(:,1)-positions(stimElecs,1),positions(:,2)-positions(stimElecs,2));
patternNo=stimElecs;
ind=setdiff(activeElecs,find(rho==0));
for e=1:length(stimElecs)
    stimElecsRel(e)=find(activeElecs==stimElecs(e));
end

params.patternInfo.ind=ind;
params.patternInfo.stimElecs=stimElecs;


%% construct covariance matrices for the non-stimulating electrodes (Equations 2,4,5,6)

Dif1= zeros(Tmax);
for i=1:Tmax
    for j=1:Tmax
        Dif1(i,j)=abs(i-j);
    end
end

if(useNaiveExtrapolate==0)
    Dif3 = zeros(size(listAmps,1),size(listAmps,1));
    for j=1:length(listAmps)
        for i=1:length(listAmps)
            Dif3(i,j)=abs(listAmps(i)-listAmps(j));
        end
    end
else
    for j=1:length(listAmps)
        for i=1:length(listAmps)
            Dif3(i,j)=min(listAmps(i),listAmps(j));
        end
    end
end
rho=rho(ind);

Difs{2}=arrayObj.difPos(ind,ind)/arrayObj.maxR;
Difs{3}=Dif3;
Difs{3}=Difs{3}/max(max(listAmps));
Dif1=Dif1/Tmax;
Difs{1}=Dif1;


Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho/arrayObj.maxR;
Diags{3}=listAmps(:,1)/max(listAmps);



%% evaluate Kernels given the constructed covariate matrices

types=[1 1 1];
factp=[0 3 6 9];

for k=1:3
    [Ker, KerD]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),types(k));
    Kers{k}=Ker;
    [a, b]=eig(Kers{k});
    Q{k}=a';
    Qt{k}=a;
    dL{k}=diag(b);
    
end
params.arrayInfo.x=x;
%% Expand params structure, including stimulation-pattern-specific fields.
%Q,Q,dL is the eigendecomposition of Ker


params.patternInfo.Q=Q;
params.patternInfo.Qt=Qt;
params.patternInfo.Kers=Kers;
params.patternInfo.dL=dL;

params.patternInfo.Art=squeeze(nanmean(TracesAll,2));
params.patternInfo.var0=var0;
params.patternInfo.rho=rho;
params.patternInfo.listAmps=listAmps;
params.patternInfo.Diags=Diags;
params.patternInfo.Difs=Difs;
params.patternInfo.breakpoints=breakpoints;
params.neuronInfo.templates = templates;

%% Compute stimulating electrode hyperparameters (these depend on stimulating electrode, therefore cannot be 'initialized')
if(useStimElectrode)
    disp('Now computing stimulating electrode hyperparameters')
    if(params.global.filterStimElectrode+params.global.extraStimElectrode>=1)
        
        Diff = zeros(size(listAmps,1),size(listAmps,1));
        for j=1:length(listAmps)
            for i=1:length(listAmps)
                Diff(i,j)=abs(listAmps(i)-listAmps(j));
            end
        end
        
        Diff=Diff/max(max(listAmps));
        params = MakeStimKernels(params,useNaiveExtrapolateStimElec,Diff,breakpoints);
        disp('Stimulating electrode hyperparameters were succesfully found')
        
    end
end

%% Do spike sorting
[spikes, Log, params]=SpikeSortingAllCases(params,TracesAll);

if(useStimElectrode>=1)
    Output.stimInfo.ActiveElectrodes=activeElecs;
    Output.stimInfo.breakpoints=params.patternInfo.breakpoints;
    if(params.global.extraStimElectrode)
        Output.stimInfo.KersStim=params.patternInfo.KersStim;
    end
else
    Output.stimInfo.ActiveElectrodes=activeElecs(ind);
end


Output.neuronInfo=params.neuronInfo;
Output.neuronInfo.spikes=spikes;
Output.stimInfo.stimElecs=stimElecs;
Output.stimInfo.listAmps=listAmps;
Output.stimInfo.var0=var0;
Output.stimInfo.Arts=params.patternInfo.Arts;
Output.stimInfo.nTrials=params.patternInfo.nTrials;
Output.stimInfo.Kers=params.patternInfo.Kers;
Output.stimInfo.rho=rho;
Output.Log=Log;
params = rmfield(params,'neuronInfo');
end
