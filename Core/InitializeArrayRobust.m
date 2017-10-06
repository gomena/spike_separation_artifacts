function [params]=InitializeArrayRobust(arrayObj,ArtStr,varargin)
% Optional inputs:
%           numElecs: 512 or 519 (default 512)
%
% Gonzalo Mena, 03/2016

% Set default params
load(ArtStr)


disp('Computing initialization hyperparameters')
nneighbors=5; %optional parameter: degree of the neighbors (in the metric induced by vicinity relations of electrodes in the array) are used for hyperparameter estimation? if nneighbor is too large computations will be large and the model will fit noise.

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

for jj=1:(nbin/2)
    if ~ischar(varargin{jj*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{jj*2-1})
        
        case 'nneighbors'
            nneighbors=varargin{jj*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end



%% Set default parameters

params.global.Tmax=55; %number of time-steps (sampling rate=20KhZ)
params.global.tarray=[0 [7:32]]; %set of time-steps to look for spikes
params.global.options=optimoptions('fminunc','Algorithm','trust-region','GradObj','on'); %optimization options
positions = arrayObj.getPositions; % Load the electrode positions
params.global.positions=positions;
params.global.maxIter=5; %maximum number of coordinate ascent iterations per amplitude of stimulation
params.global.thresEI=30; %EI strenght treshold: for spike inference, only electrodes such that for some neuron the EI at that electrode exceeds the value thresEI are considered for spike sorting.



tarray=params.global.tarray;
Tmax=params.global.Tmax;
options=params.global.options;



%% Load covariate matrices to represent localization(equations 2,4,5,6)

% Get active electrode list from the array object
activeElecs = arrayObj.getElectrodes;
DifPos=arrayObj.difPos;


pattern0 = arrayObj.center;
% distances with the stimulating electrode
[theta,rho]  = cart2pol(positions(:,1)-positions(pattern0,1),positions(:,2)-positions(pattern0,2));

rhoUnique=rho(els);
thetaUnique=theta(els);
sizes=[];
contTraces=1;


els2=intersect(els,arrayObj.getNeighbors(pattern0,nneighbors));
maxCond=length(listAmpsAll);
for l=1:maxCond
    for f=1:length(els2)
        f2=find(els2(f)==els);
        Arts(l,f,:)=Arts0(l,f2,:);
    end
end

Dif1= zeros(Tmax,Tmax);
for jj=1:Tmax
    for ii=1:Tmax
        Dif1(ii,jj)=abs(ii-jj);
    end
end


Dif1=Dif1/Tmax;
Difs{1}=Dif1;
Difs{2}=DifPos(els2,els2)/arrayObj.maxR;


listAmpsAll=unique(listAmpsAll);

Dif3 = zeros(length(listAmpsAll),length(listAmpsAll));
for jj=1:length(listAmpsAll)
    for ii=1:length(listAmpsAll)
        Dif3(ii,jj)=abs(listAmpsAll(ii)-listAmpsAll(jj));
    end
end


DifPos=arrayObj.difPos;

Dif3=Dif3(1:maxCond,1:maxCond)/max(listAmpsAll);
Difs{3}=Dif3;
Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho(els2)/arrayObj.maxR;
Diags{3}=listAmpsAll(1:maxCond)/max(listAmpsAll);



%% estimate variance (phi^2 in equation 7)
for i1=1:size(Arts,1);
    a=Arts(i1,:,:);
    vars(i1)=nanvar(a(:));
end
var0=nanmean(vars(1:5));



%% compute hyperparameters

f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1:7) -100 -100  x(end) log(var0)],Difs,setdiff([1:10],[8 9]),Diags);
g1=@(x)f1(Arts,x);
x1 = fminunc(g1,[-2 0 0 -2 0 0 -2 10],options);
x =[x1(1:7) -100 -100  x1(end)];
%x(1:3)=time parameters (log)
%x(4:6)=space parameters (log)
%x(7:9)=amplitude parameters (log)
%x(10)=log(var0);
%in all cases first component= (log)lambda, second component = (log) alpha, third component= (log) beta
%for example, x(8)=x(9)=-100 since localization is not imposed for
%amplitude-wise kernel (therefore, alpha=beta=0).



%% save remaining parameters
params.arrayInfo.var0=var0;
params.arrayInfo.x= x;
disp('Initialization hyperparameters were succesfully computed')
