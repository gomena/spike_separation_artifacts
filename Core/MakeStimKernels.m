function params = MakeStimKernels(params,useNaiveExtrapolateStimElec,Diff,breakpoints)
% MakeStimKernels creates stimulating electrode kernels, by updating the
% params structure. Those kernels are stated in equation 3 of the paper.
% if useNaiveExtrapolateStimElec is true, then naive extrapolation is
% applied
% Diff is a matrix containing the differences between stimuli
% breakpoints contains the indexes of te breakpoints in the stimulating
% electrode
%Gonzalo Mena, 09/2017

%% loading covariate matrices for the stimulating electrode (equations 3,4,5,6)
Difs=params.patternInfo.Difs;
Diags=params.patternInfo.Diags;
Art=params.patternInfo.Art;
stimElecs=params.patternInfo.stimElecs;
var0=params.patternInfo.var0;

DifsStim{1}=Difs{1};
Difs{3}=Diff;
DifsStim{2}=zeros(size(Difs{3},1));
options=params.global.options;
params.patternInfo.breakpoints=breakpoints;

DiagsStim{1}=Diags{1};
DiagsStim{2}=breakpoints{1}';
if(params.global.extraStimElectrode+params.global.filterStimElectrode>=1)
    if(~useNaiveExtrapolateStimElec)
        for k=1:length(breakpoints{1})-1
            DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))=Difs{3}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1));
            if(DiagsStim{2}(k+1)-DiagsStim{2}(k)>1)
                DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))=DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))/max(max(DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))));
            end
        end
    else
        for k=1:length(breakpoints{1})-1
            l1= breakpoints{1}(k)+1:breakpoints{1}(k+1);
            l2=breakpoints{1}(k)+1:breakpoints{1}(k+1);
            m=max(length(l1),length(l2));
            for i=1:length(l1)
                for j=1:length(l2)
                    DifsStim{2}(l1(i),l2(j))=min(i,j)/m;
                end
            end
        end
    end
    
    %% Compute kernel parameters for each of the inter-breakpoint ranges
    for br=1:length(DiagsStim{2})-1
        inter=[DiagsStim{2}(br)+1:DiagsStim{2}(br+1)];
        DifsStimTemp=DifsStim;
        DifsStimTemp{2}=DifsStim{2}(inter,inter);
        DiagsStimTemp=DiagsStim;
        DiagsStimTemp{2}=DiagsStim{2}(br:br+1);
        types=[1 2];
        x01=[3 2 1];
        x02=-5;
        x03=15;
        nVar=[3 1];
        nVarCum=cumsum([0 nVar]);
        f2 = @(Art,x)logDetKronStimElec(squeeze(Art(inter,stimElecs(1),:)),[x log(var0)],DifsStimTemp,[1:5],types,DiagsStimTemp,nVar);
        factor=1;
        cont=1;
        xSti=-11*ones(1,5);
        stiOld=xSti(4);
        xStiHist=[];
        while(xSti(4)<-10||xSti(4)>5)
            g2 = @(x)f2(Art,x);
            if(length(inter)>2)
                xSti = fminunc(g2,[x01 factor*x02 x03],options); 
            else
                xSti=[x01 x02 x03];
            end
            evalsf(cont)=g2(xSti);
            xStiHist=[xStiHist;xSti];
            factor=factor/2;
            if(mod(cont,10)==0)
                factor=1;
            end
            cont=cont+1;
            if(cont==10)
                contind=find(evalsf==min(evalsf));
                xSti=xStiHist(contind(1),:);
                break
            end
        end
        for k=1:2
            [Ker KerD]=evalKernels(DifsStimTemp{k},DiagsStimTemp{k},xSti(nVarCum(k)+1:nVarCum(k+1)),types(k));
            [a b]=eig(Ker);
            QStim{br}{k}=a';
            QtStim{br}{k}=a;
            dLStim{br}{k}=diag(b);
            KersStim{br}{k}=Ker;
            xStim(br,:)=xSti;
        end
    end
    %% save the obtained parameters: Kernels (kersstim), their corresponding eigendecomposition dLStim,QtStim,QStim, and the hyperparameters (xstim)
    params.patternInfo.KersStim=KersStim;
    params.patternInfo.dLStim=dLStim;
    params.patternInfo.QtStim=QtStim;
    params.patternInfo.QStim=QStim;
    params.patternInfo.xStim = xStim;
end