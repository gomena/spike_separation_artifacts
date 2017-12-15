function [ActionPotential] = makeActionPotential(nNeuron,spikeTime,templates,T)

offset0 = 12;
E = size(templates{1},1);
ActionPotential=zeros(E,T);

if(spikeTime == 0)
    return
end

minindex = offset0-spikeTime+1;
if(minindex >= 1)
    ActionPotentialAux = templates{nNeuron}(:,offset0-spikeTime+1:end);
else
    ActionPotentialAux = [zeros(E,1-minindex) templates{nNeuron}];
end

if(size(ActionPotentialAux,2)>=T)
    ActionPotential = ActionPotentialAux(:,1:T);
else
    ActionPotential = [ActionPotentialAux repmat(ActionPotentialAux(:,end),1,T-size(ActionPotentialAux,2))];
end
