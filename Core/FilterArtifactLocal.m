function ArtF=FilterArtifactLocal(Kers,Art,x,cmax,ind,Q,Qt,krondiaginv)
%Filter the Artifact (mean of spike substracted traces accross repetitions
%of the same stimulus amplitudes) using the model and noise level. Use for
%all electrodes but the stimulating
%Here, x(end) contains log variance sigma0.
%Gonzalo Mena, 03/2016


varsmax=length(x);
logs=0;

Kers{3}=Kers{3}(cmax,cmax);
Art=Art(cmax,:,:);
Q{3}=1;
Qt{3}=1;

prod1=KronProd(Q,Art);
prod1=krondiaginv.*prod1;
alpha=KronProd(Qt,prod1);


Filtered=exp(x(varsmax-1))*KronProd(Kers,alpha);
%Filtered=KronProd(Kers,alpha);

ArtF=reshape(Filtered,size(Art,1),length(ind),size(Art,3));
