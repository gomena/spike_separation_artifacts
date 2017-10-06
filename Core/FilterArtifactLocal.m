function ArtF=FilterArtifactLocal(Kers,Art,x,cmax,ind,Q,Qt,krondiaginv)
%  FilterArtifactLocal applies equation 10, i.e. it
%    filters the Artifact (mean of spike substracted)
%    traces accross repetitions of the same stimulus amplitudes determined by
%    cmax) using the  model and noise level. Use for all electrodes 
%    but the stimulating.
%  x is a vector of parameters
%  L,Q,Qt represent the eigenvalue decomposition of Kers (cell arrays)
%  Here, x(end) contains log variance sigma0.
%  ind contains indices of the non-stimulating electrodes
%  krondiaginv is the vectorized kronecker product of the inverses of the
%    eigenvalues of Kers
%  Art is the Artifact

%Gonzalo Mena, 09/2017.


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
ArtF=reshape(Filtered,size(Art,1),length(ind),size(Art,3));
