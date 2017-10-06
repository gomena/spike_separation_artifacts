function ArtF=FilterArtifactLocalStim(Kers,Art,x,cmax,Q,Qt,krondiaginv)
%  FilterArtifactLocalStim applies equation 10, i.e. it
%    filters the Artifact (mean of spike substracted) traces accross 
%    repetitions of the same stimulus amplitudes determined by cmax using 
%    the model and noise level. Use for stimulating electrode
%  x is a vector of parameters
%  L,Q,Qt represent the eigenvalue decomposition of Kers (cell arrays)
%  Here, x(end) contains log variance sigma0.
%  krondiaginv is the vectorized kronecker product of the inverses of the
%    eigenvalues of Kers.
%  Art is the artifact

%Gonzalo Mena, 09/17
varsmax=length(x);
a=Kers{2}(cmax,cmax);
Kers{2}=1;
Art=Art(cmax,:);
Q{2}=1;
Qt{2}=1;

prod1=KronProd(Q,Art);
prod1=krondiaginv.*prod1;

alpha=KronProd(Qt,prod1);
Filtered=a*exp(x(varsmax-1))*KronProd(Kers,alpha);
ArtF=reshape(Filtered,size(Art,1),size(Art,2));
