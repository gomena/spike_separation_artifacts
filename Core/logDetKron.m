function [val grad]=logDetKron(Art,var,Difs,varsactive,Diags)
%  logDetKron computes the objective in equation 7, for the non-stimulating
%    electrodes
%  var: vector having all variables (for each kernel + variance): 1-3
%    correspond to time kernel, 4-6 to spatial kernel, 7-9 stimulus kernel
%    10 is variance.
%    in all cases first component= (log)lambda, 
%    second component = (log) alpha, third component= (log) beta
%  Difs: array of matrices with 
%    indices of active variables (the ones that will be differentiated)
%  Diag: array of diagonals of Difs, used to compute
%    non-stationarities
%  Art: Artifact
%  varsactive: the indexes of the variables that are active for optimization.
%  References: Thesis: Scalable Inference for Structured Gaussian 
%    Process Models. Yunus Saatchi, Cambride University
%Gonzalo Mena, 09/2017

logs=0;
types=[1 1 1]; 
nvar=[3 3 3];
varsmax=length(var);
nvarcum=cumsum([0 nvar]);
if(max(exp(var))==Inf)
    val = Inf;
    grad = Inf*ones(length(var),1);
    return
end

logs=0;
vec=[];
krondiag=1;
for k=1:3
    [Ker KerD]=evalKernels(Difs{k},Diags{k},var(nvarcum(k)+1:nvarcum(k+1)),types(k));
    [a b]=eig(Ker);
    Q{k}=a';
    Qt{k}=a;
    krondiag=kron(krondiag,diag(b));
    KerDs{k}=KerD;
    Kers{k}=Ker;
    KersQ{k}=diag(Q{k}*Kers{k}*Qt{k}); 
end

krondiag=exp(var(varsmax-1))*krondiag+exp(var(varsmax));
logs=logs+sum(log(krondiag));


prod1=KronProd(Q,Art);
krondiaginv=krondiag.^(-1);
prod1=krondiaginv.*prod1;
alpha=KronProd(Qt,prod1);
alpha2=dot(Art(:),alpha);
val=alpha2+logs;
if(nargout==1)
    
    return
end
for j=intersect(varsactive,[1:varsmax-1])
    g{j}=Kers;
    diagGQ{j}=KersQ;
end


for k=1:3
    for i=1:nvar(k)
        if(~isempty(intersect(varsactive,nvarcum(k)+i)))
            g{nvarcum(k)+i}{k}=exp(var(varsmax-1))*KerDs{k}{i};
            diagGQ{nvarcum(k)+i}{k}=exp(var(varsmax-1))*diag(Q{k}*KerDs{k}{i}*Qt{k});
        end
    end
end
if(any(varsmax-1==varsactive))
g{varsmax-1}{1}=g{varsmax-1}{1}*exp(var(varsmax-1));
diagGQ{varsmax-1}{1}=exp(var(varsmax-1))*diagGQ{varsmax-1}{1};
end

if(any(varsmax==varsactive))

g{varsmax}{1}=exp(var(varsmax))*speye(size(Kers{1},1));
g{varsmax}{2}=speye(size(Kers{2},1));
g{varsmax}{3}=speye(size(Kers{3},1));
diagGQ{varsmax}{1}=diag(exp(var(varsmax))*speye(size(Kers{1},1)));
diagGQ{varsmax}{2}=diag(speye(size(Kers{2},1)));
diagGQ{varsmax}{3}=diag(speye(size(Kers{3},1)));

end
for j=1:length(varsactive)
    varac=varsactive(j);
    d=1;
    
    for k=1:length(diagGQ{varac})
        d=kron(d,diagGQ{varac}{k});
    end
    grad(j)=dot(krondiaginv,d)-alpha'*KronProd(g{varac},alpha); 
end
grad=grad';

end