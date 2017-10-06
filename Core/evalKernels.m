function [Ker KerD]=evalKernels(Dif,x,vars,type,varargin)
%  evalKernels computes the Kernel(Ker) and its gradient (KerD) for each of
%    the three params.
%  Dif: diferences between covariates (matrix)
%  x: on-stationarity vector (diagonal of Dif usually)
%    (also, for stimulating electrode, indicate
%    breakpoints)
%  vars: values of the variables, depends on the type Kernel
%    (type=1 for regular non-stimulating electrodes kernels(equation 2),
%    type=2 to encode information of stimulating electrode (equation 3)
% Gonzalo Mena, 09/2017

if(type==1)
    %non-stimulating electrodes
    lambda=exp(vars(1));
    beta=exp(vars(2));
    alpha=exp(vars(3));
    K=1;
    Dg=K*exp(-x*beta).*x.^alpha;
    DgK=Dg;
    Dga=DgK.*alpha.*log(x);
    Dga(isinf(Dga))=0;
    Dga(isnan(Dga))=0;
    Dgb=-DgK*beta.*x;
    if(lambda==Inf)
        Ker0=eye(size(Dif,1));
    else
        Ker0=(1+sqrt(3)*(Dif)*lambda).*exp(-sqrt(3)*(Dif)*lambda);
    end
    Ker=(Dg*Dg').*Ker0;
    Ker0l=-3*lambda*exp(-sqrt(3)*(Dif)*lambda).*lambda.*(Dif).^2;
    KerD{1}=(Dg*Dg').*Ker0l;
    KerD{2}=(Dg*Dgb'+Dgb*Dg').*Ker0;
    KerD{3}=(Dg*Dga'+Dga*Dg').*Ker0;
    
elseif(type==2)
    %for Stimulating electrode, include information of breakpoints in x
    lambda=exp(vars);
    lambdas=zeros(size(Dif,1));
    ind=zeros(size(Dif,1));
    for k=1:length(lambda)
        ma=zeros(size(Dif,1));
        ma(1:x(2)-x(1),1:x(2)-x(1))=lambda(k);
        lambdas=lambdas+ma;
        ind(1:x(2)-x(1),1:x(2)-x(1))=1;
        KerD{k}=-3*ma.*exp(-sqrt(3)*(Dif).*ma).*ma.*(Dif).^2;
    end
    Ker=(ind+sqrt(3)*(Dif).*lambdas).*exp(-sqrt(3)*(Dif).*lambdas);
end


