function [Apred]=ExtrapolateArtifactCondStim(Kers,Q,Qt,dL,cmax,Art,x,sigma2)
%Extrapolate artifact at condition cmax given (full) observations until
%condition cmax-1. Extrapolation is done in al electrodes but the
%stimulating
%Gonzalo Mena, 03/2016


T=size(Kers{1},1);
Kers{2}=Kers{2}(1:cmax,1:cmax);
ctest=cmax;
test{1}=[1:T];
test{2}=ctest;


train{1}=[1:T];
train{2}=setdiff([1:size(Kers{2},1)],ctest);


Ktrain{1}=Kers{1};
Ktrain{2}=Kers{2}(train{2},train{2});


Ktraintest{1}=Kers{1};
Ktraintest{2}=Kers{2}(test{2},train{2});



krondiag=1;
for k=1:2
if(k==2)
    [a b]=eig(Ktrain{k});
    Q{k}=a';
    Qt{k}=a;
    dL{k}=diag(b);
end

krondiag=kron(krondiag,dL{k});


end

krondiag=exp(x(end))*krondiag+sigma2;

prod1=KronProd(Q,Art);
krondiaginv=krondiag.^(-1);
prod1=krondiaginv.*prod1;
alpha=KronProd(Qt,prod1);

Apred=exp(x(end))*KronProd(Ktraintest,alpha);

Apred=reshape(Apred,length(ctest),T);
