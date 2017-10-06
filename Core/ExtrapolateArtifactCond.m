function [Apred]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,cmax,Art,x,sigma0)
%Extrapolate artifact at condition cmax given (full) observations until
%condition cmax-1. Extrapolation is done in al electrodes but the
%stimulating
%Gonzalo Mena, 03/2016


%extrapolate
ctest=cmax;
test{1}=[1:size(Art,3)];
test{2}=[1:size(Art,2)];
test{3}=ctest;


train{1}=[1:size(Art,3)];
train{2}=[1:size(Art,2)];
train{3}=setdiff([1:size(Art,1)],ctest);



Ktrain{1}=Kers{1};
Ktrain{2}=Kers{2};

Ktrain{3}=Kers{3}(train{3},train{3});


Ktraintest{1}=Kers{1};
Ktraintest{2}=Kers{2};
Ktraintest{3}=Kers{3}(test{3},train{3});

krondiag=1;
for k=1:3
    if(k==3)
        [a b]=eig(Ktrain{k});
        Q{k}=a';
        Qt{k}=a;
        dL{k}=diag(b);
    end
    
    krondiag=kron(krondiag,dL{k});
    
    
end

krondiag=krondiag+sigma0*exp(-x(end));
prod1=KronProd(Q,Art(train{3},:,:));
krondiaginv=krondiag.^(-1);
prod1=krondiaginv.*prod1;
alpha=KronProd(Qt,prod1);

Apred=KronProd(Ktraintest,alpha);
Apred=reshape(Apred,length(ctest),size(Art,2),size(Art,3));

