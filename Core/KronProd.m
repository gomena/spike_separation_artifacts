function [Z]=KronProd(K,y)
%computes A*y where K is the kronecker product between K{1},K{2},...K{n}
%Gonzalo Mena, 03/2016
%References:Elad Gilboa, Yunus Saatc ?and John P. Cunningham 
%Scaling Multidimensional Inference for Structured Gaussian Processes
% IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, VOL. 37, NO. 2, FEBRUARY 2015
%
Ntotal=prod(size(y));
D=length(K);
Z=y;

for i=1:D

    d=D-i+1;
    Gd=size(K{d},2);
    
X=reshape(Z,Gd,prod(size(Z))/Gd);

Z=K{d}*X;

Z=Z';
Z=Z(:);

end
