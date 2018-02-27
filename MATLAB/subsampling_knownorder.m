function [PrMAT]=subsampling_knownorder(X,Y,R,alpha,LASSO,RIDGE,MAXITER,CONV,INIT,nsamples,splitratio,asl,obl,qlambda,PTRUE)
%SUBSAMPLING Stability selection for sparse covariates regression
%Calculates the relative frequency of obtaining a non-zero component weight for spcovr.
%
%INPUT: Input parameters for sparse covariates regression
%       L_int: interval of lasso penalty tuning values (descending)
%       nsamples: number of resamples
%       splitratio: ratio of sample size resampled data to sample size
%       original data
%
%OUTPUT: PrMAT: matrix with relative frequency of non-zero component weight
%
%K. Van Deun, OCT2015; checked on 23 MARCH 2017
[I Jx]=size(X);
sample_size=round(splitratio*I);

PrMAT=zeros(Jx,R);
toosparse=0;
%reference for determ. comp.scores
[W,Px,Py,Loss,RsqX,Rsqy]=spcovr(X,Y,R,alpha,LASSO,RIDGE,MAXITER,CONV,INIT,[],asl,obl);
refmatrix=X*W;
if sum(sum(W~=0))==0
    n=nsamples;
else
    n=1;
end;
while n<nsamples
    v=randperm(I,round(splitratio*I));
    X_s=X(v,:);
    X_s=STD(X_s);
    X_s=X_s/(sqrt(size(X_s,1)-1));
    Y_s=Y(v,:);
    %Y_s=Y_s-mean(Y_s);
    [W,Px,Py,Loss,RsqX,Rsqy]=spcovr_randomized(X_s,Y_s,R,alpha,LASSO,RIDGE,MAXITER,CONV,INIT,[],asl,obl,0.5);
    if R>1
        [comporder reflex tucker tuckervector]=tuckercongruence_pr2(Px,PTRUE);
        W=W(:,comporder);
        PrMAT=PrMAT+(W~=0);
    else
        PrMAT=PrMAT+(W~=0);
    end;
    if sum(sum(abs(W)))<0.5*qlambda%do not train too long on penalties yielding too sparse results
        toosparse=toosparse+1;
        if toosparse>3 & n<5
            n=nsamples;
        end;
    end;
    n=n+1;
end;
PrMAT=PrMAT/nsamples;