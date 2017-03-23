function [PrMAT]=subsampling(X,Y,R,alpha,LASSO,RIDGE,MAXITER,CONV,INIT,nsamples,splitratio,asl,obl,q)
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
%reference for determ. comp.scores
[W,Px,Py,Loss,RsqX,Rsqy]=spcovr(X,Y,R,alpha,LASSO,RIDGE,MAXITER,CONV,INIT,[],asl,obl);
refmatrix=X*W;
for n=1:nsamples
    v=randperm(I,round(splitratio*I));
    X_s=X(v,:);
    X_s=STD(X_s);
    X_s=X_s/(sqrt(size(X_s,1)-1));
    Y_s=Y(v,:);
    [W,Px,Py,Loss,RsqX,Rsqy]=spcovr(X_s,Y_s,R,alpha,LASSO,RIDGE,MAXITER,CONV,INIT,[],asl,obl);
    if R>1
        comporder=maxcardoverlap(PrMAT,W);
        newmatrix=X_s*W;
        [perm reflex tucker tuckervector]=tuckercongruence_pr2(refmatrix(v,:),newmatrix);
        W=W(:,perm);
        PrMAT=PrMAT+(W~=0);
    else
        PrMAT=PrMAT+(W~=0);
    end;
end;
PrMAT=PrMAT/nsamples;
fname=['PrMAT_randomN250_orderind',num2str(q)];
save([fname],'PrMAT')