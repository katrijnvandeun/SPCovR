function [L L_int]=maxLambda(X,Y,R,alpha,nLambda,asl,obl,LASSOrel);
% Find maximal lambda and create suitable range of lasso values for sparse
% covariates regression.
%
% input: X, Y, r, alfa, nLambda
%
% output:
%       L: approximation of maximal lambda (yielding all zero weights)
%       L_int: interval from Lmax to Lmin with nLambda elements

[I,Jx]=size(X);
[I,Jy]=size(Y);
%weights s.t. data are normalized
J=Jx+Jy;
w1=I*J*alpha/(sum(sum(Y.^2)));
w2=I*J*(1-alpha)/(sum(sum(X.^2)));
%concatenated data
wy=sqrt(w1)*Y;
wX=sqrt(w2)*X;
Xstar=[wy wX];
sumXstar=sum(Xstar,2);
%create initial value based on only one non-zero component weight: this
%weight can be extremely large
L=max((sumXstar'*Xstar))*2;
Lmin=(1e-8)*L;%compensate for very large weight
step=(log(L)-log(Lmin))/(nLambda-1);
L_int=exp([log(L):-step:log(Lmin)]);

W=zeros(Jx,R);
teller=1;
teller2=1;%for control over excessive looping
% To avoid an interval with much too large lasso values, the initial range
% of lasso values created above is progressively changed to arrive at a refined interval
sumXsq=sum(X.^2);
while sum(sum(W~=0))==0 && teller < nLambda && teller2<101
    LASSO=L_int(teller)*LASSOrel
    RIDGE=L_int(teller)*(1-LASSOrel);
    %use spcovr analysis to obtain an indication about the size of the component
    %weights before soft thresholding; these are in the CP matrix
    [W,Px,Py,Loss,RsqX,RsqY]=spcovr(X,Y,R,alpha,LASSO,RIDGE,25,1e-4,'rational',[],asl,obl);
    CP=W.*(RIDGE+((sumXsq')*ones(1,R)))+(LASSO/2);
    sum(sum(W~=0))
    %control that initial range of lasso values is not to far of (too few
    %or too many non-zero values)
    if sum(sum(W~=0))>R*J/1000
        %check if this occurs in first run, this is initial interval not ok
        if teller==1
            L=max(max(abs(CP)))*2/LASSOrel;
            W=zeros(Jx,R);
            Lmin=(1e-4)*L;
            step=(log(L)-log(Lmin))/(nLambda-1);
            L_int=exp([log(L):-step:log(Lmin)]);
            teller=1;
            teller2=teller2+1;
        else
            L=(L_int(teller)+Lold)/2;
            Lmin=(1e-4)*L;
            step=(log(L)-log(Lmin))/(nLambda-1);
            L_int=exp([log(L):-step:log(Lmin)]);
            W=zeros(Jx,R);
            teller=1;
            teller2=teller2+1;
        end;
    elseif sum(sum(W~=0))==0
        teller=teller+round((nLambda-teller)/20);
        teller2=teller2+1;
        %         L=max(max(abs(CP)))*2/LASSOrel;
        %         Lmin=(1e-4)*L;
        %         step=(log(L)-log(Lmin))/(nLambda-1);
        %         L_int=exp([log(L):-step:log(Lmin)]);
        %    teller=teller-1;
    end;
    %teller=teller+1;
    Lold=LASSO/LASSOrel;
end;
% To start with value that yields all zeroes, take previous lasso value
L=L_int(teller);
step=(log(L)-log(Lmin))/(nLambda-1);
L_int=exp([log(L):-step:log(Lmin)]);