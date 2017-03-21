function [W,Px,Py,RsqX,Rsqy]=pcovr(X,Y,R,alpha);
% Principal Covariates Regression: Closed form (see Heij et al., 2005)
%
% Combines regression of Y on X (ON DATA AS INPUT, USUALLY STANDARDIZED) and pca of X
% maximizes alfa Rx^2 + (1-alfa) Ry^2, or equivalently, minimizes 
% alfa||Y-XWPy'||²/||Y||²+(1-alfa)||X-XWPy'||²/||X||²  s.t. T'T is diagonal
%
% input: X, Y, r, alfa
%
% output:
%       W: component weights matrix
%       Py: regression weights matrix
%       Px: loadings
%       Rx2: proportion of explained variance in X
%       Ry2: proportion of explained variance in Y


[I,Jx]=size(X);
[Iy,Jy]=size(Y);
J=Jx+Jy;
if Iy~=I, disp(' size Y and X do not match ');return;end;
eps=1e-12;

%weights s.t. data are normalized
w1=I*J*alpha/(sum(sum(Y.^2)));
w2=I*J*(1-alpha)/(sum(sum(X.^2)));

%concatenated data
wY=sqrt(w1)*Y;
wX=sqrt(w2)*X;
Xstar=[wY wX];

%SVD data
if J>I  %account for high-dimensional data
    XX=X*X';
    [Ux,Ssq]=eig(XX);
    Sxsq=diag(Ssq);
    Sxsq(Sxsq<1e-6)=1e-6;
    Sx=sqrt(Sxsq);
    invSx=Sx.^(-1);
    %invSx(Sxsq<1e-6)=0;
    Vx=X'*Ux*(diag(invSx));
    UXstar=Ux'*Xstar;
    [U S V]=svds(UXstar,R);
elseif I>=J
    XX=X'*X;
    [Vx,Ssq]=eig(XX);
    Sxsq=diag(Ssq);
    Sxsq(Sxsq<1e-6)=1e-6;
    Sx=sqrt(Sxsq);
    invSx=Sx.^(-1);
    Ux=X*Vx*(diag(invSx));
    UXstar=Ux'*Xstar;
    [U S V]=svds(UXstar,R);
end;
P=V(:,1:R);
W=Vx*diag(invSx)*U*S;
Py=P(1:Jy,:)/sqrt(w1);
Px=P(Jy+1:end,:)/sqrt(w2);

%fit
RsqX=1-sum(sum((X-X*W*Px').^2))/(sum(sum(X.^2)));
Rsqy=1-sum(sum((Y-X*W*Py').^2))/(sum(sum(Y.^2)));