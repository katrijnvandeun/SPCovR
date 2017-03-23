function [P] = rlsfast(X,T,P0);
%[P] = rls(X,T,P0);
% Function calculates P for problem (X-TP') with restriction
% diag(P'P) = 1;
%NOV 2015: Implemented by K.Van Deun: efficient variant of RLS of
%J.Westerhuis; checked on 23d MARCH 2017

[J,R] = size(P0);
X=X';
P = P0;
TT=T'*T;
for iter=1:10
    for r=1:R
        Pnotr=P;
        Pnotr(:,r)=[];
        tvec=TT(:,r);
        tvec(r)=[];
        %crosst=T'*T(:,r);
        %crosst(r)=0;
        %sumhat=sum((ones(J,1)*crosst').*P,2);
        sumhat=Pnotr*tvec;
        p = X*T(:,r) - sumhat;     
        if sum(p.^2~=0)%no division by zero
        P(:,r) = p/norm(p);
        end;
    end;
end; 