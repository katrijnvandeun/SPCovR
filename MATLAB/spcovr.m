function [W,Px,Py,Loss,RsqX,RsqY]=spcovr(X,Y,R,alpha,LASSO,RIDGE,MAXITER,CONV,initial,WINIT,asl,obl)
%SPCOVR: sparse principal covariates regression
%Objective function: min_{(W,Px,Py)}
%   alpha||y-XWPy'||²/||y||²+(1-alpha)||X-XWPx'||²/||X||²+l1|W|_1+l2|Wk|_2^2
%           such that diag(Px'Px)=1,
%with alpha a weight parameter (0<alpha<1) and l1 and l2 tuning parameters
%for the lasso and ridge penalties respectively.
%
%INPUT:
%     X: predictor data
%     Y: criterion data
%     R: number of components to be extracted
%     alpha: value of the weighing parameter
%     LASSO: value of the tuning parameter for the Lasso penalty
%     RIDGE: value of the tuning parameter for the Ridge penalty
%     MAXITER: maximum number of iterations
%     conv: iterations are stopped if diff loss is smaller than conv
%     initial: starting configuration 'rational' or 'random'
%     WINIT: user-constrained zero-weights
%     asl: active set learning yes (asl=TRUE) or no (asl=FALSE)
%     obl: oblique loadings yes (obl=1) or no (obl=0)
%OUTPUT
%     W: Component weights; XW yields the component scores
%     Px: Component loadings for reconstruction of the predictor data
%     Py: Regression weights for prediction by the component scores
%     Loss: Value of the objective function
%     RsqX: proportion of variance accounted for in X block
%     RsqY: proportion of variance accounted for in Y block
%
%Author: Katrijn Van Deun
%version 0: 02OKT2015; last revision 21MARCH2017

[I,Jx]=size(X);
[I,Jy]=size(Y);
%eps=1e-12;

%include warning about assumed scaling of X!
if sum(abs(sum(X.^2)-ones(1,Jx)))> 1e-6
    fprintf('This implementation supposes X to be scaled to sum of squares equal to one!\n')
    fprintf('Please use the spcovr_nonstandardized function!\n')
    return;
end;

%weights s.t. data are normalized
J=Jx+Jy;
w1=I*J*alpha/(sum(sum(Y.^2)));
w2=I*J*(1-alpha)/(sum(sum(X.^2)));

%concatenated data
wy=sqrt(w1)*Y;
wX=sqrt(w2)*X;
Xstar=[wy wX];

% initialize W & T: random or rational
switch initial
    case 'random'
        W=randn(Jx,R);
    case 'rational'
        [W,Px,Py,RsqX,Rsqy]=pcovr(X,Y,R,alpha);
end;
if prod(size(WINIT))~=0 | asl
    initw=1;
    W(WINIT==0)=0;
    WINIT=W;
else
    initw=0;
end;
T=X*W;
[U S Px]=svds(wX,R);
Py=randn(Jy,R);
P=[Py;Px];

% initial loss
pen1=LASSO*sum(sum(abs(W)));
sqW=W.^2;
pen2=RIDGE*(sum(sum(sqW)));
residual=sum(sum((Xstar-T*P').^2));
ssxstar=sum(sum((Xstar).^2));
Lossc=(residual+pen1+pen2)/ssxstar;
%fprintf('Init %5.0f Loss: %8.4f \n',0,Lossc)

%upfront calculation of matrices / vectors that remain constant in the
%iterative procedure
%sumXsq=sum(X.^2);
XX=wX*wX';

conv=0;
iter=1;
if initw~=1%unconstrained case
    while conv==0
        
        if obl
            Px=UpdatePx_oblique(wX,T,Px);
        else
            Px=UpdatePx_orthogonal(wX,XX,T);%slow for n>p!
        end;
        Py=pinv(T)*wy;%LU decomposition?
        P=[Py';Px];
        
        residual=sum(sum((Xstar-T*P').^2));
        Lossu=(residual+pen1+pen2)/ssxstar;
        %        fprintf('P: \t %5.0f Loss: %8.4f Diff: %8.4f \n',iter,Lossu,Lossc-Lossu)
        if Lossc-Lossu<0
            Lossc-Lossu
        end;
        
        %update W: soft thresholding
        sumPsq=sum(P.^2);
        %        W=updateW(Xstar,X,sumXsq,LASSO,RIDGE,W,P,sumPsq,R,Jx);
        W=updateW(Xstar,X,LASSO,RIDGE,W,P,sumPsq,R,Jx);
        
        %check all zero
        if prod(sum(W))==0
            if sum(W==0)==Jx
                fprintf('Strong penalty: All component weights are equal to zero \n')
            else
                fprintf('Al least one component has all weights equal to zero \n')
            end;
            Loss=Lossu;
            Py=P(1,:)/sqrt(w1);
            Px=P(2:end,:)/sqrt(w2);
            RsqX=1-(sum(sum((X-T*Px').^2))/sum(sum(X.^2)));
            RsqY=1-(sum(sum((Y-T*Py').^2))/sum(sum(Y.^2)));
            return;
        end;
        T=X*W;
        pen1=LASSO*sum(sum(abs(W)));
        sqW=W.^2;
        pen2=RIDGE*(sum(sum(sqW)));
        residual=sum(sum((Xstar-T*P').^2));
        Lossu2=(residual+pen1+pen2)/ssxstar;
        
        %      fprintf('W: \t %5.0f Loss: %8.4f Diff: %8.4f \n',iter,Lossu2,Lossu-Lossu2)
        if Lossu-Lossu2<0
            Lossu-Lossu2;
        end;
        
        %account for divergence: one/few very small w and very high py
        if max(max(abs(T)))<1e-9
            conv=1;
            Loss=Lossu;
            Py=P(1:Jy,:)/sqrt(w1);
            Px=P((Jy+1):end,:)/sqrt(w2);
            RsqX=1-(sum(sum((X-T*Px').^2))/sum(sum(X.^2)));
            RsqY=1-(sum(sum((Y-T*Py').^2))/sum(sum(Y.^2)));
            return;
        end;
        
        %       fprintf('W&P:  %4.0f Loss: %8.4f Diff: %8.4f \n',iter,Lossu2,Lossc-Lossu2)
        
        %check convergence
        if iter>1 && (Lossc-Lossu2)<CONV
            Loss=Lossu;
            Py=P(1:Jy,:)/sqrt(w1);
            Px=P((Jy+1):end,:)/sqrt(w2);
            RsqX=1-(sum(sum((X-T*Px').^2))/sum(sum(X.^2)));
            RsqY=1-(sum(sum((Y-T*Py').^2))/sum(sum(Y.^2)));
            return;
        elseif iter==MAXITER
            Loss=Lossu;
            Py=P(1:Jy,:)/sqrt(w1);
            Px=P((Jy+1):end,:)/sqrt(w2);
            RsqX=1-(sum(sum((X-T*Px').^2))/sum(sum(X.^2)));
            RsqY=1-(sum(sum((Y-T*Py').^2))/sum(sum(Y.^2)));
            return;
        end;
        iter=iter+1;
        Lossc=Lossu2;
    end;
else%constrained case or active set learning to speed up computations
    while conv==0
        
        if obl
            Px=UpdatePx_oblique(wX,T,Px);
        else
            Px=UpdatePx_orthogonal(wX,XX,T);%slow for n>p!
        end;
        Py=pinv(T)*wy;%LU decomposition?
        P=[Py';Px];
        
        residual=sum(sum((Xstar-T*P').^2));
        Lossu=(residual+pen1+pen2)/ssxstar;
        %      fprintf('P: \t %5.0f Loss: %8.4f Diff: %8.4f \n',iter,Lossu,Lossc-Lossu)
        if Lossc-Lossu<0
            fprintf('Increase in loss function value after updating the loadings!!! \n')
            Lossc-Lossu
        end;
        
        %update W: soft thresholding
        sumPsq=sum(P.^2);
        W=updateW_constrained(Xstar,X,LASSO,RIDGE,W,P,sumPsq,R,Jx,WINIT);
        %        W=updateW_constrained(Xstar,X,sumXsq,LASSO,RIDGE,W,P,sumPsq,R,Jx,WINIT);
        %active set learning here means that a weight that becomes zero is
        %no longer updated; asl initialized after 10th iteration
        if iter>10 && asl
            WINIT(W==0)=0;
        end;
        
        %check all zero
        if prod(sum(W))==0
            if sum(W==0)==Jx
                fprintf('Strong penalty: All component weights are equal to zero \n')
            else
                fprintf('Al least one component has all weights equal to zero \n')
            end;
            Loss=Lossu;
            Py=P(1,:)/sqrt(w1);
            Px=P(2:end,:)/sqrt(w2);
            RsqX=1-(sum(sum((X-T*Px').^2))/sum(sum(X.^2)));
            RsqY=1-(sum(sum((Y-T*Py').^2))/sum(sum(Y.^2)));
            return;
        end;
        T=X*W;
        pen1=LASSO*sum(sum(abs(W)));
        sqW=W.^2;
        pen2=RIDGE*(sum(sum(sqW)));
        residual=sum(sum((Xstar-T*P').^2));
        Lossu2=(residual+pen1+pen2)/ssxstar;
        
        %     fprintf('W: \t %5.0f Loss: %8.4f Diff: %8.4f \n',iter,Lossu2,Lossu-Lossu2)
        if Lossu-Lossu2<0
            fprintf('Increase in loss function value after updating the weights!!! \n')
            Lossu=Lossu2
        end;
        
        %account for divergence: one/few very small w and very high py
        if max(max(abs(T)))<1e-9
            conv=1;
            Loss=Lossu;
            Py=P(1:Jy,:)/sqrt(w1);
            Px=P((Jy+1):end,:)/sqrt(w2);
            RsqX=1-(sum(sum((X-T*Px').^2))/sum(sum(X.^2)));
            RsqY=1-(sum(sum((Y-T*Py').^2))/sum(sum(Y.^2)));
            return;
        end;
        
        %    fprintf('W&P:  %4.0f Loss: %8.4f Diff: %8.4f \n',iter,Lossu2,Lossc-Lossu2)
        
        %check convergence
        if iter>1 & (Lossc-Lossu2)<CONV
            Loss=Lossu;
            Py=P(1:Jy,:)/sqrt(w1);
            Px=P((Jy+1):end,:)/sqrt(w2);
            RsqX=1-(sum(sum((X-T*Px').^2))/sum(sum(X.^2)));
            RsqY=1-(sum(sum((Y-T*Py').^2))/sum(sum(Y.^2)));
            return;
        elseif iter==MAXITER
            Loss=Lossu;
            Py=P(1:Jy,:)/sqrt(w1);
            Px=P((Jy+1):end,:)/sqrt(w2);
            RsqX=1-(sum(sum((X-T*Px').^2))/sum(sum(X.^2)));
            RsqY=1-(sum(sum((Y-T*Py').^2))/sum(sum(Y.^2)));
            return;
        end;
        iter=iter+1;
        Lossc=Lossu2;
    end;
end;

%function W=updateW(Xstar,X,sumXsq,LASSO,RIDGE,W,P,sumPsq,R,Jx)
function W=updateW(Xstar,X,LASSO,RIDGE,W,P,sumPsq,R,Jx)
for r=1:R
    T=X*W;
    RES=Xstar-T*P';
    sumPRES=RES*P(:,r);
    updorder=randperm(Jx);%make coord.desc. indep. of order of variables
    for j2=1:Jx
        j=updorder(j2);%%%%%%%
        wold=W(j,r);
        %        CP=sum(sumPRES.*X(:,j))+(wold*(sumXsq(j)*sumPsq(r)));
        CP=sum(sumPRES.*X(:,j))+(wold*sumPsq(r));
        wols=sign(CP).*(abs(CP)-LASSO/2);
        w=(wols)./(RIDGE+sumPsq(r));
        %        w=(wols)./(RIDGE+(sumXsq(j)*sumPsq(r)));
        if abs(CP)<LASSO/2
            w=0;
            sumPRES=sumPRES+X(:,j)*sumPsq(r)*wold;%sumPsq ~= 1
            %           sumPRES=sumPRES+X(:,j)*sumPsq(r)*wold;%sumPsq ~= 1
        else
            sumPRES=sumPRES+X(:,j)*sumPsq(r)*(wold-w);
            %            sumPRES=sumPRES+X(:,j)*sumPsq(r)*(wold-w);
        end;
        W(j,r)=w;
    end;
end;

function W=updateW_constrained(Xstar,X,LASSO,RIDGE,W,P,sumPsq,R,Jx,WINIT)
%function W=updateW_constrained(Xstar,X,sumXsq,LASSO,RIDGE,W,P,sumPsq,R,Jx,WINIT)
for r=1:R
    T=X*W;
    RES=Xstar-T*P';
    sumPRES=RES*P(:,r);
    updorder=randperm(Jx);%make coord.desc. indep. of order of variables
    sumPsqr=sumPsq(r);
    for j2=1:Jx
        j=updorder(j2);%%%%%%%
        if WINIT(j,r)~=0
            wold=W(j,r);
            %            CP=sum(sumPRES.*X(:,j))+(wold*(sumXsq(j)*sumPsq(r)));
            CP=sum(sumPRES.*X(:,j))+(wold*sumPsqr);
            wols=sign(CP).*(abs(CP)-LASSO/2);
            %            w=(wols)./(RIDGE+(sumXsq(j)*sumPsq(r)));%*******
            w=(wols)./(RIDGE+sumPsqr);%*******
            if abs(CP)<LASSO/2
                w=0;
                %                sumPRES=sumPRES+X(:,j)*sumPsq(r)*wold;
                sumPRES=sumPRES+X(:,j)*wold*sumPsqr;
            else
                %                sumPRES=sumPRES+X(:,j)*sumPsq(r)*(wold-w);
                sumPRES=sumPRES+X(:,j)*(wold-w)*sumPsqr;
            end;
            W(j,r)=w;
        end;
    end;
end;

function Px=UpdatePx_oblique(wX,T,Px)
Px=rlsfast(wX,T,Px);

function Px=UpdatePx_orthogonal(wX,XX,T)
%update Px st orthogonal: heavy computation for large J
%use 'Kernel trick': via eig of T'wXwX'T
K1=wX'*T;
K=T'*XX*T;
[V,Ssq]=eig(K);
S=diag(((diag(Ssq)).^(-0.5)));
Px=K1*V*S*V';