%Script to generate the data for the simulation study

%FIXED parameters
I=100;%nr of observations
J=200;%nr of variables
Przero=0.8;%proportion of zeroes
R=2;%nr of components
betas=[1; -0.02];%regression weights
nrreplics=20;%nr of replicate data sets
MAXITER=200;%Stopping criterion: maximum number of iterations
CONV=1e-6;%Stopping criterion: difference in loss/tolerance

%FACTORS
VAFx=[4/(2*J) 0.40 0.7];%variance accounted for by the components in covariates
VAFy=[0.80 0.5 0.02];%variance accounted for by the components in outcome
VAFr = [0.10 0.90;0.50 0.50;0.90 0.10];% perc. VAF in X accounted for per component

%GENERATION

%set seed
rng(171017)

teller=1;
RESULT=[];%to check if data are properly generated
for vx=1:length(VAFx)
    for vy=1:length(VAFy)
        for vr=1:size(VAFr,1)
            for n=1:nrreplics
                X=randn(I,J);%training data
                Xout=randn(I,J);%test data
                %center and scale data
                X=X-ones(I,1)*mean(X);
                X=X./(sqrt(I-1)*ones(I,1)*std(X));
                Xout=Xout-ones(I,1)*mean(Xout);
                Xout=Xout./(ones(I,1)*std(Xout));
                
                %obtain initial weights using pca decomposition
                [U S V]=svds(X,R);
                W=V;
                %impose desired level of sparseness
                z=randperm(J*R,round(Przero*J*R));
                W(z)=0;
                T=X*W;
                
                %rescale columns of W in accordance with desired VAFr
                for r=1:R
                    W(:,r)=sqrt(VAFr(vr,r)/sum(T(:,r).^2))*W(:,r);
                end;
                T=X*W;
                Y=T*betas;%initial outcome
                
                %derive weights and loadings from pcovr model with given
                %zero weights
                [W,P,Py,Loss,RsqX,Rsqy]=spcovr(X,Y,R,0.0001,0,0.01,MAXITER,CONV,'rational',W,0,0);
                T=X*W;
                %rescale columns of W in accordance with desired VAFr
                for r=1:R
                    W(:,r)=sqrt(VAFr(vr,r)/sum(T(:,r).^2))*W(:,r);
                end;
                %note RsqX and Rsqy both > 0 as data X are of full rank
                
                %create low-rank data and re-analyze to obtain data with
                %perfect fit to PCovR model
                XTRUE=X*W*P';
                XTRUEout=Xout*W*P';
                XTRUE=XTRUE./(sqrt(I-1)*ones(I,1)*std(XTRUE));
                XTRUEout=XTRUEout./(sqrt(I-1)*ones(I,1)*std(XTRUEout));
                YTRUE=XTRUE*W*betas;
                [W,P,Py,Loss,RsqX,Rsqy]=spcovr(XTRUE,YTRUE,R,0.0001,0,0.01,MAXITER,CONV,'rational',W,0,0);
                T=XTRUE*W;
                %rescale columns of W in accordance with desired VAFr
                for r=1:R
                    W(:,r)=sqrt(VAFr(vr,r)/sum(T(:,r).^2))*W(:,r);
                end;
                T=XTRUE*W;
                XTRUE=T*P';
                XTRUEout=XTRUEout*W*P';
                YTRUE=T*betas;
                Tout=XTRUEout*W;
                %rescale columns of Wout in accordance with desired VAFr
                Wout=W;
                for r=1:R
                    Wout(:,r)=sqrt(VAFr(vr,r)/sum(Tout(:,r).^2))*W(:,r);
                end;
                Tout=XTRUEout*Wout;
                YTRUEout=Tout*betas;
                
                %add noise with desired noise level
                ssqXtrue=sum(sum(XTRUE.^2));
                Ex=randn(I,J);
                Ex=Ex-ones(I,1)*mean(Ex);
                Exout=randn(I,J);
                Exout=Exout-ones(I,1)*mean(Exout);
                Ey=randn(I,1);
                Ey=Ey-mean(Ey);
                Eyout=randn(I,1);
                Eyout=Eyout-mean(Eyout);
                ssqEx=sum(sum(Ex.^2));
                ssqEy=sum(sum(Ey.^2));
                ssqYtrue=sum(sum(YTRUE.^2));
                
                %Rescale noise to desired level
                fx=sqrt(ssqXtrue*(1-VAFx(vx))/(VAFx(vx)*ssqEx));
                fy=sqrt(ssqYtrue*(1-VAFy(vy))/(VAFy(vy)*ssqEy));
                X=XTRUE+fx*Ex;
                Xout=XTRUEout+fx*Exout;
                Y=YTRUE+fy*Ey;
                Yout=YTRUEout+fy*Eyout;
                xhat1=T(:,1)*P(:,1)';
                ssq1=sum(sum(xhat1.^2))/(sum(sum((T*P').^2)));
                xhat2=T(:,2)*P(:,2)';
                ssq2=sum(sum(xhat2.^2))/(sum(sum((T*P').^2)));
                check1=(sum(sum(XTRUE.^2))/sum(sum(X.^2)));
                check2= (sum(sum(YTRUE.^2))/sum(sum(Y.^2)));
                
                RESULT=[RESULT;VAFx(vx) VAFy(vy) VAFr(vr,:) check1 check2 ssq1 ssq2];
                
                %store data as .txt (also to be analyzed with R packages)
                dlmwrite(['SIMDATA/X',num2str(teller),'.txt'],X,'delimiter','\t')
                dlmwrite(['SIMDATA/Y',num2str(teller),'.txt'],Y,'delimiter','\t')
                dlmwrite(['SIMDATA/Xout',num2str(teller),'.txt'],Xout,'delimiter','\t')
                dlmwrite(['SIMDATA/Yout',num2str(teller),'.txt'],Yout,'delimiter','\t')
                dlmwrite(['SIMDATA/WTRUE',num2str(teller),'.txt'],W,'delimiter','\t')
                dlmwrite(['SIMDATA/PTRUE',num2str(teller),'.txt'],P,'delimiter','\t')
                dlmwrite(['SIMDATA/TTRUE',num2str(teller),'.txt'],T,'delimiter','\t')
                teller=teller+1;
            end;
        end;
    end;
end;