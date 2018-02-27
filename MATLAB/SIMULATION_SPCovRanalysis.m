%Script for the analysis of the simulated data

%***********************************
%       SET FIXED PARAMETERS SPCOVR
%***********************************
R=2;%Number of components
MAXITER=500;%Maximal number of iterations for SPCovR procedure
LASSOrel=0.99;%proportional weight of lasso vs ridge
NLAMBDA=500;%number of values for tuning parameters
%ALPHA=.99;%value for the tuning parameter as in the paper: \alphaR²x + (1-\alpha)R²y
CONV=1e-6;%convergence criterion defined as difference in loss
INIT='rational';%starting configuration: 'random' or 'rational'
ASL=0;%active set learning
OBL=0;%oblique loadings with OBL=1 and orthogonal with OBL=0

%***********************************
%     SET PARAMETERS SUBSAMPLING
%***********************************
%Based on Meinshausen & Buhlmann (2010)
nsamples=100;%number of resamples
splitratio=0.5;%proportional sample size for resamples
EV=1;%error rate: expected number of false positives
pi_thr=0.9;%probability trhreshold

ALPHAv=[0.01 0.50 0.99];
%load PERFORMANCE_a.txt
%PERFORMANCE=PERFORMANCE_a; %to retake the calculations from where the
%process was interrupted
PERFORMANCE=[];
for tellersim=(1):(3*3*3*20)
    for a=1:1%length(ALPHAv)  %change according to desired alpha; we ran three scrips in parallel to speed up calculations
        
        %***********************************
        %       LOAD DATA+PRE-PROCESS
        %***********************************
        X=dlmread(['SIMDATA/X',num2str(tellersim),'.txt'],'\t');
        Xout=dlmread(['SIMDATA/Xout',num2str(tellersim),'.txt'],'\t');
        [I,Jx]=size(X);
        %center and scale to sum of squares equal to one per variable
        X=STDss1(X);
        Xout=STDss1(Xout);
        y=dlmread(['SIMDATA/Y',num2str(tellersim),'.txt'],'\t');
        y=STDss1(y);
        yout=dlmread(['SIMDATA/Yout',num2str(tellersim),'.txt'],'\t');
        yout=STDss1(yout);
        WTRUE=dlmread(['SIMDATA/WTRUE',num2str(tellersim),'.txt'],'\t');
        PTRUE=dlmread(['SIMDATA/PTRUE',num2str(tellersim),'.txt'],'\t');
        TTRUE=dlmread(['SIMDATA/TTRUE',num2str(tellersim),'.txt'],'\t');
        
        ALPHA=ALPHAv(a);
        
        %******************************************
        %          STABILITY SELECTION PROCEDURE
        %******************************************
        
        %1. Calculate suitable range of lasso tuning parameter values (See lasso function of MATLAB)
        L=2*R*sqrt(I*Jx)/sqrt(sum(sum(WTRUE~=0)));
        Lmin=(1e-6)*L;
        nLambda=100;
        step=(log(L)-log(Lmin))/(nLambda-1);
        L_int=exp([log(L):-step:log(Lmin)]);
        %Initialize zero/non-zero status loading matrix with all zero values
        SEL=zeros(Jx,R);
        
        %loop over lambda values, from largest to smallest, until the number of
        %non-zero coefficients is >= q
        teller=1;
        Pr1=[];
        Pr2=[];
        qlambda=sum(sum(WTRUE~=0));
        q=0;
        STEP=0;
        while (q<0.8*qlambda & STEP<3)
            SEL_old=SEL;%we need the matrix BEFORE the while loop is interrupted,
            %this is SEL_old (in order to maintain the expected number of false
            %positives conservatively under control)
            %2. Calculate probability matrices for each of the lasso tuning
            %parameter values
            LASSO=LASSOrel*L_int(teller);
            RIDGE=L_int(teller)*(1-LASSOrel);
            [PrMAT]=subsampling_knownorder(X,y,R,ALPHA,LASSO,RIDGE,MAXITER,CONV,INIT,nsamples,splitratio,ASL,OBL,qlambda,PTRUE);
            Pr1=[Pr1 PrMAT(:,1)];
            Pr2=[Pr2 PrMAT(:,2)];
            SEL=SEL+(PrMAT>=pi_thr);
            q=sum(sum(SEL~=0))
            if q>(qlambda*1.2)%control q s.t. not too far from qlambda
                STEP=STEP+1;
                q=q-sum(sum(SEL~=0));
                SEL=SEL_old;
                Pr1(:,end)=[];
                Pr2(:,end)=[];
                if oldteller==0
                    L=2*L_int(teller);
                else
                    L=L_int(oldteller);
                end;
                if oldteller<5
                    Lmin=L_int(oldteller+1);
                else
                    Lmin=L_int(teller);%(1e-4)*L;
                end;
                step=(log(L)-log(Lmin))/(nLambda-1);
                L_int=exp([log(L):-step:log(Lmin)]);
                teller=0;
            end;
            if q<0.80*qlambda
                oldteller=teller;
                teller=teller+1+round(teller/10);%speed up
            else
                oldteller=teller;
                teller=teller+1;
            end;
            if teller>nLambda
                STEP=STEP+1;
                L=LASSO/LASSOrel;
                Lmin=(1e-4)*L;
                step=(log(L)-log(Lmin))/(nLambda-1);
                L_int=exp([log(L):-step:log(Lmin)]);
                teller=1;
            end;
        end;
        [x i1]=sort(max(Pr1'),'descend');
        SEL2=zeros(Jx,R);
        SEL2(i1(1:sum(WTRUE(:,1)~=0)),1)=1;
        [x i2]=sort(max(Pr2'),'descend');
        SEL2(i2(1:sum(WTRUE(:,2)~=0)),2)=1;
        SEL=SEL2;
        %***********************************
        %       CALCULATION PERFORMANCE
        %***********************************
        %Recovery W
        %1. false positives
        k=min(sum(sum(WTRUE==0 & SEL~=0)),sum(sum(WTRUE==0 & SEL(:,[2 1])~=0)));
        nrfalsepos=k;
        %2. false negatives
        k=min(sum(sum(WTRUE~=0 & SEL==0)),sum(sum(WTRUE~=0 & SEL(:,[2 1])==0)));
        nrfalseneg=k;
        %3. Tucker congruence
        %3.1. Derive non-zero wjr and regr. weights
        [W,Px,Py,Loss,RsqX,Rsqy]=spcovr(X,y,R,ALPHA,0,RIDGE,MAXITER,CONV,INIT,SEL,ASL,OBL);
        [perm reflex tucker tuckervector]=tuckercongruence_pr2(TTRUE,X*W);
        [perm reflex tuckerW tuckervectorW]=tuckercongruence_pr2(WTRUE,W);
        [perm reflex tuckerP tuckervectorP]=tuckercongruence_pr2(PTRUE,Px);
        %Prediction error
        %2. Predict y and calculate error
        T=Xout*W;
        yhat=Xout*W*Py';
        dev=yhat-yout;
        prederr=sum(dev.^2)/sum(yout.^2)
        
        %store performance measures
        PERFORMANCE=[PERFORMANCE;nrfalsepos nrfalseneg tucker prederr ALPHA RsqX Rsqy tuckerP tuckerW]
        dlmwrite(['SIMDATA/PERFORMANCE_a.txt'],PERFORMANCE,'delimiter','\t')%performance_b and performance_c for the other alpha values
    end;
end;