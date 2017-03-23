%***********************************
%       LOAD DATA+PRE-PROCESS
%***********************************
load ../DATA/TIVD3_rev
X=TIVBlockD3_std';
%scale to sum-of-squares to allow more efficient calculations of the SPCovR
%estimates
X=X/(sqrt(size(X,1)-1));
load ../DATA/TIVtiter;
y=log2(M);
y=y-mean(y);
%External validation: 2007 season
load ../DATA/TIVD3_2007_rev
DATA2007=TIVBlockD3_std;
DATA2007=DATA2007/(sqrt(size(DATA2007,2)-1));
load ../DATA/TIVtiter2007
y2007=log2(M);
y2007=y2007-mean(y2007);
%annotation
fid=fopen('../DATA/annotation.txt');
AFFYID_all=textscan(fid,'%s%s%s');
AFFYID=AFFYID_all{1};
GENEID=AFFYID_all{2};
ID=AFFYID_all{3};
fclose(fid);

%SELECT 54 675 shared probesets
X=X(:,1:size(DATA2007,1));

%***********************************
%       SET PARAMETERS SPCOVR
%***********************************
R=2;%Number of components
MAXITER=500;%Maximal number of iterations for SPCovR procedure 
LASSOrel=0.95;%proportional weight of lasso vs ridge
NLAMBDA=500;%number of values for tuning parameters
ALPHA=.99;%value for the tuning parameter as in the paper: \alphaR²x + (1-\alpha)R²y
CONV=1e-6;%convergence criterion defined as difference in loss
INIT='rational';%starting configuration: 'random' or 'rational'
ASL=0;%active set learning
OBL=0;%oblique loadings with OBL=1 and orthogonal with OBL=0

%***********************************
%     SET PARAMETERS SUBSAMPLING
%***********************************
%Based on Meinshausen & Buhlmann (2010)
nsamples=250;%number of resamples
splitratio=0.9;%proportional sample size for resamples
EV=1;%error rate: expected number of false positives
pi_thr=0.9;%probability trhreshold


%set state of random generator
c=clock;
rng(randi(1e6,1)+c(6))



%******END OF USER INPUT*********


ALPHA=1-ALPHA;%The (s)pcovr procedure is implemented with \alphaR²y + (1-\alpha)R²x

%******************************************
%          STABILITY SELECTION PROCEDURE
%******************************************

%1. Calculate suitable range of lasso tuning parameter values (See lasso function of MATLAB)
[L L_inta]=maxLambda(X,y,R,ALPHA,NLAMBDA,ASL,OBL,LASSOrel);
%in trial runs, these values for the lasso yielded low selection
%probabilities, also taking L_inta(490). Hence refinement:
L=L_inta(400);
Lmin=(1e-1)*L;
nLambda=100;
step=(log(L)-log(Lmin))/(nLambda-1);
L_int=exp([log(L):-step:log(Lmin)]);

%Calculate number of non-zero loadings given E(V), probability threshold
%and number of variables and components
Jx=size(X,2);
qlambda=R*sqrt(EV*Jx*(2*pi_thr-1));
q=0;
teller=1;

%Initialize zero/non-zero status loading matrix with all zero values
SEL=zeros(Jx,R);

%*****
%This code may be useful to restart somewhere in loop over lasso values or
%to recalculate the final zero/non-zero status of the loadings
% for i=1:38
%     eval(['load PrMAT_randomN250_orderindDEBUG',num2str(i),'.mat'])
%     SEL=SEL+(PrMAT>=pi_thr);
%     q=sum(sum(SEL~=0));
%     teller=teller+1;
% end;
% SEL_old=SEL;
% RIDGE=L_int(38)*(1-LASSOrel);
%*******

%loop over lambda values, from largest to smallest, until the number of
%non-zero coefficients is >= q
while q<qlambda
    SEL_old=SEL;%we need the matrix BEFORE the while loop is interrupted, 
    %this is SEL_old (in order to maintain the expected number of false
    %positives conservatively under control)
    %2. Calculate probability matrices for each of the lasso tuning
    %parameter values
    LASSO=LASSOrel*L_int(teller);
    RIDGE=L_int(teller)*(1-LASSOrel);
    [PrMAT]=subsampling(X,y,R,ALPHA,LASSO,RIDGE,MAXITER,CONV,INIT,nsamples,splitratio,ASL,OBL,teller);
    SEL=SEL+(PrMAT>=pi_thr);
    q=sum(sum(SEL~=0))
    teller=teller+1;
end;


%***************************************
%           SOLUTION: POST-PROCESSING
%***************************************

%use output stability selection to find constrained spcovr solution and
%take higher value for ridge given the very low sample sizes for these data
LASSO=0;
k=find(SEL_old<1);
SEL_old(k)=0;
[W,Px,Py,Loss,RsqX,Rsqy]=spcovr(X,y,R,ALPHA,LASSO,10*RIDGE,MAXITER,CONV,INIT,SEL_old,ASL,OBL);
T=X*W;
yhat=T*Py';
%VAF for: reported in the paper in Table 1
RsqX 
Rsqy

%2007 data
T2007=DATA2007'*W(1:54675,:);
yhat2007=T2007*Py';
Rsqy2007=(corr(yhat2007,y2007))^2
Expvary2007=1-(sum((y2007-yhat2007).^2)/sum(y2007.^2))

%check gene names of selected genes and save their affymetrix probeset
%identifiers
for r=1:R
    GENEID(W(:,r)~=0)
    affies=AFFYID(W(:,r)~=0)
    ['save ../DATA/affyidentifiers_SPCOVR_pc',num2str(r),'.txt affies']
end;