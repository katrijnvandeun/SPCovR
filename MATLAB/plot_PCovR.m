%REQUIRES fig.m and exportfig.m to create eps figures
%http://www.mathworks.com/matlabcentral/fileexchange/30736

%**************************************
%       LOAD DATA+PRE-PROCESS TITERS
%**************************************
load ../DATA/TIVD3_rev
X=TIVBlockD3_std';
X=X(:,1:54675);%select probe sets also available for test data
load ../DATA/TIVtiter;
y=log2(M);
y=y-mean(y);
%External validation: 2007 season
load ../DATA/TIVD3_2007_rev
DATA2007=TIVBlockD3_std;
load ../DATA/TIVtiter2007
y2007=log2(M);
y2007=y2007-mean(y2007);%


%***********************************
%       SET PARAMETERS PCOVR
%***********************************
R=2;
MAXITER=250;
ALPHA=[.99];
ALPHAv=[0.01:.01:.99];

%*********************************
%      SCREE PLOT
%*********************************

%Create Figure 1 in the paper
[~, s, ~]=svds(X,10);
ssq=sum(sum(X.^2));
s=diag(s);
propve=s.^2/ssq;
bar(propve)
xlabel('Component')
ylabel('Variance Accounted For')
h=gcf()
%exportfig(h,'PlotPVE.EPS','Format','eps','Resolution',300)
pve=sum(s(1:2).^2/ssq)

%*********************************
%      SOLUTION OVER ALPHAv
%*********************************

RX=[];
RY=[];
RY2007=[];
for a=1:length(ALPHAv)
    alpha=ALPHAv(a);
    [W,Px,Py,RsqX,Rsqy]=pcovr(X,y,R,alpha);
    RX(a)=RsqX;
    RY(a)=Rsqy;
    T2007=DATA2007'*W;
    yhat2007=T2007*Py';
    RY2007(a)=(corr(yhat2007,y2007))^2;
end;

%Create Figure 2 in the paper: 2 lines are commented out because these
%require to download the fig.m and exportfig.m functions
%fig('units','centimeters','width',15,'height',15,'font','Arial')
plot(1-ALPHAv,RX,'b:','LineWidth',2)
axis([-0.05 1.1 0 1.05])
hold on
plot(1-ALPHAv,RY,'r:','LineWidth',2)
plot(1-ALPHAv,RY2007,'k:','LineWidth',2)
RXpls=0.1160108;%obtained with Script_sgcca_spls.R
RYpls=1;
RY2007pls=0.5465989;
plot(0,RXpls(1),'b.','MarkerSize',16)
plot(0,RYpls(1),'r.','MarkerSize',16)
plot(0,RY2007pls(1),'k.','MarkerSize',16)
h=gcf()
text([0 0 0]-0.035,[RXpls RYpls RY2007pls]-0.025,['PLS';'PLS';'PLS'])
text([0.5 0.5 0.5],[RX(50) RY(50) RY2007(50)]+0.025,{'VAF X';'R²y';'R²y2007'})
xlabel('\alpha')
ylabel('Fit estimated-observed')
%exportfig(h,'../../PlotPCOVR.EPS','Format','eps','Resolution',300)