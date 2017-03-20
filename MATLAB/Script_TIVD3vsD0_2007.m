load ../DATA/TIVRMA_2007.txt;
DATA=TIVRMA_2007;

%Here below are subj and time labels, in the same order as in the RMA
%pre-processed gene expression data TIVRMA.txt
Subject={'12';'12';'12';'16';'16';'16';'18';'18';'18';'21';'21';'21';'23';'23';'23';'29';'29';'29';'32';'32';'32';'35';'35';'35';'39';'39';'39'};

Time={'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7'};

% %make two blocks of data
%D3
D3=strmatch('D3',Time);
Subj_D3=Subject(D3,:);
%D0
D0=strmatch('D0',Time);
Subj_D0=Subject(D0,:);

%find subjects with data on D0, D3 and D7
matchD0=[];
matchD3=[];
for i=1:length(D3)
    match0=strmatch(Subj_D3(i,:),Subj_D0,'exact');
    matchD0=[matchD0;D0(match0,:)];
    matchD3=[matchD3;D3(i,:)];
end;
%CHECK!
[Subject(matchD0,:) Subject(matchD3,:)]%are the subject ids the same in each row?
[Time(matchD0,:)  Time(matchD3,:)]%are the right time points selected?

TIVBlockD0=DATA(:,matchD0);
TIVBlockD3=DATA(:,matchD3);

%Obtain difference scores wrt baseline (D0)
TIVBlockD3=TIVBlockD3-TIVBlockD0;

%pre-processing: centering + scaling to unit SS per gene
TIVBlockD3_std=(STD(TIVBlockD3')');

%save data
save ../DATA/TIVD3_2007_rev.mat TIVBlockD3_std

%for spls in R
save ../DATA/TIVD3_2007_rev.txt TIVBlockD3_std -double -ascii