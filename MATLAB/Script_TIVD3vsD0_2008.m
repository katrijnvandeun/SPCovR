load ../DATA/TIVRMA.txt;
DATA=TIVRMA;

%Here below are subj and time labels, in the same order as in the RMA
%pre-processed gene expression data TIVRMA.txt
Subject={'2';'2';'16';'16';'16';'29';'29';'29';'3';'3';'3';'32';'32';'32';'35';'35';'35';'38';'38';'38';'39';'39';'39';...
    '4';'4';'42';'42';'42';'43';'43';'43';'44';'44';'44';'46';'46';'46';'47';'47';'47';'48';'48';'48';'51';'51';'53';'53';...
    '53';'63';'63';'63';'65';'65';'65';'68';'68';'68';'70';'70';'72';'72';'72';'73';'73';'73';'74';'74';'74';'78';'78';'78';...
    '80';'80';'80';'83';'83';'83';'85';'85';'85'};

Time={'D0';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';...
    'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D0';'D3';'D7';'D0';'D3';...
    'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';...
    'D7';'D0';'D3';'D7'};

% make two blocks of data
%D3
D3=strmatch('D3',Time);
Subj_D3=Subject(D3,:);
%D0
D0=strmatch('D0',Time);
Subj_D0=Subject(D0,:);

%find subjects with data on D0 and D3
matchD0=[];
matchD3=[];
for i=1:length(D3)
    match0=strmatch(Subj_D3(i,:),Subj_D0,'exact');
    matchD0=[matchD0;D0(match0,:)];
    matchD3=[matchD3;D3(i,:)];
end;
%CHECK!
[Subject(matchD0,:) Subject(matchD3,:)] %are the subject ids the same in each row?
[Time(matchD0,:)  Time(matchD3,:)] %are the right time points selected?

%Collect the expression data for day 0 and day 3
TIVBlockD0=DATA(:,matchD0);
TIVBlockD3=DATA(:,matchD3);

%Obtain difference scores wrt baseline (D0)
TIVBlockD3=TIVBlockD3-TIVBlockD0;

%pre-processing: centering + scaling to unit SS per gene
TIVBlockD3_std=(STD(TIVBlockD3')');

%save data as a matlab file for spcovr analyses
save ../DATA/TIVD3_rev.mat TIVBlockD3_std

%for spls in R: in .txt format
save ../DATA/TIVD3_rev.txt TIVBlockD3_std -double -ascii