%script to obtain HAI titer for TIV data
load ../DATA/TIVtiters.txt
DATA=TIVtiters;

%We have to match the titers to the gene expression data: this means
%selecting the subjects in exactly the same way
Subject={'2';'2';'16';'16';'16';'29';'29';'29';'3';'3';'3';'32';'32';'32';'35';'35';'35';'38';'38';'38';'39';'39';'39';...
    '4';'4';'42';'42';'42';'43';'43';'43';'44';'44';'44';'46';'46';'46';'47';'47';'47';'48';'48';'48';'51';'51';'53';'53';...
    '53';'63';'63';'63';'65';'65';'65';'68';'68';'68';'70';'70';'72';'72';'72';'73';'73';'73';'74';'74';'74';'78';'78';'78';...
    '80';'80';'80';'83';'83';'83';'85';'85';'85'};

Time={'D0';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';...
    'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D0';'D3';'D7';'D0';'D3';...
    'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';...
    'D7';'D0';'D3';'D7'};

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
    matchD3=[matchD3;D3(i,:)];
end;
titers=DATA(:,matchD3);%selects the right subjects

%calculate antibody titers as done by Nakaya et al. (2011)
m1=titers(2,:)./titers(1,:);
m2=titers(4,:)./titers(3,:);
m3=titers(6,:)./titers(5,:);
M=max([m1;m2;m3])';

save ../DATA/TIVtiter M

save ../DATA/TIVtiter.txt M -double -ascii
