%script to obtain HAI titer for TIV data 2007
load ../DATA/TIVtiters_2007.txt
DATA=TIVtiters_2007;

%We need to select the titers corresponding to the selected gene expression
%data
Subject={'12';'12';'12';'16';'16';'16';'18';'18';'18';'21';'21';'21';'23';'23';'23';'29';'29';'29';'32';'32';'32';'35';'35';'35';'39';'39';'39'};

Time={'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7';'D0';'D3';'D7'};

% %make two blocks of data
%D3
D3=strmatch('D3',Time);
Subj_D3=Subject(D3,:);
%D0
D0=strmatch('D0',Time);
Subj_D0=Subject(D0,:);
%find subjects with data on D3 and D0
matchD3=[];
for i=1:length(D3)
    match=strmatch(Subj_D3(i,:),Subj_D0,'exact');
    matchD3=[matchD3;D3(i,:)];
end;

%derive HA index as described by Nakaya et al.
titers=DATA(:,matchD3);
m1=titers(2,:)./titers(1,:);
m2=titers(4,:)./titers(3,:);
m3=titers(6,:)./titers(5,:);
M=max([m1;m2;m3])';

save ../DATA/TIVtiter2007 M

save ../DATA/TIVtiter2007.txt M -double -ascii