function STDDATA=STD(DATA)
%STD standardizes the DATA per column, over the rows

[n m]=size(DATA);
v = ones(n,n);
DATAc = DATA-(v*DATA/n); 
CP=v*(DATAc.^2);
STDDATA = sqrt(n-1)*DATAc./(CP.^(0.5));