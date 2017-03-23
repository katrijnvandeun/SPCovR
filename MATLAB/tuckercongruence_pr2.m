function [perm reflex tucker tuckervector]=tuckercongruence_pr2(matrix1,matrix2)
%calculate Tucker's coefficient of congruence between columns but after
%accounting for permutational freedom and reflections. All possible permutations
%are checked; reflections are accounted for by taking the absolute value of
%
%INPUT
%   matrix1: reference matrix
%   matrix2: matrix to compare to the reference
%OUTPUT
%   perm: indices to permute the columns of matrix2 with 
%   reflex: vector of 1/-1 to reflect the vectors of matrix2
%   tucker: tucker coefficient of congruence for matrix2 with optimal
%   permutation and reflection
%   tuckervector: tucker's congruence per vector
%
%Code checked on March 23 2017

[m n]=size(matrix1);
INDIC_MAT=perms([1:n]);
TUCK=[];
for i=1:size(INDIC_MAT,1)
    matrix2_perm=matrix2(:,INDIC_MAT(i,:));
    teller=1;
    tuckerr=[];
    for r=1:n
        vec1=matrix1(:,r);
        vec2=matrix2_perm(:,r);
        cp=vec1'*vec2;
        var1=vec1'*vec1;
        var2=vec2'*vec2;
        if var1>0 & var2>0 %numerical problem with perfect dist.comp.
            tuckerr(teller)=(trace(cp))/(sqrt((trace(var1))*(trace(var2))));
            teller=teller+1;
        elseif var2==0
            tuckerr(teller)=0;
            teller=teller+1;
        end;
    end;
    tucker(i)=mean(abs(tuckerr));
    TUCK=[TUCK;tuckerr];
end;
k=find(tucker==max(tucker));
k=k(1);
perm=INDIC_MAT(k,:);
TUCK(k,:);
reflex=sign(TUCK(k,:));
tucker=max(tucker);
tuckervector=TUCK(k,:);