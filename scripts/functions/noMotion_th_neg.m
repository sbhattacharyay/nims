function nm_A = noMotion_th_neg(A,th)
%   Finds percentage of values in each row that are below threshold
dimA = size(A);
nm_A = NaN(dimA(1),1);
for i = 1:dimA(1)
    currRow = A(i,:);
    if sum(~isnan(currRow)) == 0
        nm_A(i)=NaN;
        continue
    else
        row_pi=(sum(currRow > th))/(sum(~isnan(currRow)));
        nm_A(i)=row_pi;
    end
end
end

