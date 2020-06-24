function lambda_A = fit_exp(A,th)
%Fit the row values above the threshold to a decaying exponential
dimA = size(A);
lambda_A = NaN(dimA(1),1);

for i = 1:dimA(1)
    currRow = A(i,:);
    if sum(~isnan(currRow)) == 0
        lambda_A(i)=NaN;
        continue
    else
        above_th = currRow(currRow > th)';
        above_th = above_th - min(above_th);
        pd = fitdist(above_th,'Exponential');
        lambda_A(i)=pd.mu;
    end
end
end