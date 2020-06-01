function filledMatrix = wKNN_impute(data)
%% Weighted KNN Function

corrs =corrcoef(data,'Rows','pairwise');
distances = zeros(length(data));
for i = 1:length(data)
    v_1 = data(i,:);
    for j = 1:length(data)
        v_2= data(j,:);
        if i == j
            continue
        else
            d_ij = abs(v_1 - v_2);
            missVar = isnan(d_ij);
            temp_corrs = corrs(~missVar,~missVar);
            distances(i,j) = sqrt(d_ij(~missVar)*temp_corrs*d_ij(~missVar)');
            distances(j,i) = sqrt(d_ij(~missVar)*temp_corrs*d_ij(~missVar)');
        end
    end
end
filledMatrix = data;
[missRow, ~] = find(isnan(data));
for i = 1:length(missRow)
    [missRow, missCol] = find(isnan(data));
    varMissing = missCol(i);
    currPat = missRow(i);
    missCol(i) = 0;
    missRow(i) = 0;
    tempDistances = distances;
    tempDistances(missRow(find(missCol == varMissing)),:)=NaN;
    tempDistances(:,missRow(find(missCol == varMissing)))=NaN;
    [B,I] = mink(tempDistances(currPat,:),8);
    B(find(I == currPat))=[];
    I(find(I == currPat))=[];
    weights = 1./(B.^2);
    weights = weights./sum(weights);
    imputeValue = dot(weights,data(I,varMissing));
    filledMatrix(currPat,varMissing)=imputeValue;
end
end