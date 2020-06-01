function val = norm_calc(cell_exp,ref_cell)
val_vect=cellfun(@(x,y) norm(x(~isnan(x)&~isnan(y))-y(~isnan(x)&~isnan(y))), ...
    cell_exp,ref_cell);
val = norm(val_vect);
end