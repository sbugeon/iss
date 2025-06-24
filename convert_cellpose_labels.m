% Convert the .png output of CellPose to cell boundary coordinates.

function DapiBound = convert_cellpose_labels(CellMap)

%% convert the mask to boundaries
NB = max(max(CellMap));
DapiBound=cell(NB,1);
[S, iSorted] = sort(CellMap(:));
markers = [ find(diff(S)); numel(CellMap)];
for h = 1:numel(markers)-1
    idces = iSorted(markers(h)+1:markers(h+1));
    if length(idces)>10
        idces = idces(round(linspace(1,length(idces),30)));
    end
    [row,col] = ind2sub(size(CellMap),idces);
    B = boundary(row ,col,0.01);
    DapiBound{h,1} = [col(B),row(B)];
end
