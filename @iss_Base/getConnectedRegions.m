function o = getConnectedRegions(o)

Pos = o.TileInitialRawPosYX;
L = size(Pos,1);
Sz = o.TileSz*0.91;
D=[];
for i=1:L
    for j=1:L
       D(i,j,:) =  Pos(i,:) - Pos(j,:);
    end
end

NonOver = abs(D(:,:,1))>Sz | abs(D(:,:,2))>Sz;

Connected = bwconncomp(~NonOver,4);
C = {};
for i=1:length(Connected.PixelIdxList)
    [row col] = ind2sub(size(NonOver),Connected.PixelIdxList{i});
    if any(~ismember(row,col))
        C{i} = [];
    else
        C{i} = unique([row;col]);
    end
end
C = C(cell2mat(cellfun(@(x) ~isempty(x),C,'UniformOutput',false)));
o.TileConnectedID = fliplr(C);

% plot results


