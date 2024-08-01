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
    C{i} = unique([row;col]);
end

o.TileConnectedID = fliplr(C);

% plot results


