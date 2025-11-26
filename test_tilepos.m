load('H:\SBI010-S10-2025-11-20\Round_00\output\oExtract.mat')
t_save_value = sub2ind([max(o.TilePosYX(:,1)),max(o.TilePosYX(:,2))],...
    o.TilePosYX(:,1),o.TilePosYX(:,2));

t_save_value2 = sub2ind([max(o.TileInitialPosYX(:,1)),max(o.TileInitialPosYX(:,2))],...
    o.TileInitialPosYX(:,1),o.TileInitialPosYX(:,2));


dd = flipud(o.TilePosYX);
t_save_value0 = sub2ind([max(dd(:,1)),max(dd(:,2))],...
    dd(:,1),dd(:,2));

% o.EmptyTiles(:) = 1;
% o.EmptyTiles(15) = 0;
[~,S] = sort(t_save_value2);
Reorder = t_save_value(S);
o.RegInfo.Method = 'Fft';
[o, VerticalPairs, vShifts, HorizontalPairs, hShifts] = get_Fft_shifts(o);

xypos = max(xypos) - xypos;
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';
figure
for i=NonemptyTiles
    [y,x] = ind2sub([nY nX], i); 
    text(x,y,num2str(i))
end
xlim([0 7])
ylim([0 9])
set(gca,'YDir','reverse')