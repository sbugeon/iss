function rename_tiles(o)

TDcopy = strrep(o.TileDirectory,'\tiles','\tiles - Copie');

t_save_value = sub2ind([max(o.TilePosYX(:,1)),max(o.TilePosYX(:,2))],...
    o.TilePosYX(:,1),o.TilePosYX(:,2));

t_save_value2 = sub2ind([max(o.TileInitialPosYX(:,1)),max(o.TileInitialPosYX(:,2))],...
    o.TileInitialPosYX(:,1),o.TileInitialPosYX(:,2));

[~,S] = sort(t_save_value2);
Reorder = t_save_value(S);

% change names
D = dir(TDcopy);
   Dn = {D(:).name};
for i=1:length(Reorder)
   TileOldN = ['t',num2str(i),'.tif'];
   TileNewN = ['t',num2str(Reorder(i)),'.tif'];
   ThisT = find(endsWith(Dn,TileOldN));
   for j=1:length(ThisT)
        oldF = fullfile( TDcopy,Dn{ThisT(j)});
        newF = fullfile( o.TileDirectory,strrep(Dn{ThisT(j)},TileOldN,TileNewN));
        copyfile(oldF, newF);
   end
end