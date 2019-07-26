function Image3D = load_3D(o,r,y,x,c)
%given round,y,x index of tile and channel c, this returns the full 3D image.

Image3D = zeros(o.TileSz,o.TileSz,o.nZ);
for z = 1:o.nZ
    Image3D(:,:,z) = imread(o.TileFiles{r,y,x,c}, z);
end
end