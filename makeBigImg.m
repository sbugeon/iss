function BigI = makeBigImg(BigI,LocalI,r,c,t,o,NoOff)
if nargin<7
    NoOff = 1;
end
if NoOff
    Offset = o.TilePixelValueShift;
else
     Offset = 0;
end
    if c==1
        Offset = 0;
    end
 if isempty(LocalI)
     LocalI = imread(o.TileFiles{r,t}, c) - Offset;
 else
     LocalI = LocalI - Offset;
 end
    BigI = BigI + LocalI;
end