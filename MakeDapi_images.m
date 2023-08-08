%% For each tile, first finds a coarse shift between the DAPI and Anchor
% Round using anchor channel

% load and store ref images for the anchor round/channel
rr = o.ReferenceRound;
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';
RefImages = zeros(o.TileSz, o.TileSz, nY, nX, 'uint16');

for t=NonemptyTiles(:)'
    [y,x] = ind2sub([nY nX], t);
    if mod(t,10)==0; fprintf('Loading tile %d anchor image\n', t); end
    Im = imread(o.TileFiles{rr,y,x}, o.AnchorChannel)  - o.TilePixelValueShift+1;
    Im = imresize(Im,[o.TileSz o.TileSz]); % in case the images are not squared
    if o.RegSmooth
        RefImages(:,:,t) = imfilter(Im, fspecial('disk', o.RegSmooth));
    else
        RefImages(:,:,t) = Im;
    end
end

% Align new reference round to old one
AnyShift = [ o.DapiRound ~= o.ReferenceRound , o.GadRound ~= o.ReferenceRound,o.GcampRound~= o.ReferenceRound];
Roundshift = [o.DapiRound,o.GadRound,o.GcampRound]; ChannelShift = [o.DapiChannel,o.GadChannel,o.GcampChannel];
RoundshiftU = unique(Roundshift(AnyShift));
WithinTileShift = nan(nTiles,2,o.nRounds+o.nExtraRounds);
o.MaxRoundShift = 100;
TileShift={};
TileCC={};
for t=NonemptyTiles
    [y,x] = ind2sub([nY nX], t);
    RefFFTStore={}; % cache to save time
    
    % match all rounds for this single tile
    for r = RoundshiftU
        
        ThisCh = find(Roundshift == r);
        Im = imread(o.TileFiles{r,y,x},o.AnchorChannel) - o.TilePixelValueShift+1;
        %         Im = imadjust(Im);
        Im = imresize(Im,[o.TileSz o.TileSz]);
        if o.RegSmooth
            MyTile = imfilter(Im, fspecial('disk', o.RegSmooth));
        else
            MyTile = Im;
        end
        
        % align to same tile in reference round, loading and cacheing fft if not
        % already there
        o.RegCorrThresh=[0.01 0.2];
        [shift, cc, ~] = ImRegFft2_SB(RefImages(:,:,t), MyTile, o.RegCorrThresh, o.RegMinSize);
        cc
        I1 = RefImages(:,:,t);I2 = MyTile;
        TileShift{r}(:,t) = shift;
        TileCC{r}(t)=cc;
    end
end
%% Then perform point cloud registration to find more accurate transformation (including rotations)
if ~isempty(TileShift) && any(any(isnan(TileShift{r}),1))
    error('\n Tile shift was not found for one or more tiles')
else
    %%
    AllI=cell(nTiles,3);
    for t=NonemptyTiles
        [y,x] = ind2sub([nY nX], t);
        RefFFTStore={}; % cache to save time
        
        % match all rounds for this single tile
        for r = RoundshiftU
            
            ThisCh = find(Roundshift == r);
            Im = imread(o.TileFiles{r,y,x},o.AnchorChannel);
            %         Im = imadjust(Im);
            Im = imresize(Im,[o.TileSz o.TileSz]);
            if o.RegSmooth
                MyTile = imfilter(Im, fspecial('disk', o.RegSmooth));
            else
                MyTile = Im;
            end
            
            I1 = RefImages(:,:,t);I2 = MyTile;
            shift = TileShift{r}(:,t);
            Pthresh = 90;Nsp=Inf;
            while Nsp>15000
                o.DetectionThresh = 'auto';
                o.DetectionThresh = prctile(double(I1(:)),Pthresh);
                x0 = o.detect_spots(I1,t,o.AnchorChannel,r);
                o.DetectionThresh = prctile(double(I2(:)),Pthresh);
                y0 = o.detect_spots(I2,t,o.AnchorChannel,r);
                Pthresh = Pthresh+0.5;
                Nsp = size(x0,1);
            end
            o.TileCentre = [0 0];
            o.PcIter = 1000;
            o.PcDist = 15;
            [D,xM] = PointCloudRegisterSB(y0, x0,o,-shift);
            
            tform = affine2d([D, transpose([0 0 1])]);
            invtform = invert(tform);
            sameAsInput = affineOutputView(size(I1),invtform,'BoundsStyle','SameAsInput');
            
            B = imwarp(I2',invtform,'OutputView',sameAsInput);
            [shift, cc, f] = ImRegFft2_SB(imadjust(I1') , imadjust(B-o.TilePixelValueShift), [0.01 0.2], o.RegMinSize);
            B2 = imtranslate(B,fliplr(shift));
            for c = 1:length(ThisCh)
                Im = imread(o.TileFiles{r,y,x},ChannelShift(ThisCh(c)));
                I2 = imresize(Im,[o.TileSz o.TileSz]);
                B = imwarp(I2',invtform,'OutputView',sameAsInput);
                B2 = imtranslate(B,fliplr(shift));
                AllI{t,ThisCh(c)} =B2';
            end
            close all
        end
        fprintf('\n-------\n');
    end
    %% Finally make new large images and save them
    o2 = load(fullfile(o.OutputDirectory,'oCall_spots_OMP.mat'));
    o.TileOrigin(:,:,o.ReferenceRound) = o2.o.TileOrigin(:,:,o.ReferenceRound);
    try
        I = imread(o.TileFiles{o.ReferenceRound,1}, o.AnchorChannel);
    catch
        NewTileDirectory = strrep(o.OutputDirectory,'output','tiles');
        o.TileFiles = cellfun(@(x) strrep(x,o.TileDirectory,NewTileDirectory),o.TileFiles,'UniformOutput',false);
    end

    AnchorOrigin = round(o.TileOrigin(:,:,o.ReferenceRound));
    MaxTileLoc = max(AnchorOrigin);
    BigAnchorIm = zeros(ceil((MaxTileLoc + o.TileSz)), 'uint16');
    BigAnchorIm_norm = double(BigAnchorIm);
    BigDapiIm = BigAnchorIm ;
    BigGcampIm = BigAnchorIm ;
    BigGadIm = BigAnchorIm ;
    BigDapiIm_norm = double(BigAnchorIm) ;
    BigGcampIm_norm = double(BigAnchorIm) ;
    BigGadIm_norm = double(BigAnchorIm) ;
    for t=NonemptyTiles
        MyOrigin = AnchorOrigin(t,:);
        
        if mod(t,10)==0; fprintf('Loading tile %d DAPI image\n', t); end
        if ~isfinite(MyOrigin(1)); continue; end
        
        LocalAnchorIm = imread(o.TileFiles{o.ReferenceRound,t}, o.AnchorChannel);
        BigAnchorIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) ...
            = BigAnchorIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) + LocalAnchorIm - o.TilePixelValueShift;
        BigAnchorIm_norm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) ...
            = BigAnchorIm_norm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) + ones(size(LocalAnchorIm));
        
        BigDapiIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) = makeBigImg(BigDapiIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)),AllI{t,1},o.DapiRound,o.DapiChannel,t,o);
        if ~isempty(AllI{t,1})
            BigDapiIm_norm(floor(MyOrigin(1))+(1:o.TileSz), ...
                floor(MyOrigin(2))+(1:o.TileSz)) = makeBigImg(BigDapiIm_norm(floor(MyOrigin(1))+(1:o.TileSz), ...
                floor(MyOrigin(2))+(1:o.TileSz)),ones(size(AllI{t,1})),o.DapiRound,o.DapiChannel,t,o);
        end
        
        BigGadIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) = makeBigImg(BigGadIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)),AllI{t,2},o.GadRound,o.GadChannel,t,o,0);
        if ~isempty(AllI{t,2})
            BigGadIm_norm(floor(MyOrigin(1))+(1:o.TileSz), ...
                floor(MyOrigin(2))+(1:o.TileSz)) = makeBigImg(BigGadIm_norm(floor(MyOrigin(1))+(1:o.TileSz), ...
                floor(MyOrigin(2))+(1:o.TileSz)),ones(size(AllI{t,2})),o.DapiRound,o.DapiChannel,t,o);
        end
        BigGcampIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) = makeBigImg(BigGcampIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)),AllI{t,3},o.GcampRound,o.GcampChannel,t,o,0);
        if ~isempty(AllI{t,2})
            BigGcampIm_norm(floor(MyOrigin(1))+(1:o.TileSz), ...
                floor(MyOrigin(2))+(1:o.TileSz)) = makeBigImg(BigGcampIm_norm(floor(MyOrigin(1))+(1:o.TileSz), ...
                floor(MyOrigin(2))+(1:o.TileSz)),ones(size(AllI{t,3})),o.DapiRound,o.DapiChannel,t,o);
        end
        
    end
    if all(cellfun(@isempty,AllI(:,1)))
        BigDapiIm_norm = BigAnchorIm_norm;
    end
    if all(cellfun(@isempty,AllI(:,2)))
        BigGadIm_norm = BigAnchorIm_norm;
    end
    if all(cellfun(@isempty,AllI(:,3)))
        BigGcampIm_norm = BigAnchorIm_norm;
    end
    BigDapiIm = BigDapiIm./(uint16(BigDapiIm_norm));
    BigAnchorIm = BigAnchorIm./(uint16(BigAnchorIm_norm));
    BigGcampIm = BigGcampIm./(uint16(BigGcampIm_norm));
    BigGadIm = BigGadIm./(uint16(BigGadIm_norm));

    imwrite(BigDapiIm, fullfile(o.OutputDirectory, 'background_image_fixed.tif'));
    imwrite(BigGcampIm, fullfile(o.OutputDirectory, 'tdTomato_image_fixed.tif'));
    imwrite(BigGadIm, fullfile(o.OutputDirectory, 'EGFP_image_fixed.tif'));
    imwrite(BigAnchorIm, fullfile(o.OutputDirectory, 'Anchor_image_fixed.tif'));

end