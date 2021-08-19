function [SpotColors, PointCorrectedLocalYX] = get_spot_colors_grid(o, pf,...
    xy, ImSz, SpotNo, SpotLocation, IncludeGT, Filter)

if nargin<7 || isempty(IncludeGT)
    IncludeGT = false;
end
if nargin<8 || isempty(Filter)
    Filter = true;
end

fprintf('loading channel/round images...');
if SpotLocation == false
    %Find tile that the point is on and local centered coordinates in reference round
    t = o.get_local_tile([xy(2),xy(1)]);
else
    t = o.([pf,'LocalTile'])(SpotNo);
end
LocalYX = [xy(2),xy(1)]-o.TileOrigin(t,:,o.ReferenceRound);
if IncludeGT
    o.UseRounds = [o.UseRounds,o.gtRounds];
end
nRounds = max(o.UseRounds);
[RoundTile,~] = get_SpotTileEachRound(o,[xy(2),xy(1)],t);
load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
[SpotColor,PointCorrectedLocalYX] = get_spot_colors(o,LocalYX,t,...
    RoundTile,AllBaseLocalYX);
if SpotLocation==true
    if max(max(abs(double(o.([pf,'SpotColors'])(SpotNo,:,o.UseRounds(o.UseRounds<=o.nRounds)))...
            -SpotColor(:,:,o.UseRounds(o.UseRounds<=o.nRounds)))))>1e-10
        warning('Spot Color found is different from than in o object');
    end
end

if ~Filter
    t_rawdata_round = zeros(max(o.UseRounds));
    for r=o.UseRounds
        t_rawdata_round(r) = str2double(cell2mat(extractBetween(o.TileFiles{r,RoundTile(r)},'_t','.tif')));
    end
end

SpotColors = zeros((ImSz*2+1)^2,o.nBP,nRounds);
for r=1:nRounds 
    if ~Filter
        imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);  %raw data file name for round r
        % construct a Bio-Formats reader with the Memoizer wrapper
        bfreader = loci.formats.Memoizer(bfGetReader(), 0);
        bfreader.setId(imfile);
        % get some basic image metadata
        [nSeries, nSerieswPos, ~, nZstacks, ~, ~] = get_ome_tilepos(bfreader);
        scene = nSeries/nSerieswPos;
        bfreader.close();
        
        % a new reader per worker
        bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
        % use the memo file cached before
        bfreader.setId(imfile);
        bfreader.setSeries(scene*t_rawdata_round(r)-1);
    end
    for b=1:o.nBP
        rbYX = round(PointCorrectedLocalYX(1,:,r,b));
        y0 = rbYX(1);
        x0 = rbYX(2);
        if y0>o.TileSz || y0<1 || x0>o.TileSz || x0<1
            continue;
        end
        y1 = max(1,y0 - ImSz);
        y2 = min(o.TileSz,y0 + ImSz);
        x1 = max(1,x0 - ImSz);
        x2 = min(o.TileSz,x0 + ImSz);
        YBaseIm_Ind = 1:y2-y1+1;
        XBaseIm_Ind = 1:x2-x1+1;
        %If overlapping with tile to left or below, then set pixels to the
        %left or below as nan.
        BaseIm = nan(ImSz*2+1,ImSz*2+1);
        if x1==1
            XBaseIm_Ind = ImSz*2+1-x2+1:ImSz*2+1;
        end
        if y1==1
            YBaseIm_Ind = ImSz*2+1-x2+1:ImSz*2+1;
        end
        if Filter
            BaseIm(YBaseIm_Ind,XBaseIm_Ind) = int32(imread(o.TileFiles{r,t}, b, 'PixelRegion',...
                {[y1 y2], [x1 x2]}))-o.TilePixelValueShift;
            if o.SmoothSize
                SE = fspecial3('ellipsoid',o.SmoothSize);
                BaseImSm = imfilter(BaseIm, SE);
            else
                BaseImSm = BaseIm;
            end
        else
            I = cell(nZstacks,1);
            for z = 1:nZstacks
                iPlane = bfreader.getIndex(z-1, b-1, 0)+1;
                I{z} = bfGetPlane(bfreader, iPlane, x1, y1, 2*ImSz+1, 2*ImSz+1);
                %I{z} = bfGetPlane(bfreader, iPlane);
            end
            % focus stacking
            BaseIm(YBaseIm_Ind,XBaseIm_Ind) = o.fstack_modified(I(o.FirstZPlane:end));
            BaseImSm = BaseIm;
            %BaseImSm = BaseImSm(y1:y2,x1:x2);
        end
        SpotColors(:,b,r) = BaseImSm(:);
    end
end
fprintf('done\n');
end

