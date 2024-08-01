function  o = extract_and_filter_extraTiles(o,RefRounds)

t_save_value0 = sub2ind([max(o.TilePosYX(:,1)),max(o.TilePosYX(:,2))],...
        o.TilePosYX(:,1),o.TilePosYX(:,2));
    
UsedEmptyTiles = false;   %If running for only a few tiles, this will change to true.

for r = o.nRounds+o.nExtraRounds+1:length(o.FileBase)
    if r == o.AnchorRound; ExtractScale = o.ExtractScaleAnchor;
    else; ExtractScale = o.ExtractScale; end
    
    imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);
    
    % construct a Bio-Formats reader with the Memoizer wrapper
    bfreader = loci.formats.Memoizer(bfGetReader(), 0);
    
    if exist(imfile)>0
        % initiate reader
        bfreader.setId(imfile);
    else
        fprintf(['Data not currently in InputDirectory. Waiting for raw data for round ' num2str(r) '... ']);
        %Wait for file to exist
        if r==1
            MaxTime = o.MaxWaitTime1;   %Don't wait long if first round
        else
            MaxTime = o.MaxWaitTime;   %Max amount of time to wait in seconds
        end
        count = 0;
        while exist(imfile)==0
            pause(1);
            count = count+1;
            if count >= MaxTime
                fprintf('Not found\n');
                error(sprintf(strcat('No file named:\n  ',imfile,'\ncreated in allowed time')));
            end
        end
        fprintf('Data found\n');
        %Wait for file to stop loading
        OldBytes = 0;
        NewBytes = 0.00001;
        while NewBytes>OldBytes
            pause(5);
            fileinfo = dir(imfile);
            OldBytes = NewBytes;
            NewBytes = fileinfo.bytes;
        end
        % initiate reader
        bfreader.setId(imfile);
    end
    
    % get some basic image metadata
    [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize,xyposRaw] = ...
        get_ome_tilepos(bfreader);
    if isempty(xypos) || (nSeries>1 && size(xypos, 1)==1)
        if r == 1
            warning('first round xypos empty - using values from initial manual input')
            assert(~isempty(o.TileInitialPosYX), 'xypos unavailable')
            xypos = o.TileInitialPosXY;
            xyposOld = xypos;
        else
            warning('xypos empty - using values from previous round')
            xypos = xyposOld;
        end
        nSerieswPos = size(xypos,1);
    else
        xyposOld = xypos;
    end
    
    
    scene = nSeries/nSerieswPos;
    
    bfreader.close();
    
    %New filter
    if strcmpi(o.ExtractR1, 'auto')
        o.ExtractR1 = round(0.5/pixelsize);  %Gives value of 3 for pixelsize = 0.1669 of most data tested
    end
    if strcmpi(o.ExtractR2, 'auto')
        o.ExtractR2 = o.ExtractR1*2;
    end
    h = -hanning(o.ExtractR2*2+1);
    h = -h/sum(h);
    h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1) = ...
        h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1)+hanning(o.ExtractR1*2+1)/sum(hanning(o.ExtractR1*2+1));
    SE = ftrans2(h');
    SE = single(gpuArray(SE));
        o = o.get_TilePos(xypos,nSeries);
    if r == 1
        if isempty(o.TilePosYX)
            o = o.get_TilePos(xypos,nSeries);
        end
        nSerieswPos =  max(o.TilePosYX(:,1)) * max(o.TilePosYX(:,2));
        if isempty(o.AutoThresh)
            o.AutoThresh = zeros(nSerieswPos,nChannels,o.nRounds+o.nExtraRounds);
        end
        if isempty(o.HistCounts)
            o.HistValues = -o.TilePixelValueShift:1:2^16-o.TilePixelValueShift;  %Entire range of possible pixel values
            o.HistCounts = zeros(length(o.HistValues),nChannels,o.nRounds);
        end
        if isempty(o.nPixelsOutsideTiffRange)
            o.nPixelsOutsideTiffRange = zeros(nSerieswPos,nChannels,o.nRounds+o.nExtraRounds);
        end
        if isempty(o.PixelsOutsideTiffRangeExtractScale)
            o.PixelsOutsideTiffRangeExtractScale = nan(nSerieswPos,nChannels,o.nRounds+o.nExtraRounds);
        end
        
        
        %New filter
        if strcmpi(o.ExtractR1, 'auto')
            o.ExtractR1 = round(0.5/pixelsize);  %Gives value of 3 for pixelsize = 0.1669 of most data tested
        end
        if strcmpi(o.ExtractR2, 'auto')
            o.ExtractR2 = o.ExtractR1*2;
        end
        h = -hanning(o.ExtractR2*2+1);
        h = -h/sum(h);
        h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1) = ...
            h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1)+hanning(o.ExtractR1*2+1)/sum(hanning(o.ExtractR1*2+1));
        SE = ftrans2(h');
        SE = single(gpuArray(SE));
        
        MaxY = max(o.TilePosYX(:,1));
        MaxX = max(o.TilePosYX(:,2));
        if (min(size(o.EmptyTiles))==1 && min(o.EmptyTiles)~=0) || ...
                (min(size(o.EmptyTiles)==[MaxY, MaxX]) && min(o.EmptyTiles(:))==0 && ...
                max(o.EmptyTiles(:)==1))
            UsedEmptyTiles = true;
            if min(o.EmptyTiles(:))==0
                o.EmptyTiles = find(o.EmptyTiles(:)==0);
            end
            EmptyTilesOrig = o.EmptyTiles;
        end
        
    end
    
    % set up filename grid for this round
    fName = cell(nSerieswPos ,1);
    
    %Set top hat structuring elements
    if strcmpi(o.DapiR,'auto')
        o.DapiR = round(8/pixelsize);
    end
    DapiSE = strel('disk', o.DapiR);
    
    %Tile index in nd2 file different to index in o.EmptyTiles
    t_save_value = sub2ind([max(o.TilePosYX(:,1)),max(o.TilePosYX(:,2))],...
        o.TilePosYX(:,1),o.TilePosYX(:,2));
    
    %Get auto value for extract scale
    if strcmpi(o.ExtractScaleTile, 'auto')
        %Work out extract scale from middle tile by default
        if UsedEmptyTiles
            o.ExtractScaleTile = o.EmptyTiles(1);
            ExtractTileIdx = find(t_save_value == o.ExtractScaleTile);
        else
            [~, ExtractTileIdx]=ismember(round(mean(o.TilePosYX,1)),o.TilePosYX,'rows');
            if ExtractTileIdx == 0
                ExtractTileIdx = round(size(o.TilePosYX,1)/2);
            end
            o.ExtractScaleTile = t_save_value(ExtractTileIdx);
        end
    else
        ExtractTileIdx = find(t_save_value == o.ExtractScaleTile);
    end
    if (r==min(setdiff(1:o.nRounds+o.nExtraRounds,o.AnchorRound)) && strcmpi(o.ExtractScale, 'auto')) ||...
            (r==o.AnchorRound && strcmpi(o.ExtractScaleAnchor, 'auto'))
        
        o = o.get_extract_scale(nChannels,nZstacks,imfile,scene,gather(SE),DapiSE,r,ExtractTileIdx,o.ExtractScaleTile);
        if r==o.AnchorRound; ExtractScale = o.ExtractScaleAnchor;
        else; ExtractScale = o.ExtractScale; end
    end
    
    
    %parfor t = 1:nSerieswPos
    for t_index = 1:length(t_save_value)
        t = t_save_value(t_index);
        fName{t} = fullfile(o.TileDirectory, ...
            [o.FileBase{r}, '_t', num2str(t), '.tif']);
        
        if exist(fName{t}, 'file')
            fprintf('Round %d tile %d already done.\n', r, t);
            continue;
        elseif min(size(o.EmptyTiles))==1 && ~ismember(t,o.EmptyTiles)
            UsedEmptyTiles = true;
            EmptyTilesOrig = o.EmptyTiles;
            continue;
        end
        
        
        % a new reader per worker
        bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
        % use the memo file cached before
        bfreader.setId(imfile);
        
        bfreader.setSeries(scene*t_index-1);
        if strcmpi(ExtractScale, 'auto')
            if r==o.AnchorRound
                error(['Some tiles in imaging rounds already exist, but o.ExtractScaleAnchor = auto.'...
                    '\nThis will result in different scalings used for different tiles.'...
                    '\nIf tiles up to this point were obtained with a manual value of o.ExtractScaleAnchor,'...
                    '\nset o.ExtractScaleAnchor to this value and rerun.'...
                    '\nIf tiles up to this point were obtained with auto value,'...
                    'delete tile %d in round %d and rerun, tile location:\n%s.'],...
                    o.ExtractScaleTile,o.AnchorRound,fName{o.ExtractScaleTile});
            else
                r0 = min(setdiff(1:o.nRounds+o.nExtraRounds,o.AnchorRound));
                error(['Some tiles in imaging rounds already exist, but o.ExtractScale = auto.'...
                    '\nThis will result in different scalings used for different tiles.'...
                    '\nIf tiles up to this point were obtained with a manual value of o.ExtractScale,'...
                    '\nset o.ExtractScale to this value and rerun.'...
                    '\nIf tiles up to this point were obtained with auto value,'...
                    ' delete tile %d in round %d and rerun, tile location:\n%s.'],...
                    o.ExtractScaleTile,r0,fullfile(o.TileDirectory,[o.FileBase{r0},...
                    '_t', num2str(o.ExtractScaleTile), '.tif']));
            end
            
        else
            for c = 1:nChannels
                % read z stacks
                I = cell(nZstacks,1);
                for z = 1:nZstacks
                    iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
                    I{z} = bfGetPlane(bfreader, iPlane);
                end
                
                % focus stacking
                I_mod = o.fstack_modified(I(o.FirstZPlane:end));
                if o.StripHack
                    [I_mod,BadColumns] = o.StripHack_raw(I_mod);
                end
                
                % tophat
                if (c == o.DapiChannel  && r == o.DapiRound) ||...
                        ((c == o.GcampChannel || c == o.GadChannel) &&...
                        (r == o.DapiRound || r == o.AnchorRound))
                    if o.RawDAPIGad
                        IFS = I_mod;
                    else
                        IFS = imtophat(I_mod, DapiSE);
                    end
                    if o.StripHack
                        %Set faulty columns all to be 0 so peaks can't be found there.
                        IFS(:,BadColumns) = 0;
                    end
                else
                    I_mod = single(padarray(I_mod,(size(SE)-1)/2,'replicate','both'));
                    IFS = convn(gpuArray(I_mod),SE,'valid');
                    clearvars I_mod I  %Free up GPU memory
                    IFS = gather(IFS*ExtractScale);
                    if o.StripHack
                        %Set faulty columns all to be 0 so peaks can't be found there.
                        IFS(:,BadColumns) = 0;
                    end
                    if c ~= o.AnchorChannel && r == o.AnchorRound
                        IFS = uint16(IFS+o.TilePixelValueShift);
                    else
                        %Determine auto thresholds
                        o.AutoThresh(t,c,r) = median(abs(IFS(:))*o.AutoThreshMultiplier);
                        if ismember(r,1:o.nRounds)
                            %Get histogram data
                            IFS = int32(IFS);
                            o.HistCounts(:,c,r) = o.HistCounts(:,c,r)+histc(IFS(:),o.HistValues);
                            if o.StripHack
                                %Dont count BadColumns
                                ZeroIndex = find(o.HistValues==0);
                                o.HistCounts(ZeroIndex,c,r) = o.HistCounts(ZeroIndex,c,r)-length(BadColumns)*o.TileSz;
                            end
                        end
                        IFS = IFS+o.TilePixelValueShift;
                        nPixelsOutsideRange = sum(sum(IFS>uint16(inf)));
                        if nPixelsOutsideRange>o.nPixelsOutsideTiffRangeThresh
                            MaxValue = double((max(IFS(IFS>uint16(inf)))-o.TilePixelValueShift))/ExtractScale;
                            NewScaling = double(uint16(inf))/MaxValue;
                            o.nPixelsOutsideTiffRange(t,c,r) = nPixelsOutsideRange;
                            o.PixelsOutsideTiffRangeExtractScale(t,c,r) = NewScaling;
                            if r==o.AnchorRound
                                ErrorFile = fullfile(o.OutputDirectory, ['oExtract-Error_with_tile',num2str(t),'_round',num2str(r)]);
                                save(ErrorFile, 'o', '-v7.3');
                                error(['Round %d, tile %d, channel %d: %d pixels have reached limit of uint16 range.'...
                                    '\nCurrent value of o.ExtractScaleAnchor = %.4f is too high.'...
                                    ' Needs to be below %.4f.\nDelete all anchor round tiles and run again with o.ExtractScaleAnchor = %.4f.'...
                                    '\nBefore running the whole thing again, you can also look at the image directly by running'...
                                    '\nview_filtering(o,round,tile) with new value of o.ExtractScaleAnchor.'...
                                    '\nProgress up to this point saved as:\n%s.mat'],...
                                    r,t,c,nPixelsOutsideRange,o.ExtractScaleAnchor,NewScaling,0.85*NewScaling,ErrorFile);
                            else
                                ErrorFile = fullfile(o.OutputDirectory, ['oExtract-Error_with_tile',num2str(t),'_round',num2str(r)]);
                                save(ErrorFile, 'o', '-v7.3');
                                error(['Round %d, tile %d, channel %d: %d pixels have reached limit of uint16 range.'...
                                    '\nCurrent value of o.ExtractScale = %.4f is too high.'...
                                    'Needs to be below %.4f.\nDelete all tiles (excluding anchor round) and run again with o.ExtractScale = %.4f.'...
                                    '\nBefore running the whole thing again, you can also look at the image directly by running'...
                                    '\nview_filtering(o,round,tile) with new value of o.ExtractScale.'...
                                    '\nProgress up to this point saved as:\n%s.mat'],...
                                    r,t,c,nPixelsOutsideRange,o.ExtractScale,NewScaling,0.85*NewScaling,ErrorFile);
                            end
                        elseif nPixelsOutsideRange>0
                            o.nPixelsOutsideTiffRange(t,c,r) = nPixelsOutsideRange;
                            MaxValue = double((max(IFS(IFS>uint16(inf)))-o.TilePixelValueShift))/ExtractScale;
                            o.PixelsOutsideTiffRangeExtractScale(t,c,r) = double(uint16(inf))/MaxValue;
                            warning('Round %d, tile %d, channel %d: %d pixels have reached limit of uint16 range',...
                                r,t,c,nPixelsOutsideRange);
                        end
                        IFS = uint16(IFS);
                    end
                end
                
                % write stack image
                imwrite(IFS,...
                    fullfile(o.TileDirectory,...
                    [o.FileBase{r}, '_t', num2str(t), '.tif']),...
                    'tiff', 'writemode', 'append');
                fprintf('Round %d tile %d colour channel %d finished.\n', r, t, c);
            end
            % if the DAPI Round has less channels than other rounds,
            % the missing channels will be replaced by the Anchor channel
            if nChannels < length(o.bpLabels) && r == o.DapiRound
                for c = (nChannels + 1) : 7
                    imwrite(IFS,...
                        fullfile(o.TileDirectory,...
                        [o.FileBase{r}, '_t', num2str(t), '.tif']),...
                        'tiff', 'writemode', 'append');
                end
                fprintf(['\n There was less channels than codes for DAPI round: ',num2str(r),....
                    ' , replacing missing channels with Anchor'])
%                 nChannels = length(o.bpLabels);
            elseif nChannels < length(o.bpLabels) && r ~= o.DapiRound
                error(['\n There were less channels than codes for round: ',num2str(r)])
            end
        end
        bfreader.close();
        
    end
    
    % find closest tile in 1st acquisition 
    id = r - (o.nRounds+o.nExtraRounds);
    RR = RefRounds(id);
    
    % find closest tile
    for i=1:length(t_save_value)
    D = pdist([xyposRaw(i,:);o.xypos]);
    Z = squareform(D); [~,Same] = min(Z(1,2:end));
    TileNb = t_save_value0(Same);
    
    % replace bad tile with better z-focused tile 
    copyfile(fullfile(o.TileDirectory,[o.FileBase{RR}, '_t', num2str(TileNb), '.tif']),...
        fullfile(o.TileDirectory,[o.FileBase{RR}, '_t', num2str(TileNb), '_wrongZ.tif']));
     copyfile(fullfile(o.TileDirectory,[o.FileBase{r}, '_t', num2str(t_save_value(i)), '.tif']),...
        fullfile(o.TileDirectory,[o.FileBase{RR}, '_t', num2str(TileNb), '.tif']));
    end
    
end




end
