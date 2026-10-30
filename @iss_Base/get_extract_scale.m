function o = get_extract_scale(o,nChannels,nZstacks,imfile,scene,SE,DapiSE,r,t_index,t)
%% o = get_extract_scale(o,nChannels,nZstacks,imfile,scene,SE,DapiSE,r,t_index,t)
%
% This finds the scale by which to multiply each filtered image. This is
% based on the max intensity of tile t round r across all color channels.
% It also writes the images for this tile and round to o.TileDirectory

%% Determine channel offset
if  o.IgnoreFirstChannel
    channelOffset = 1; % Skip channel 1 (e.g., DAPI)
else
    channelOffset = 0;
end

fName = fullfile(o.TileDirectory, ...
    [o.FileBase{r}, '_t', num2str(t), '.tif']);

if exist(fName, 'file')
    if r~=o.AnchorRound
        fprintf('Tile %d, from which o.ExtractScale derived already exists.\n',o.ExtractScaleTile);
    else
        fprintf('Tile %d, from which o.ExtractScaleAnchor derived already exists.\n',o.ExtractScaleTile);
    end
    
else
    % a new reader per worker
    bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
    % use the memo file cached before
    bfreader.setId(imfile);
    
    bfreader.setSeries(scene*t_index-1);
    
    if r~=o.AnchorRound; fprintf('Obtaining ExtractScale for imaging rounds from tile %d... '...
            ,o.ExtractScaleTile);
    else; fprintf('Obtaining ExtractScale for anchor round from tile %d... '...
            ,o.ExtractScaleTile);
    end
    
    MaxPixelValue = 0;
    RawImage_all = zeros(o.TileSz,o.TileSz,nChannels);
    IFS_all = zeros(o.TileSz,o.TileSz,nChannels);
    
    for c=1:nChannels
        file_c = c + channelOffset; % Map loop channel to file channel
        
        I = cell(nZstacks,1);
        for z = 1:nZstacks
            iPlane = bfreader.getIndex(z-1, file_c-1, 0)+1;
            I{z} = bfGetPlane(bfreader, iPlane);
        end
        
        % focus stacking
        I_mod = o.fstack_modified(I(o.FirstZPlane:end));
        
        if o.StripHack
            [I_mod,BadColumns] = o.StripHack_raw(I_mod);
        end
        RawImage_all(:,:,c) = I_mod;
        
        % Check against physical file channels
        if (file_c == o.DapiChannel  && r == o.DapiRound) ||...
                ((file_c == o.GcampChannel || file_c == o.GadChannel) &&...
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
            IFS_all(:,:,c) = IFS;
        else
            I_mod = single(padarray(I_mod,(size(SE)-1)/2,'replicate','both'));
            IFS = convn(I_mod,SE,'valid');
            if o.StripHack
                %Set faulty columns all to be 0 so peaks can't be found there.
                IFS(:,BadColumns) = 0;
            end
            if max(IFS(:))>MaxPixelValue
                MaxPixelValue = max(IFS(:));
            end
            IFS_all(:,:,c) = IFS;
        end
        clearvars I_mod I IFS %Free up memory
    end
    
    if r==o.AnchorRound
        o.ExtractScaleAnchor = o.ExtractScaleNorm/MaxPixelValue;
        % Exclude special channels from scaling if they map correctly
        skipChannels = [];
        if isprop(o, 'DapiChannel') && ~isempty(o.DapiChannel), skipChannels = [skipChannels, o.DapiChannel - channelOffset]; end
        if isprop(o, 'GadChannel') && ~isempty(o.GadChannel), skipChannels = [skipChannels, o.GadChannel - channelOffset]; end
        if isprop(o, 'GcampChannel') && ~isempty(o.GcampChannel), skipChannels = [skipChannels, o.GcampChannel - channelOffset]; end
        
        validScaleIdx = setdiff(1:nChannels, skipChannels);
        IFS_all(:,:,validScaleIdx) = IFS_all(:,:,validScaleIdx) * o.ExtractScaleAnchor;
    else
        o.ExtractScale = o.ExtractScaleNorm/MaxPixelValue;
        IFS_all = IFS_all*o.ExtractScale;
    end
    
    if r~=o.AnchorRound; fprintf('Value is %.2f\n', o.ExtractScale);
    else; fprintf('Value is %.2f\n', o.ExtractScaleAnchor);end
    
    for c=1:nChannels
        file_c = c + channelOffset;
        
        if (file_c == o.DapiChannel  && r == o.DapiRound) ||...
                ((file_c == o.GcampChannel || file_c == o.GadChannel) &&...
                (r == o.DapiRound || r == o.AnchorRound))
            I_mod = uint16(IFS_all(:,:,c));
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
            if file_c ~= o.AnchorChannel && r == o.AnchorRound
                IFS = uint16(IFS_all(:,:,c)+o.TilePixelValueShift);
            else
                %Determine auto thresholds
                IFS = IFS_all(:,:,c);
                o.AutoThresh(t,c,r) = median(abs(IFS(:)))*o.AutoThreshMultiplier;
                if ismember(r,1:o.nRounds)
                    %Get histogram data
                    IFS = int32(IFS);
                    o.HistCounts(:,c,r) = o.HistCounts(:,c,r)+histc(IFS(:),o.HistValues);
                end
                if o.StripHack
                    %Set faulty columns all to be 0 so peaks can't be found there.
                    IFS(:,BadColumns) = 0;
                end
                IFS = uint16(IFS+o.TilePixelValueShift);
            end
        end
        
        imwrite(IFS,...
            fullfile(o.TileDirectory,...
            [o.FileBase{r}, '_t', num2str(t), '.tif']),...
            'tiff', 'writemode', 'append');
        fprintf('Round %d tile %d colour channel %d finished.\n', r, t, c);
    end
    
    if nChannels < length(o.bpLabels) && r == o.DapiRound
        for c = (nChannels + 1) : length(o.bpLabels)
            imwrite(IFS,...
                fullfile(o.TileDirectory,...
                [o.FileBase{r}, '_t', num2str(t), '.tif']),...
                'tiff', 'writemode', 'append');
        end
        fprintf(['\n There was less channels than codes for DAPI round: ',num2str(r),....
            ' , replacing missing channels with Anchor'])
    end
end

%Save o object so keep ExtractScale values
save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');

end