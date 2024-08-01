function o = get_TilePos(o, xypos, nTiles)
% This function takes the tile coordinates and derives the Y,X index of
% each tile. Adds o.TileInitialPosYX and o.TilePosYX to o object.
% o: iss object
% xypos: [X,Y] coordinate of each tile as given from metadata.
% nTiles: number of tiles

o.TilePosYX = zeros(size(xypos));
o.TileInitialRawPosYX = xypos;
%%%%%%%% find regions that are connected, split the ones that aren't
o = o.getConnectedRegions;

if strcmp(o.RawFileExtension,'.czi')
    B = max(xypos);
    
    xypos0 = B - xypos;
end

% xypos0 = xypos;
for i = 1:length(o.TileConnectedID)
    
    SubTiles = o.TileConnectedID{i};
    xypos = xypos0(SubTiles,:);
    nTiles = length(SubTiles);
    
    % find x and y grid spacing as median of distances that are about
    % right
    dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
    if max(abs(dx(:)))==0
        xStep = o.MicroscopeStepSize;
    else
        xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
    end
    dy = xypos(:,2)-xypos(:,2)'; % all pairs of y distances
    
    if max(abs(dy(:)))<o.MicroscopeStepSize/3
        yStep = o.MicroscopeStepSize;
    else
        yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));
    end
    
    
    % find index for each tile
    %     if isempty(o.TileInitialPosYX)
    if nTiles==1
        o.TileInitialPosYX = [1,1];
    else
        o.TileInitialPosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));
    end
    %     end
    
    % only consider tiles with index less than number of tiles and only 1
    % away from another tile index.
    uniqueY = unique(o.TileInitialPosYX(:,1),'sorted');
    uniqueX = unique(o.TileInitialPosYX(:,2),'sorted');
    Y = uniqueY(diff([0;uniqueY])==1 & uniqueY<=nTiles);
    X = uniqueX(diff([0;uniqueX])==1 & uniqueX<=nTiles);
    MaxY = max(Y);
    MaxX = max(X);
    if isempty(MaxX)
        MaxX = NaN;
    end
    if isempty(MaxY)
        MaxY = NaN;
    end
    
    %Sometimes get Nan, if only one Nan, then check if all tiles
    %arranged along only one direction i.e. Nan should be 1.
    if max(isnan([MaxX,MaxY])) && nanmax(MaxX,MaxY)==nTiles
        if isnan(MaxX); MaxX=1; else; MaxY=1; end
    end
    
    if MaxY*MaxX ~= nTiles & strcmp(o.RawFileExtension,'.czi')
        error('Need to write a fix for czi')
    end
    if MaxY*MaxX ~= nTiles
        %Hack if one of the dimensions is correct
        Ydiff = unique(diff(Y));
        Xdiff = unique(diff(X));
        if length(Ydiff)==1 && Ydiff==1 && mod(nTiles/MaxY,1)==0
            MaxX = nTiles/MaxY;
        elseif length(Xdiff)==1 && Xdiff==1 && mod(nTiles/MaxX,1)==0
            MaxY = nTiles/MaxX;
        else
            % save workspace if error for debugging
            save(fullfile(o.OutputDirectory, 'ExtractError_TilePosErrorWorkspace'));
            error('Number of tiles (%d) is not equal to maximum Y position (%d) multiplied by maximum X position (%d)'...
                , nTiles, MaxY, MaxX);
        end
    end
    
    if strcmp(o.RawFileExtension,'.nd2')
        %Once found correct MaxY and MaxX, we know what the answer should be.
        o.TilePosYX = zeros(size(o.TileInitialPosYX));
        TilePosY = flip(repelem(1:MaxY,MaxX));
        o.TilePosYX(:,1) = TilePosY;
        TilePosX = repmat([flip(1:MaxX),1:MaxX],1,ceil(MaxY/2));
        o.TilePosYX(1:nTiles,2) = TilePosX(1:nTiles);
    elseif strcmp(o.RawFileExtension,'.czi') %%%%%%%%%%%%
        M = max(o.TilePosYX);
        o.TilePosYX(SubTiles,:) = o.TileInitialPosYX + [M(1) 0];
        dd = o.TileInitialPosYX + [M(1) 0];
        if any(isnan(dd(:)))
            error
        end
    end
    
end

end

