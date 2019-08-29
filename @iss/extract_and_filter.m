function o = extract_and_filter(o)
% create tiff files for each tile that are top-hat filtered versions of
% original czi files

    o.TileFiles = cell(o.nRounds+o.nExtraRounds,1,1,1); % 1,1,1 because we don't yet know how many tiles

    for r = 1:o.nRounds+o.nExtraRounds       
        imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);

        % construct a Bio-Formats reader with the Memoizer wrapper
        bfreader = loci.formats.Memoizer(bfGetReader(), 0);
        % initiate reader
        bfreader.setId(imfile);

        % get some basic image metadata
        [nSeries, nSerieswPos, nChannels, o.nZ, xypos, o.XYpixelsize,o.Zpixelsize] = ...
            get_ome_tilepos(bfreader);
        if isempty(xypos) || size(xypos, 1)==1
            if r == 1
                warning('first round xypos empty - using values from initial manual input')
                assert(~isempty(o.TileInitialPosXY), 'xypos unavailable')
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
        
        % find x and y grid spacing as median of distances that are about
        % right
        dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
        xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
        dy = xypos(:,1)-xypos(:,1)'; % all pairs of y distances
        yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));
        
        % find coordinates for each tile
        TilePosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));
        if size(xypos, 1)==1
            TilePosYX = [1,1];
        end
        o.TilePosYXC = zeros(nSerieswPos*nChannels,3);

        % set up filename grid for this round
        fName = cell(nSerieswPos*nChannels,1);
        
        Index = 1;
        %parfor t = 1:nSerieswPos  
        for t = 1:nSerieswPos  
                       
            % a new reader per worker
            bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % use the memo file cached before
            bfreader.setId(imfile);

            bfreader.setSeries(scene*t-1);
            for c = 1:nChannels
                
                fName{Index} = fullfile(o.TileDirectory, ...
                    [o.FileBase{r}, '_t', num2str(t),'c', num2str(c), '.tif']);  
                
                if exist(fName{Index}, 'file')
                    fprintf('Round %d tile %d already done.\n', r, t);
                    o.TilePosYXC(Index,:) = [TilePosYX(t,:),c];          %Think first Z plane is the highest
                    o.TileFiles{r,o.TilePosYXC(Index,1), o.TilePosYXC(Index,2),o.TilePosYXC(Index,3)} = fName{Index};
                    Index = Index+1;
                    continue;
                end
                                                                        
                %TopHat SE
                if c == o.DapiChannel && r == o.ReferenceRound    
                        %SE = strel3D_2(20,10);       % I.e. set to 8 microns for DAPI
                        SE = get_3DSE(o.DapiR1YX,o.DapiR1Z,o.DapiR2YX,o.DapiR2Z);       
                else
                        %SE = strel3D_2(3,3);    %I.e. Set to 1 micron
                        SE = get_3DSE(o.ExtractR1YX,o.ExtractR1Z,o.ExtractR2YX,o.ExtractR2Z);
                end

                I = zeros(o.TileSz,o.TileSz,o.nZ); 
                for z = 1:o.nZ
                    iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
                    I(:,:,z) = bfGetPlane(bfreader, iPlane);
                end       
                
                %Make noise white first by divding amplitude of FT
                tic
                FT = fftn(I);
                Norm_FT = FT ./ abs(FT);
                %filter = fspecial3('gaussian',size(I),2);%DESCRIBE BETTER!!!!
                %Shiftfilter = fftshift(filter);     %Shift for FT so centered on 0
                %FT_filter = fftn(Shiftfilter);          
                %NormFT_filter = FT_filter ./ abs(FT);
                %Final_FT = Norm_FT .* NormFT_filter;
                %IFS = ifftn(Final_FT);
                
                I = ifftn(Norm_FT);
                I = padarray(I,(size(SE)-1)/2,'replicate','both');
                IFS = convn(I,SE,'valid');
                IFS = IFS*o.ExtractScale;
                toc

                % tophat the 3D image
                %IFS = imtophat(I, SE);
                
                %Pad by replicating edges
                %tic
                %I = padarray(I,(size(SE)-1)/2,'replicate','both');
                %IFS = convn(I,SE,'valid');
                %Offset = 0;                %To remove some of the negative numbers
                                            %Not set to min(IFS) as some images have 
                                            %much lower so can't have universal threshold
                                          
                %IFS = (IFS+Offset)*100;      %scale so get more info when rounded to int
                %toc
                
                %tic
                %IFS = imfilter(I,SE,'replicate','same','conv');
                %toc
                %IFS = I;
                
                %Append each z plane to same tiff image
                for z = 1:o.nZ
                    imwrite(uint16(IFS(:,:,z)),...  %Not sure if uint16 is correct, wasnt working without
                            fullfile(o.TileDirectory,...
                            [o.FileBase{r}, '_t', num2str(t),'c', num2str(c), '.tif']),...
                            'tiff', 'writemode', 'append');
                end

                o.TilePosYXC(Index,:) = [TilePosYX(t,:),c];          %Think first Z plane is the highest
                o.TileFiles{r,o.TilePosYXC(Index,1), o.TilePosYXC(Index,2),o.TilePosYXC(Index,3)} = fName{Index};
                fprintf('Round %d tile %d colour channel %d finished.\n', r, t, c);                                               
                Index = Index+1; 
               
            end
            bfreader.close();
            
        end
        
    
    o.EmptyTiles = cellfun(@isempty, squeeze(o.TileFiles(o.ReferenceRound,:,:,1)))*0;

    end
end

function SE = get_3DSE(r1YX,r1Z,r2YX,r2Z)
    % structuring element for convlolution filtering
    % Positive inner circle radius r1 and negative outer annulus radius r2. Overall sums to
    % zero. 
    SE = zeros(r2YX*2+1,r2YX*2+1,r2Z*2+1);
    SE(r2YX+1-r1YX:r2YX+1+r1YX,r2YX+1-r1YX:r2YX+1+r1YX,r2Z+1-r1Z:r2Z+1+r1Z) = fspecial3('ellipsoid',[r1YX,r1YX,r1Z]);
    SE = SE - fspecial3('ellipsoid',[r2YX,r2YX,r2Z]);
end