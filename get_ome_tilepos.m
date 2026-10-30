function [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize,xyposRaw, imageSize] =...
    get_ome_tilepos(reader)
% [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] =...
%     get_ome_tilepos(reader)
% input: OMERO reader or microscopy image file that contains OMERO metadata
% Xiaoyan, 2017


% get reader if input is file name (character)
if ischar(reader)
    reader = bfGetReader(reader);
end

nSeries = reader.getSeriesCount();
nChannels = reader.getSizeC;

omeMeta = reader.getMetadataStore();
nZstacks = omeMeta.getPixelsSizeZ(0).getValue();
pixelsize = double(omeMeta.getPixelsPhysicalSizeX(0).value());

xypos  = [];
for i = 1:nSeries
    
    % Extract image width (X) and height (Y) in pixels for the current series
    % Note:OME-XML metadata index for pixels is usually 0 if shared, or matches series depending on format,
    % but reader.setSeries(i-1) ensures reader methods target the correct series.
    reader.setSeries(i-1);
    imgWidth = reader.getSizeX();
    imgHeight = reader.getSizeY();
    imageSize = [imgWidth, imgHeight];

    if isempty(omeMeta.getPlanePositionX(i-1, 0)) || isempty(omeMeta.getPlanePositionY(i-1, 0))
        fprintf('Could not load tile positions for series %d\n', i);
    else
        posX = double(omeMeta.getPlanePositionX(i-1, 0).value());
        posY = double(omeMeta.getPlanePositionY(i-1, 0).value());
        xypos = [xypos; posX, posY];
        fprintf('Loaded tile positions for series %d\n', i);
    end
end

if ~isempty(xypos)
    xyposRaw = xypos/pixelsize;
    xypos = bsxfun(@minus, xypos, min(xypos,[],1))/pixelsize;
end

nSerieswPos = size(xypos,1);

% close
reader.close();

end
