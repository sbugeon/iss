clear
SliceName = 'Round-';
MainF =['D:\ISS\446003\Slide5\',SliceName];

load(fullfile(MainF,'output\oExtract.mat'))
oOut = o.OutputDirectory;

% process the different regions (if there are some)
for j = 1:length(o.TileConnectedID)
   
    o.OutputDirectory = fullfile(oOut,num2str(j));
    load(fullfile(o.OutputDirectory,'oExtract'));
        
    %% add dapi channel info
    o.DapiChannel = 1;
    % o.DapiRound = 8;
    o.DapiRound = 9;
    o.AnchorChannel = 3;
    o.ReferenceChannel = 3;
    o.FileBase{9} = strcat(SliceName,'08');
    
    o.nExtraRounds = 2;         %Treat DAPI and Anchor channel as extra round
    
    o.RawDAPIGad = 1; % if set to 1, the DAPI images will not be filtered
    %% set ExtractScale with some arbitrary values
    o.ExtractScale = 9.3;
    o.ExtractScaleAnchor = 8.6;
    
    if size(o.AutoThresh,4) ~= o.nRounds + o.nExtraRounds
        o.AutoThresh = cat(3,o.AutoThresh,zeros(size(o.AutoThresh,1),size(o.AutoThresh,2),1));
    end
    %% extract and filter dapi images
    o = o.extract_and_filter;
    save(fullfile(o.OutputDirectory,'oExtractDAPI.mat'),'o','-v7.3')
    %%
    MakeDapi_images
end



