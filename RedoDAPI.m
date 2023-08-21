clear
SliceName = 'B3S4';
SliceList = {[SliceName,'_Slice001'],[SliceName,...
    '_Slice002'],[SliceName,'_Slice003'],[SliceName,'_Slice004']};
MainF =repmat({['\\zaru\Subjects\SB038\ISS\',SliceName]},1,4);

for iSlice = 1:4
    SliceNb = SliceList{iSlice};

    load(fullfile(MainF,SliceNb,'output\oExtract.mat'))
    %% add dapi channel info
    o.DapiChannel = 1;
    % o.DapiRound = 8;
    o.DapiRound = 9;
    
    o.FileBase{9} = strcat(SliceNb,'_DAPI');
   
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



