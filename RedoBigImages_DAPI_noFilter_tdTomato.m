MainF = '\\zaru\Subjects\ISS\Sara\1130494\Slide3\';
SliceL = {'Slice_001','Slice_002','Slice_003','Slice_004','Slice_005'};


for iSlice=1:length(SliceL)
    % load(fullfile(o.OutputDirectory,'oExtract.mat'))
    SliceNb = SliceL{iSlice};
    load(fullfile(MainF,SliceNb,'output\oExtract.mat'))
    
    %% Parameters that should be checked before each run
    
    o.InitialShiftChannel = 7;      %Channel to use to find initial shifts between rounds
    o.ReferenceRound = 8;           %Global coordinate system is built upon o.ReferenceRound and
    o.ReferenceChannel = 7;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
    o.AnchorChannel=7;
    o.RawFileExtension = '.nd2';    %Format of raw data
    o.LogToFile = 0;                %Set to 1 if you want to save command window to txt file, else set to 0.
    
    o.GadChannel = 4;
    o.GadRound = 8;
    o.GcampChannel = 6;
    o.GcampRound = 8;
    o.DapiChannel = 1;
    % o.DapiRound = 8;
    o.DapiRound = 9;
    
    o.FileBase{9} = strcat(SliceNb,'_DAPI');
    %% File Names
    %CHECK BEFORE EACH RUN
    %Folder path of raw data
%     o.FileBase{9} = strcat(SliceNb,'_DAPI');
    o.RawFileExtension = '.nd2';
    %Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
    %the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
    o.CodeFile = strcat('F:\Data-Analysis\SB010\ISS\codebook\codebook_7rounds.txt');
    %% extract and filter
    
    %parameters
    o.nRounds = 7;
    o.nExtraRounds = 2;         %Treat Anchor channel as extra round
    o.FirstBaseChannel = 1;
    o.bpLabels = {'0', '1', '2', '3','4','5','6'}; %order of bases
    
    %These specify the dimensions of the filter. R1 should be approximately the
    %radius of the spot and R2 should be double this.
    o.ExtractR1 = 'auto';
    o.ExtractR2 = 'auto';
    
    %     o.ExtractScale = 2.8258;
    %     o.ExtractScale = 'auto';
    o.TilePixelValueShift = 15000;
    
    %Max time (seconds) to wait for raw .nd2 files to be obtained
    o.MaxWaitTime1 = 60;      %Less time for round 1 incase name is wrong
    o.MaxWaitTime = 21600;
    
    %%%%%%
    o.ExtractScale = 9.3;
    o.ExtractScaleAnchor = 8.6;
    %%%%%
    
    %run code
    if size(o.AutoThresh,4) ~= o.nRounds + o.nExtraRounds
        o.AutoThresh = cat(3,o.AutoThresh,zeros(size(o.AutoThresh,1),size(o.AutoThresh,2),1));
    end
    %%
%     dd = dir(o.TileDirectory); dd = {dd(:).name};
%     dd = dd(contains(dd,'Anchor'));
%     
%     mkdir(fullfile(o.TileDirectory,'old_anchor'))
%     for ff = 1:length(dd)
%         copyfile(fullfile(o.TileDirectory,dd{ff}),fullfile(o.TileDirectory,'old_anchor'))
%         delete(fullfile(o.TileDirectory,dd{ff}))
%     end
    
    
    o = o.extract_and_filter;
    % save(fullfile(o.OutputDirectory,'oExtractDAPI.mat'),'o','-v7.3')
    %%
    MakeDapi_images_fixedV3
    
    % DapiImg = imread(fullfile( o.OutputDirectory,'background_image_fixed.tif'));
    % GadImg = imread(fullfile( o.OutputDirectory,'Gad_image_fixed.tif'));
    % GcampImg = imread(fullfile( o.OutputDirectory,'Gcamp_image_fixed.tif'));
    % AnchorImg = imread(fullfile( o.OutputDirectory,'Anchor_image_fixed.tif'));
    %
    % MultiImg = cat(3,imadjust(GadImg), imadjust(AnchorImg),...
    %     imadjust(DapiImg));
end



