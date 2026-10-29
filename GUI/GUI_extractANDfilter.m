function GUI_extractANDfilter(app)
    o = iss_OMP;
   
    
    
    o.AnchorChannel = app;
    o.AnchorRound = 8;            %Channel that has most spots in o.AnchorRound
    o.GadChannel = 4;
    o.GadRound = 8;
    o.GcampChannel = 5;
    o.GcampRound = 8;
    o.DapiChannel = 1;
    o.DapiRound = 8;             %Channel in o.AnchorRound that contains Dapi images
    
    o.InitialShiftChannel = 5;      %Channel to use to find initial shifts between rounds
    o.ReferenceRound = o.AnchorRound;           %Global coordinate system is built upon o.ReferenceRound and
    o.ReferenceChannel = o.AnchorChannel;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
    o.RawFileExtension = '.nd2';    %Format of raw data
    o.LogToFile = 0;                %Set to 1 if you want to save command window to txt file, else set to 0.
    o.StripHack = true;
    
    o.TileInitialPosYX = [3,3;3,2;3,1;2,1;2,2;2,3;1,3;1,2;1,1]; % hack for this experiment tile positions
    %% File Names
    %CHECK BEFORE EACH RUN
    o.InputDirectory = MainFolder{iSlice};     %Folder path of raw data
    
    o.TileSz = 3200;        %Dimension of tile in pixels
    o.nBP = 7;              %Number of Channels
    o.nRounds = 7;          %Number of Imaging Rounds
    o.nExtraRounds = 1;     %Treat Anchor channel as extra round
    
    %FileBase{r} is the file name of the raw data of round r in o.InputDirectory
    o.FileBase = cell(1,1);
    o.FileBase{1} = strcat(SliceNb,'R0');
    o.FileBase{2} = strcat(SliceNb,'R1');
    o.FileBase{3} = strcat(SliceNb,'R2');
    o.FileBase{4} = strcat(SliceNb,'R3');
    o.FileBase{5} = strcat(SliceNb,'R4');
    o.FileBase{6} = strcat(SliceNb,'R5');
    o.FileBase{7} = strcat(SliceNb,'R6');
    o.FileBase{8} = strcat(SliceNb,'anchordapi');

    o.TileDirectory = fullfile(o.InputDirectory,SliceNb,'\tiles');
    mkdir(o.TileDirectory);
    o.OutputDirectory = fullfile(o.InputDirectory,SliceNb,'\output');
    mkdir(o.OutputDirectory);
    
    o.RawFileExtension = '.nd2';
    %Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
    %the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
%     o.CodeFile = 'C:\Users\bugeon\Documents\data_coppaFISH\codebook_73g_ctx.txt';
     o.CodeFile = 'C:\Users\bugeon\Documents\codebook_BaptisteHumanChimp.txt';%%%%%%%%%%%%
    %% Logging
    if o.LogToFile
        if isempty(o.LogFile)
            o.LogFile = fullfile(o.OutputDirectory,'Log.txt');
        end
    end
    %% extract and filter
    
    %parameters
    o.FirstBaseChannel = 1;
%     o.bpLabels = {'0', '2', '1','3', '4','5','6'}; %order of bases %%%%%%%%%%%%%%%
    o.bpLabels = {'0', '1', '2','3', '4','5','6'}; %order of bases %%%%%%%%%%%%%%%
    %These specify the dimensions of the filter. R1 should be approximately the
    %radius of the spot and R2 should be double this.
    o.ExtractR1 = 'auto';
    o.ExtractR2 = 'auto';
    
    %     o.ExtractScale = 2.8258;
    o.ExtractScale = 'auto';
    o.TilePixelValueShift = 15000;
    
    %Max time (seconds) to wait for raw .nd2 files to be obtained
    o.MaxWaitTime1 = 60;      %Less time for round 1 incase name is wrong
    o.MaxWaitTime = 21600;
    
    %     o.ExtractScale = 12.9230;
    %run code
    o = o.extract_and_filter;
    
    save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');