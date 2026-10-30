function app = GUI_extractANDfilter(app)
o = iss_OMP;

Vars = {'AnchorChannel','AnchorRound','GadChannel','GadRound','GcampChannel','GcampRound',...
    'DapiChannel','DapiRound','InitialShiftChannel','InputDirectory','nBP','nRounds','nExtraRounds','IgnoreFirstChannel'};
for i=1:length(Vars)
    o.(Vars{i}) = app.([Vars{i},'EditField']).Value;
end

SliceNb = app.SliceprefixEditField.Value;
RoundSuffix = app.RoundsuffixEditField.Value;
RoundSuffix = split( RoundSuffix, ',');
if length(RoundSuffix) ~= o.nRounds + o.nExtraRounds
    error('\n Number of total rounds is not equal to the number of provided suffixes')
end

%FileBase{r} is the file name of the raw data of round r in o.InputDirectory
o.FileBase = cell(1,1);
for i=1:length(RoundSuffix)
    o.FileBase{i} = strcat(SliceNb,RoundSuffix{i});
end

OutDir = app.OutputDirectoryEditField.Value;

o.TileDirectory = fullfile(OutDir,SliceNb,'\tiles');
mkdir(o.TileDirectory);
o.OutputDirectory = fullfile(OutDir,SliceNb,'\output');
mkdir(o.OutputDirectory);

o.ReferenceRound = o.AnchorRound;           %Global coordinate system is built upon o.ReferenceRound and
o.ReferenceChannel = o.AnchorChannel;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
o.LogToFile = 1;                %Set to 1 if you want to save command window to txt file, else set to 0.
o.StripHack = true;

%Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
%the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
o.CodeFile = 'C:\Users\bugeon\Documents\codebook_73g_OtrAdra_Cholinergic_IEG_final.txt';
%% Logging
if o.LogToFile
    if isempty(o.LogFile)
        o.LogFile = fullfile(o.OutputDirectory,'Log.txt');
    end
end
%% extract and filter

%parameters
o.FirstBaseChannel = 1;
o.bpLabels = {'0', '1', '2','3', '4','5','6'}; %order of bases %%%%%%%%%%%%%%%

o.ExtractR1 = 'auto';
o.ExtractR2 = 'auto';
o.ExtractScale = 'auto';
o.TilePixelValueShift = 15000;

o.MaxWaitTime1 = 0;      %Less time for round 1 incase name is wrong
o.MaxWaitTime = 0;
o = o.extract_and_filter;

save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');

app.o = o;