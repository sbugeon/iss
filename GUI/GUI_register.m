function app = GUI_register(app,o)

TileOverlap = 0.1;
D = o.TileSz * (1- TileOverlap);

% register
%parameters
%Anchor spots are detected in register2
o.DetectionRadius = 2;
o.SmoothSize = 0;
o.IsolationRadius1 = 4;
o.IsolationRadius2 = 14;

o.DetectionThresh = 'auto';
o.ThreshParam = 5;
o.MinThresh = 10;
o.minPeaks = 1;
o.InitalShiftAutoMinScoreParam = 2;   %a lower value will make it quicker but more likely to fail

%paramaters to find shifts between overlapping tiles
o.RegMinScore = 'auto';
o.RegStep = [5,5];
o.RegSearch.South.Y = -1.05*D:o.RegStep(1):-0.95*D;
o.RegSearch.South.X = -0.05*o.TileSz:o.RegStep(2):0.05*o.TileSz;
o.RegSearch.East.Y = -0.05*o.TileSz:o.RegStep(1):0.05*o.TileSz;
o.RegSearch.East.X = -1.05*D:o.RegStep(2):-0.95*D;

o.RegWidenSearch = [50,50];

%If a channel or round is faulty, you can ignore it by selecting only the
%good ones in o.UseChannels and o.UseRounds.
o.UseChannels = 1:o.nBP;
o.UseRounds = 1:o.nRounds;

o.RegMethod='Fft';
o = o.register2;
save(fullfile(o.OutputDirectory, 'oRegister'), 'o', '-v7.3');

app.o = o;