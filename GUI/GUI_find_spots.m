function app = GUI_find_spots(app,o)

% find spots

%Search paramaters
o.FindSpotsMinScore = 'auto';
o.FindSpotsStep = [5,5];
%FindSpotsSearch can either be a 1x1 struct or a o.nRounds x 1 cell of
%structs - have a different range for each round:
o.FindSpotsSearch = struct();
o.FindSpotsSearch.Y = -150:o.FindSpotsStep(1):150;
o.FindSpotsSearch.X = -150:o.FindSpotsStep(2):150;
%Make WidenSearch larger if you think you have a large shift between rounds
o.FindSpotsWidenSearch = [50,50];

o.PcDist = 3;
o.PointCloudMethod = 1;     %1 or 2, set to 2 if no anchor round.

o.DetectionThresh='auto';
o.IsolationThresh ='auto';
%run code

%     nMatches<o.PcMinSpots | o.AllBaseSpotNo<o.PcMinSpots;
% o.PcMinSpotsScaling
o = o.find_spots2;
save(fullfile(o.OutputDirectory, 'oFind_spots'), 'o', '-v7.3');

app.o = o;