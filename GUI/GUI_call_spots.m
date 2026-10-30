function app = GUI_call_spots(app,o)

% call spots
%run code
o.CallSpotsCodeNorm = 'WholeCode';      %Other alternative is 'Round'
o = o.call_spots;

%OMP
o.ompInitialNeighbThresh = 5;  % Increase to use less memory. Keep below 10.
o = o.call_spots_omp;
save(fullfile(o.OutputDirectory, 'oCall_spots_OMP'), 'o', '-v7.3');

app.o = o;