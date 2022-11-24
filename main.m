clc; clear
fileName = 'Berea';

% Pressure difference (Pa)
inletPressure = 1;
outletPressure = 0;

network = Network(fileName);

tic
network.calculateNetworkProperties(inletPressure, outletPressure);
fprintf('Porosity of the model(in percentage) is: %3.5f \n', network.Porosity);
toc 

tic
fprintf('Drainage Start\n');
network.PrimaryDrainage(inletPressure, outletPressure);
toc

tic
fprintf('Imbibition Start\n'); 
network.ScoendaryImbibition(inletPressure, outletPressure); 
toc