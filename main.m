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

% tic
% network.pressureDistribution_singlePhaseFlow (inletPressure, outletPressure);
% toc

% tic
% network.calculateAbsolutePermeability(inletPressure, outletPressure);
% fprintf('Permeability of the model(in mD) is: %3.5f \n', network.absolutePermeability);
% toc

tic
fprintf('Drainage Start\n');
network.PrimaryDrainage(inletPressure, outletPressure);
toc

% tic
% fprintf('Imbibition Start\n');
% network.ScoendaryImbibition(inletPressure, outletPressure);
% network.ScoendaryImbibition_Spon(inletPressure, outletPressure);
% network.ScoendaryImbibition_Forced(inletPressure, outletPressure);
% toc