clc; clear
% fileName = 'CARB';
fileName = 'cylindrical_101010';

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
% network.pressureDistribution_Cylindrical (inletPressure, outletPressure);
% toc

% tic
% network.calculateAbsolutePermeability(inletPressure, outletPressure);
% fprintf('Permeability of the model(in mD) is: %3.5f \n', network.absolutePermeability);
% toc

tic
soluteConcentration = 1;
poreVolumeInjected = 1* network.poreVolume;
% network.calculateReactiveTransport_SinglePhase_circle(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
network.calculateReactiveTransport_SinglePhase(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
% network.calculateReactiveTransport_Edgar(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
% network.calculateReactiveTransport_SinglePhaseDiffusion(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
toc
