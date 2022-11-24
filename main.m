clc; clear
fileName = 'simplest'; 

% Pressure difference (Pa)
inletPressure = 1;
outletPressure = 0;

network = Network(fileName);
tic
network.calculateNetworkProperties(inletPressure, outletPressure);
% networkInfo(network)
% visualized(network) 
toc

tic
soluteConcentration = 1;
poreVolumeInjected = network.poreVolume/50; 
simulationVolume = network.poreVolume/10;
fprintf('=============================== Diffusion Start =====================================\n'); 
network.calculateReactiveTransport_SinglePhaseDiffusion(inletPressure, outletPressure, soluteConcentration, simulationVolume, poreVolumeInjected);
% Raoof
fprintf('=============================== Diffusion Start =====================================\n'); 
network.calculateReactiveTransport_SinglePhaseDiffusion_R(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
% fprintf('=============================== Desorption Start =====================================\n'); 
% network.calculateReactiveTransport_SinglePhaseDesorption(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
toc

% Visualization  
fileName = network.vtkWriter_glyph('init',0);
