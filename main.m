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

% tic
% soluteConcentration = 1;
% poreVolumeInjected = network.poreVolume/15; 
% fprintf('=============================== Diffusion Start =====================================\n'); 
% network.calculateReactiveTransport_SinglePhaseDiffusion(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
% fprintf('=============================== Desorption Start =====================================\n'); 
% network.calculateReactiveTransport_SinglePhaseDesorption(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
% toc

% Visualization 
% fileName = network.vtkWriter('init',0);
fileName = network.vtkWriter_glyph('init',0);
 
% tic
fprintf('=============================== Drainage Start =====================================\n'); 
network.PrimaryDrainage(inletPressure, outletPressure);
% % Plot Pc & Kr 
% PLOTDRAIN(network);
% toc
%  
% tic
fprintf('=============================== Imbibition Start =====================================\n');  
% network.contactAngleDistribution();
network.ScoendaryImbibition(inletPressure, outletPressure); 
% % Plot Pc & Kr 
% PLOTDRAIN_IMB(network);
% toc