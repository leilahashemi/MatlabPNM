clc; clear
fileName = 'cylindrical_101010';

% Pressure difference (Pa)
inletPressure = 1;
outletPressure = 0;

network = Network(fileName);

tic
network.calculateNetworkProperties(inletPressure, outletPressure);
fprintf('Porosity of the model(in percentage) is: %3.5f \n', network.Porosity);
toc 

tic
soluteConcentration = 1;
poreVolumeInjected = network.poreVolume/10; 
network.calculateReactiveTransport_SinglePhaseDiffusion(inletPressure, outletPressure, soluteConcentration, poreVolumeInjected);
toc

tic
network.vtkOutput();
network.vtkOutput_old();
network.vtpOutput();
toc