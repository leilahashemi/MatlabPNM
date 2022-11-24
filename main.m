clc; clear
fileName = 'Berea';

% Pressure difference (Pa)
inletPressure = 1;
outletPressure = 0;

network = Network(fileName);

tic
network.calculateNetworkProperties(inletPressure, outletPressure);
networkInfo(network)
toc 
