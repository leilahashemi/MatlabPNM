function visualized(network)
%% Fluid info
fprintf('===============================Ntework Information=====================================\n');
fprintf('Water Viscosity:                        %f \n', network.waterViscosity); 
fprintf('Water Viscosity:                        %f \n', network.oilViscosity); 
fprintf('Interfacial Tension:                    %f \n', network.sig_ow); 
fprintf('======================================================================================\n\n'); 
%% Plot 3D
figure('name','PorePressure')
N = nan(network.numberOfNodes ,4); % coordinate & radius
M = nan(network.numberOfNodes ,3); % color
for i = 1:network.numberOfNodes
    
    N (i,1:3)=[network.Nodes{i}.x_coordinate,network.Nodes{i}.y_coordinate,network.Nodes{i}.z_coordinate];
    N (i,4)=network.Nodes{i}.radius*10000000;
    if network.Nodes{i}.waterPressure > 0.9
        M(i,:) = [0,0,0.1];
    elseif network.Nodes{i}.waterPressure <= 0.9 && network.Nodes{i}.waterPressure > 0.8
        M(i,:) = [0,0,0.2];
    elseif network.Nodes{i}.waterPressure <= 0.8 && network.Nodes{i}.waterPressure > 0.7
        M(i,:) = [0,0,0.3];
    elseif network.Nodes{i}.waterPressure <= 0.7 && network.Nodes{i}.waterPressure > 0.6
        M(i,:) = [0,0,0.4];
    elseif network.Nodes{i}.waterPressure <= 0.6 && network.Nodes{i}.waterPressure > 0.5
        M(i,:) = [0,0,0.5];
    elseif network.Nodes{i}.waterPressure <= 0.5 && network.Nodes{i}.waterPressure > 0.4
        M(i,:) = [0,0,0.6];
    elseif network.Nodes{i}.waterPressure <= 0.4 && network.Nodes{i}.waterPressure > 0.3
        M(i,:) = [0,0,0.7];
    elseif network.Nodes{i}.waterPressure <= 0.3 && network.Nodes{i}.waterPressure > 0.2
        M(i,:) = [0,0,0.8];
    elseif network.Nodes{i}.waterPressure <= 0.2 && network.Nodes{i}.waterPressure > 0.1
        M(i,:) = [0,0,0.9];
    else
        M(i,:) = [0,0,1]; 
    end
end
scatter3(N(:,1),N(:,2),N(:,3),N(:,4),M,'filled');
hold on
% figure('name',' Links')
B = nan(network.numberOfLinks ,4);
D = nan(network.numberOfLinks ,3);
for i = 1:network.numberOfLinks
    pore1Index = network.Links{i}.pore1Index;
    pore2Index = network.Links{i}.pore2Index;
    B (i,4)=network.Links{i}.radius*5000000;
    if network.Links{i}.isInlet
        B (i,1:3)=[network.Nodes{pore2Index}.x_coordinate,network.Nodes{pore2Index}.y_coordinate,network.Nodes{pore2Index}.z_coordinate];
    elseif network.Links{i}.isOutlet   
        B (i,1:3)=[network.Nodes{pore1Index}.x_coordinate+network.Links{i}.linkLength/2,network.Nodes{pore1Index}.y_coordinate,network.Nodes{pore1Index}.z_coordinate];
    else
        B (i,1:3)=[(network.Nodes{pore1Index}.x_coordinate+network.Nodes{pore2Index}.x_coordinate)/2,...
            (network.Nodes{pore1Index}.y_coordinate+network.Nodes{pore2Index}.y_coordinate)/2,...
            (network.Nodes{pore1Index}.z_coordinate+network.Nodes{pore2Index}.z_coordinate)/2]; 
    end
    if network.Links{i}.waterPressure > 0.9
        D(i,:) = [0,0,0.1];
    elseif network.Links{i}.waterPressure <= 0.9 && network.Links{i}.waterPressure > 0.8
        D(i,:) = [0,0,0.2];
    elseif network.Links{i}.waterPressure <= 0.8 && network.Links{i}.waterPressure > 0.7
        D(i,:) = [0,0,0.3];
    elseif network.Links{i}.waterPressure <= 0.7 && network.Links{i}.waterPressure > 0.6
        D(i,:) = [0,0,0.4];
    elseif network.Links{i}.waterPressure <= 0.6 && network.Links{i}.waterPressure > 0.5
        D(i,:) = [0,0,0.5];
    elseif network.Links{i}.waterPressure <= 0.5 && network.Links{i}.waterPressure > 0.4
        D(i,:) = [0,0,0.6];
    elseif network.Links{i}.waterPressure <= 0.4 && network.Links{i}.waterPressure > 0.3
        D(i,:) = [0,0,0.7];
    elseif network.Links{i}.waterPressure <= 0.3 && network.Links{i}.waterPressure > 0.2
        D(i,:) = [0,0,0.8];
    elseif network.Links{i}.waterPressure <= 0.2 && network.Links{i}.waterPressure > 0.1
        D(i,:) = [0,0,0.9];
    else
        D(i,:) = [0,0,1]; 
    end
end
scatter3(B(:,1),B(:,2),B(:,3),B(:,4),D,'filled');
end