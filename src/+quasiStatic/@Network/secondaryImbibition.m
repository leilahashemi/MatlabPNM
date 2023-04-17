function secondaryImbibition(network)

%counter for invaded elements
numOfLinks_SnapOff = 0;
numOfLinks_PistoneLike = 0;
numOfLinks_LayerCollapse = 0;
numOfNodes_SnapOff = 0;
numOfNodes_PoreBodyFilling = 0;
numOfNodes_LayerCollapse = 0;
network.thresholdPressure_Pa = zeros(network.numberOfLinks, 14);
Pc_imb = network.Pc_drain_max_Pa;
Pc_min = Pc_imb;

% Calculating throat Snap-Off & Pistone-Like displacement & layer collapse
for i = 1:network.numberOfLinks
    
    if network.Links{i}.occupancy == 'B' % if the throat is nonWetting filled
        network.Links{i}.calculateThresholdPressurePistonLike_Imbibition (network.Pc_drain_max_Pa);
        network.Links{i}.calculateThresholdPressureSnapOff (network.Pc_drain_max_Pa);
        %         network.Links{i}.imbThresholdPressure_SnapOff
        %         network.Links{i}.calculateThresholdPressureSnapOff_Patzek (network.Pc_drain_max_Pa);
        %         network.Links{i}.imbThresholdPressure_SnapOff
        %         network.Links{i}.calculateThresholdPressureSnapOff_Valvatne (network.Pc_drain_max_Pa);
        %         network.Links{i}.imbThresholdPressure_SnapOff
        if network.Links{i}.isInlet
            network.thresholdPressure_Pa(i,1) = -1;
        elseif network.Links{i}.isOutlet
            network.thresholdPressure_Pa(i,1) = 1;
        end
        network.thresholdPressure_Pa(i,2) = network.Links{i}.imbThresholdPressure_PistonLike;
        network.thresholdPressure_Pa(i,3) = network.Links{i}.imbThresholdPressure_SnapOff;
    end
end

% Calculating Pore Snap-Off & Pore-Body Filling displacement & layer collapse
for i = 1:network.numberOfNodes
    if network.Nodes{i}.occupancy == 'B' % if the throat is nonWetting filled
        network.Nodes{i}.calculateThresholdPressurePoreBodyFilling (network);
        network.Nodes{i}.calculateThresholdPressurePistonLike_Imbibition (network.Pc_drain_max_Pa);
        network.Nodes{i}.calculateThresholdPressureSnapOff (network.Pc_drain_max_Pa);
        if network.Nodes{i}.isInlet
            network.thresholdPressure_Pa(i,8) = -1;
        elseif network.Nodes{i}.isOutlet
            network.thresholdPressure_Pa(i,8) = 1;
        end
        network.thresholdPressure_Pa(i,9) = network.Nodes{i}.imbThresholdPressure_PistonLike;
        network.thresholdPressure_Pa(i,10) = network.Nodes{i}.imbThresholdPressure_SnapOff;
    end
end
t = 0;
network.ImbibitionData = zeros(100,12);

t = t+1;
network.ImbibitionData(t,:) = ...
    [network.DrainageData(end,1),network.DrainageData(end,2), network.DrainageData(end,3),network.DrainageData(end,4),0,0,0,0,0,0,0,Pc_min];
if network.flowVisualization
    network.IO.visualization(network,'SI',t);
end

invaded_Element = zeros(2*(network.numberOfLinks+network.numberOfNodes), 11);
e = 0;
network.sequence = zeros(2*network.numberOfLinks, 11);
percList = -1000000*ones(network.numberOfNodes+network.numberOfLinks,1);
poreVolumeInjected = 0;
[~, NodeL, LinkL, cluster_A_nums, cluster_B_nums] = clustering_nonWettingPhase(network);
inv = false;

% Invasion & Percolation List
for i = 1:network.numberOfLinks
    
    node1Index = network.Links{i}.pore1Index;
    node2Index = network.Links{i}.pore2Index;
    
    if (any(LinkL(i) == cluster_A_nums(:)) || any(LinkL(i) == cluster_B_nums(:)))
        
        if network.Links{i}.isInlet
            
            if any(network.Links{i}.imbThresholdPressure_PistonLike)
                
                percList(i) = network.Links{i}.imbThresholdPressure_PistonLike ;
            end
        elseif network.Links{i}.isOutlet
            
            if  network.Nodes{node1Index}.occupancy == 'A' && any(network.Links{i}.imbThresholdPressure_PistonLike)
                
                percList(i) = network.Links{i}.imbThresholdPressure_PistonLike ;
            elseif  network.Nodes{node1Index}.occupancy == 'B' && ...
                    any(network.Links{i}.imbThresholdPressure_SnapOff) % if the throat is non circular
                
                percList(i) = network.Links{i}.imbThresholdPressure_SnapOff;
            end
        else
            if  (network.Nodes{node1Index}.occupancy == 'A' && network.Nodes{node2Index}.occupancy == 'B') || ...
                    (network.Nodes{node1Index}.occupancy == 'B' && network.Nodes{node2Index}.occupancy == 'A') && ...
                    any(network.Links{i}.imbThresholdPressure_PistonLike)
                
                
                percList(i) = network.Links{i}.imbThresholdPressure_PistonLike ;
                
            elseif network.Nodes{node1Index}.occupancy == 'B' &&...
                    network.Nodes{node2Index}.occupancy == 'B' &&...
                    any(network.Links{i}.imbThresholdPressure_SnapOff)% if the throat is non circular
                
                percList(i) = network.Links{i}.imbThresholdPressure_SnapOff;
            end
        end
    end
end
a = network.numberOfLinks;
for i = 1:network.numberOfNodes
    if (any(NodeL(i) == cluster_A_nums(:)) || any(NodeL(i) == cluster_B_nums(:)))
        
        filledThroats = 0;
        for j = 1:network.Nodes{i}.connectionNumber
            if (network.Links{network.Nodes{i}.connectedLinks(j)}.occupancy == 'A')
                filledThroats = filledThroats + 1;
            end
        end
        
        if filledThroats ~= 0 &&  any(network.Nodes{i}.imbThresholdPressure_PoreBodyFilling)
            
            percList(a+i) = network.Nodes{i}.imbThresholdPressure_PoreBodyFilling;
            
        elseif filledThroats == 0 && any(network.Nodes{i}.imbThresholdPressure_SnapOff)% if the node is non circular
            
            percList(a+i) = network.Nodes{i}.imbThresholdPressure_SnapOff; % snap off threshold pressure
        end
        
    end
end


[~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = clustering_nonWettingPhase(network);
[~, NodeL_W, LinkL_W,cluster_A,cluster_B] = clustering_wettingPhase(network);

% Saturation calculation before invading
Pc_interval_beforeInvading = (Pc_imb - max(percList))/network.NoOfPc_interval;
while (max(percList)) < Pc_imb -1
    Pc_imb = Pc_imb - Pc_interval_beforeInvading;
    t = t+1;
    
    invaded = numOfLinks_SnapOff + numOfLinks_PistoneLike + ...
        numOfNodes_SnapOff + numOfNodes_PoreBodyFilling + numOfNodes_LayerCollapse;
    
    % Updating saturations and conductances
    calculateConductance_and_Saturation_Imbibition(network, Pc_imb,NodeL, NodeL_W, LinkL, LinkL_W, cluster_A_nums, cluster_A, cluster_B_nums, cluster_B);
    
    % Relative Permeability Calculation
    if network.calculateRelativePermeability
        calculatePressureDistribution_TwoPhases(network, network.inletPressure_Pa, network.outletPressure_Pa);
        
        [Kr_w_imb, Kr_nw_imb] =...
            calculateRelativePermeability_Imbibition(network, network.outletPressure_Pa, LinkL, LinkL_W, cluster_A_nums, cluster_A);
    else
        Kr_w_imb =NaN;
        Kr_nw_imb=NaN;
    end
    network.ImbibitionData(t,:) = ...
        [network.wettingPhaseSaturation,Pc_imb,Kr_w_imb, Kr_nw_imb,invaded, ...
        numOfLinks_SnapOff,numOfLinks_PistoneLike, ...
        numOfLinks_LayerCollapse,numOfNodes_SnapOff, ...
        numOfNodes_PoreBodyFilling,numOfNodes_LayerCollapse,Pc_min];
    if network.flowVisualization
        network.IO.visualization(network,'SI',t);
    end
end

NoOfPc_interval = Pc_imb /network.NoOfPc_interval;

while (~isempty(cluster_A_nums) || ~isempty(cluster_B_nums)) && Pc_imb >-99999
    
    press = 1;
    deltaS = 0;
    
    % Percolation Section
    if (max(percList)) >= Pc_imb
        pressure = 1;
    else
        pressure = 0;
    end
    
    while pressure == 1 && (max(percList)) > Pc_imb  && deltaS <= network.deltaS_input
        
        inv = true;
        % Descending sorting of threshold pressures
        [PcTh, ix] = max(percList(1:end));
        if PcTh ~= -1000000
            indexElement = ix(1);
            Pc_min = Pc_imb;
            if Pc_min > percList(indexElement)
                Pc_min = percList(indexElement);
            end
            
            [~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = clustering_nonWettingPhase(network);
            [~, NodeL_W, LinkL_W,cluster_A,cluster_B] = clustering_wettingPhase(network);
            
            % if the first element is a throat
            if indexElement <= network.numberOfLinks
                
                linkIndex = indexElement;
                node1Index = network.Links{linkIndex}.pore1Index;
                node2Index = network.Links{linkIndex}.pore2Index;
                if any(LinkL(linkIndex) == cluster_A_nums(:)) || any(LinkL(linkIndex) == cluster_B_nums(:))
                    
                    if network.Links{linkIndex}.isInlet
                        
                        if network.Links{linkIndex}.imbThresholdPressure_PistonLike >= Pc_imb
                            
                            network.Links{linkIndex}.occupancy = 'A';
                            network.Links{linkIndex}.nonWettingLayerExistance()
                            network.Links{linkIndex}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                            poreVolumeInjected = poreVolumeInjected + network.Links{linkIndex}.volume;
                            numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;
                            network.Links{linkIndex}.isInvaded = true;
                            e = e+1;
                            invaded_Element(e,1:4) = [linkIndex, percList(linkIndex), ...
                                network.Links{linkIndex}.imbThresholdPressure_PistonLike, network.Links{linkIndex}.imbThresholdPressure_SnapOff];
                            
                            if  network.Nodes{node2Index}.occupancy == 'B'
                                % Updating pore body filling of the pore
                                network.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (network);
                                if any(network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling) && ...
                                        Pc_imb <= network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling
                                    network.Nodes{node2Index}.occupancy = 'A';
                                    network.Nodes{node2Index}.nonWettingLayerExistance()
                                    network.Nodes{node2Index}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                                    
                                    poreVolumeInjected = poreVolumeInjected + network.Nodes{node2Index}.volume;
                                    numOfNodes_PoreBodyFilling = numOfNodes_PoreBodyFilling + 1;
                                    network.Nodes{node2Index}.isInvaded = true;
                                    e = e+1;
                                    invaded_Element(e,5:8) = [node2Index, network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling, ...
                                        network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling, network.Nodes{node2Index}.imbThresholdPressure_SnapOff];
                                    percList(network.numberOfLinks+node2Index) = -1000000;
                                    
                                    for j=1:network.Nodes{node2Index}.connectionNumber
                                        if network.Nodes{node2Index}.connectedLinks(j)~=linkIndex
                                            percList(network.Nodes{node2Index}.connectedLinks(j))=...
                                                network.Links{network.Nodes{node2Index}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                                        end
                                    end
                                end
                            end
                        end
                        
                    elseif network.Links{linkIndex}.isOutlet
                        if network.Links{linkIndex}.imbThresholdPressure_PistonLike >= Pc_imb && network.Nodes{node1Index}.occupancy == 'A'
                            
                            network.Links{linkIndex}.occupancy = 'A';
                            network.Links{linkIndex}.nonWettingLayerExistance()
                            network.Links{linkIndex}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                            
                            poreVolumeInjected = poreVolumeInjected + network.Links{linkIndex}.volume;
                            numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;
                            network.Links{linkIndex}.isInvaded = true;
                            e = e+1;
                            invaded_Element(e,1:4) = [linkIndex, percList(linkIndex), ...
                                network.Links{linkIndex}.imbThresholdPressure_PistonLike, network.Links{linkIndex}.imbThresholdPressure_SnapOff];
                            
                        elseif network.Links{linkIndex}.imbThresholdPressure_SnapOff >= Pc_imb && network.Nodes{node1Index}.occupancy == 'B'
                            
                            network.Links{linkIndex}.occupancy = 'A';
                            network.Links{linkIndex}.nonWettingLayerExistance()
                            network.Links{linkIndex}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                            numOfLinks_SnapOff = numOfLinks_SnapOff + 1;
                            poreVolumeInjected = poreVolumeInjected + network.Links{linkIndex}.volume;
                            network.Links{linkIndex}.isInvaded = true;
                            e = e+1;
                            invaded_Element(e,1:4) = [linkIndex, percList(linkIndex), ...
                                network.Links{linkIndex}.imbThresholdPressure_PistonLike, network.Links{linkIndex}.imbThresholdPressure_SnapOff];
                            percList(network.numberOfLinks+node1Index) = -1000000;
                            
                            % Updating pore body filling of the pore
                            network.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (network);
                            if any(network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling) && ...
                                    Pc_imb <= network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling
                                
                                percList(network.numberOfLinks+node1Index) = network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling;
                            end
                        end
                        
                    elseif ~network.Links{linkIndex}.isOutlet && ~network.Links{linkIndex}.isInlet
                        
                        if network.Nodes{node1Index}.occupancy == 'A' || network.Nodes{node2Index}.occupancy == 'A'
                            if network.Links{linkIndex}.imbThresholdPressure_PistonLike >= Pc_imb
                                
                                network.Links{linkIndex}.occupancy = 'A';
                                network.Links{linkIndex}.nonWettingLayerExistance()
                                network.Links{linkIndex}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                                poreVolumeInjected = poreVolumeInjected + network.Links{linkIndex}.volume;
                                numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;
                                network.Links{linkIndex}.isInvaded = true;
                                e = e+1;
                                invaded_Element(e,1:4) = [linkIndex, percList(linkIndex), ...
                                    network.Links{linkIndex}.imbThresholdPressure_PistonLike, network.Links{linkIndex}.imbThresholdPressure_SnapOff];
                                
                                if network.Nodes{node1Index}.occupancy == 'B'
                                    
                                    % Updating pore body filling of the pore
                                    network.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (network);
                                    if any(network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling) && ...
                                            Pc_imb <= network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling
                                        
                                        network.Nodes{node1Index}.occupancy = 'A';
                                        network.Nodes{node1Index}.nonWettingLayerExistance()
                                        network.Nodes{node1Index}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                                        poreVolumeInjected = poreVolumeInjected + network.Nodes{node1Index}.volume;
                                        numOfNodes_PoreBodyFilling = numOfNodes_PoreBodyFilling + 1;
                                        network.Nodes{node1Index}.isInvaded = true;
                                        e = e+1;
                                        invaded_Element(e,5:8) = [node1Index, network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling, ...
                                            network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling, network.Nodes{node1Index}.imbThresholdPressure_SnapOff];
                                        
                                        for j=1:network.Nodes{node1Index}.connectionNumber
                                            if network.Nodes{node1Index}.connectedLinks(j)~=linkIndex
                                                percList(network.Nodes{node1Index}.connectedLinks(j))=...
                                                    network.Links{network.Nodes{node1Index}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                                            end
                                        end
                                    end
                                    percList(network.numberOfLinks+node1Index) = -1000000;
                                end
                                if  network.Nodes{node2Index}.occupancy == 'B'
                                    
                                    % Updating pore body filling of the pore
                                    network.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (network);
                                    if any(network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling) && ...
                                            Pc_imb <= network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling
                                        
                                        network.Nodes{node2Index}.occupancy = 'A';
                                        network.Nodes{node2Index}.nonWettingLayerExistance()
                                        network.Nodes{node2Index}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                                        poreVolumeInjected = poreVolumeInjected + network.Nodes{node2Index}.volume;
                                        numOfNodes_PoreBodyFilling = numOfNodes_PoreBodyFilling + 1;
                                        network.Nodes{node2Index}.isInvaded = true;
                                        e = e+1;
                                        invaded_Element(e,5:8) = [node2Index, network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling, ...
                                            network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling, network.Nodes{node2Index}.imbThresholdPressure_SnapOff];
                                        
                                        for j=1:network.Nodes{node2Index}.connectionNumber
                                            if network.Nodes{node2Index}.connectedLinks(j)~=linkIndex
                                                percList(network.Nodes{node2Index}.connectedLinks(j))=...
                                                    network.Links{network.Nodes{node2Index}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                                            end
                                        end
                                    end
                                    percList(network.numberOfLinks+node2Index) = -1000000;
                                end
                            end
                        elseif network.Links{linkIndex}.imbThresholdPressure_SnapOff >= Pc_imb
                            
                            network.Links{linkIndex}.occupancy = 'A';
                            network.Links{linkIndex}.nonWettingLayerExistance()
                            network.Links{linkIndex}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                            numOfLinks_SnapOff = numOfLinks_SnapOff + 1;
                            
                            % Updating pore body filling of the pore
                            network.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (network);
                            if any(network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling) && ...
                                    Pc_imb <= network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling
                                
                                network.Nodes{node1Index}.occupancy = 'A';
                                network.Nodes{node1Index}.nonWettingLayerExistance()
                                network.Nodes{node1Index}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                                poreVolumeInjected = poreVolumeInjected + network.Nodes{node1Index}.volume;
                                numOfNodes_PoreBodyFilling = numOfNodes_PoreBodyFilling + 1;
                                network.Nodes{node1Index}.isInvaded = true;
                                e = e+1;
                                invaded_Element(e,5:8) = [node1Index,network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling, ...
                                    network.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling, network.Nodes{node1Index}.imbThresholdPressure_SnapOff];
                                
                                for j=1:network.Nodes{node1Index}.connectionNumber
                                    if network.Nodes{node1Index}.connectedLinks(j)~=linkIndex
                                        percList(network.Nodes{node1Index}.connectedLinks(j))=...
                                            network.Links{network.Nodes{node1Index}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                                    end
                                end
                            end
                            
                            % Updating pore body filling of the pore
                            network.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (network);
                            if any(network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling) && ...
                                    Pc_imb <= network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling
                                
                                network.Nodes{node2Index}.occupancy = 'A';
                                network.Nodes{node2Index}.nonWettingLayerExistance()
                                network.Nodes{node2Index}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                                poreVolumeInjected = poreVolumeInjected + network.Nodes{node2Index}.volume;
                                numOfNodes_PoreBodyFilling = numOfNodes_PoreBodyFilling + 1;
                                network.Nodes{node2Index}.isInvaded = true;
                                e = e+1;
                                invaded_Element(e,5:8) = [node2Index, network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling, ...
                                    network.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling, network.Nodes{node2Index}.imbThresholdPressure_SnapOff];
                                
                                for j=1:network.Nodes{node2Index}.connectionNumber
                                    if network.Nodes{node2Index}.connectedLinks(j)~=linkIndex
                                        percList(network.Nodes{node2Index}.connectedLinks(j))=...
                                            network.Links{network.Nodes{node2Index}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                                    end
                                end
                            end
                            percList(network.numberOfLinks+node1Index) = -1000000;
                            percList(network.numberOfLinks+node2Index) = -1000000;
                        end
                    end
                end
                
                % if the first element is a pore
            else
                nodeIndex = indexElement-network.numberOfLinks;
                if any(NodeL(nodeIndex) == cluster_A_nums(:)) || any(NodeL(nodeIndex) == cluster_B_nums(:))
                    
                    filledThroats = 0;
                    for j = 1:network.Nodes{nodeIndex}.connectionNumber
                        
                        if network.Links{network.Nodes{nodeIndex}.connectedLinks(j)}.occupancy == 'A'
                            
                            filledThroats = filledThroats + 1;
                        end
                    end
                    
                    if filledThroats ~= 0 &&  any(network.Nodes{nodeIndex}.imbThresholdPressure_PoreBodyFilling) && ...
                            network.Nodes{nodeIndex}.imbThresholdPressure_PoreBodyFilling >= Pc_imb
                        
                        network.Nodes{nodeIndex}.occupancy = 'A'; % make the pore wettingPhase type
                        network.Nodes{nodeIndex}.nonWettingLayerExistance()
                        network.Nodes{nodeIndex}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                        poreVolumeInjected = poreVolumeInjected + network.Nodes{nodeIndex}.volume;
                        numOfNodes_PoreBodyFilling = numOfNodes_PoreBodyFilling + 1;
                        network.Nodes{nodeIndex}.isInvaded = true;
                        e = e+1;
                        invaded_Element(e,5:8) = [nodeIndex, percList(network.numberOfLinks+nodeIndex), ...
                            network.Nodes{nodeIndex}.imbThresholdPressure_PoreBodyFilling, network.Nodes{nodeIndex}.imbThresholdPressure_SnapOff];
                        
                        for j=1:network.Nodes{nodeIndex}.connectionNumber
                            percList(network.Nodes{nodeIndex}.connectedLinks(j)) = -1000000;
                            if network.Links{network.Nodes{nodeIndex}.connectedLinks(j)}.occupancy == 'B'
                                percList(network.Nodes{nodeIndex}.connectedLinks(j))=...
                                    network.Links{network.Nodes{nodeIndex}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                            end
                        end
                        
                    elseif filledThroats == 0 && any(network.Nodes{nodeIndex}.imbThresholdPressure_SnapOff) && ...% if the node is non circular
                            network.Nodes{nodeIndex}.imbThresholdPressure_SnapOff >= Pc_imb
                        
                        network.Nodes{nodeIndex}.occupancy = 'A'; % make the pore wettingPhase type
                        network.Nodes{nodeIndex}.nonWettingLayerExistance()
                        network.Nodes{nodeIndex}.calculateThresholdPressureLayerCollapse(network.Pc_drain_max_Pa);
                        poreVolumeInjected = poreVolumeInjected + network.Nodes{nodeIndex}.volume;
                        numOfNodes_SnapOff = numOfNodes_SnapOff + 1;
                        network.Nodes{nodeIndex}.isInvaded = true;
                        e = e+1;
                        invaded_Element(e,5:8) = [nodeIndex, percList(network.numberOfLinks+nodeIndex), ...
                            network.Nodes{nodeIndex}.imbThresholdPressure_PoreBodyFilling, network.Nodes{nodeIndex}.imbThresholdPressure_SnapOff];
                        
                        for j=1:network.Nodes{nodeIndex}.connectionNumber
                            percList(network.Nodes{nodeIndex}.connectedLinks(j)) = -1000000;
                            if network.Links{network.Nodes{nodeIndex}.connectedLinks(j)}.occupancy == 'B'
                                percList(network.Nodes{nodeIndex}.connectedLinks(j))=...
                                    network.Links{network.Nodes{nodeIndex}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                            end
                        end
                    end
                end
            end
            
            percList(indexElement) = -1000000;
            deltaS = poreVolumeInjected /network.networkPoreVolume_m3 ;
            if deltaS > network.deltaS_input
                press = 0;
            end
            if max(percList)>= Pc_imb
                pressure = 1;
            else
                pressure = 0;
            end
        end
    end
    
    if Pc_imb < 0 % forced imbibition
        % Updating Pc collapse of the layers
        for ii = 1:network.numberOfNodes
            
            if any(network.Nodes{ii}.nonWettingLayerExist) && any(network.Nodes{ii}.imbThresholdPressure_LayerCollapse(1,:))...
                    && (any(NodeL(ii) == cluster_A_nums(:)) || any(NodeL(ii) == cluster_B_nums(:)) )
                
                % Updating Pc of layer collapse
                % Cheking layer collapse
                for jj = 1:4
                    if ~isnan(network.Nodes{ii}.imbThresholdPressure_LayerCollapse(1,j)) && ...
                            network.Nodes{ii}.imbThresholdPressure_LayerCollapse(1,j) > Pc_imb
                        
                        network.Nodes{ii}.nonWettingLayerExist(1,j) = nan;
                        
                        numOfNodes_LayerCollapse = numOfNodes_LayerCollapse + 1;
                    end
                end
            end
        end
        for ii = 1:network.numberOfLinks
            
            if any(network.Links{ii}.nonWettingLayerExist) && any(network.Links{ii}.imbThresholdPressure_LayerCollapse(1,:))...
                    && (any(LinkL(ii) == cluster_A_nums(:))|| any(LinkL(ii) == cluster_B_nums(:)) )
                
                % Updating Pc of layer collapse
                % Cheking layer collapse
                for jj = 1:4
                    
                    if ~isnan(network.Links{ii}.imbThresholdPressure_LayerCollapse(1,j)) && ...
                            network.Links{ii}.imbThresholdPressure_LayerCollapse(1,j) > Pc_imb
                        
                        network.Links{ii}.nonWettingLayerExist(1,j) = nan;
                        numOfLinks_LayerCollapse = numOfLinks_LayerCollapse + 1;
                    end
                end
            end
        end
    end
    
    if inv
        
        invaded = numOfLinks_SnapOff + numOfLinks_PistoneLike + ...
            numOfNodes_SnapOff + numOfNodes_PoreBodyFilling + numOfNodes_LayerCollapse;
        
        t = t+1;
        Pc_imb = Pc_min;
        % Updating saturations and conductances
        calculateConductance_and_Saturation_Imbibition(network, Pc_imb,NodeL, NodeL_W, LinkL, LinkL_W, cluster_A_nums, cluster_A, cluster_B_nums, cluster_B);
        
        % Relative Permeability Calculation
        if network.calculateRelativePermeability
            calculatePressureDistribution_TwoPhases(network, network.inletPressure_Pa, network.outletPressure_Pa);
            [Kr_w_imb, Kr_nw_imb] =...
                calculateRelativePermeability_Imbibition(network, network.outletPressure_Pa, LinkL, LinkL_W, cluster_A_nums, cluster_A);
            
    else
        Kr_w_imb =NaN;
        Kr_nw_imb=NaN;
        end
        network.ImbibitionData(t,:) = ...
            [network.wettingPhaseSaturation,Pc_imb,Kr_w_imb, Kr_nw_imb,invaded, ...
            numOfLinks_SnapOff,numOfLinks_PistoneLike, ...
            numOfLinks_LayerCollapse,numOfNodes_SnapOff, ...
            numOfNodes_PoreBodyFilling,numOfNodes_LayerCollapse,Pc_min];
        if network.flowVisualization
            network.IO.visualization(network,'SI',t);
        end
        [~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = clustering_nonWettingPhase(network);
        [~, NodeL_W, LinkL_W,cluster_A,cluster_B] = clustering_wettingPhase(network);
    end
    inv = false;
    if press ~= 0
        Pc_imb = Pc_imb - NoOfPc_interval;
    end
    
end

[~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = clustering_nonWettingPhase(network);
[~, NodeL_W, LinkL_W,cluster_A,cluster_B] = clustering_wettingPhase(network);

% Saturation calculation before invading

if  Pc_imb ~= 0
    Pc_imb = 0;
    t = t+1;
    
    invaded = numOfLinks_SnapOff + numOfLinks_PistoneLike + ...
        numOfNodes_SnapOff + numOfNodes_PoreBodyFilling + numOfNodes_LayerCollapse;
    
    % Updating saturations and conductances
    calculateConductance_and_Saturation_Imbibition(network, Pc_imb,NodeL, NodeL_W, LinkL, LinkL_W, cluster_A_nums, cluster_A, cluster_B_nums, cluster_B);
    
    % Relative Permeability Calculation
    if network.calculateRelativePermeability
        calculatePressureDistribution_TwoPhases(network, network.inletPressure_Pa, network.outletPressure_Pa);
        
        [Kr_w_imb, Kr_nw_imb] =...
            calculateRelativePermeability_Imbibition(network, network.outletPressure_Pa, LinkL, LinkL_W, cluster_A_nums, cluster_A);
    else
        
        Kr_w_imb =NaN;
        Kr_nw_imb=NaN;
    end
    network.ImbibitionData(t,:) = ...
        [network.wettingPhaseSaturation,Pc_imb,Kr_w_imb, Kr_nw_imb,invaded, ...
        numOfLinks_SnapOff,numOfLinks_PistoneLike, ...
        numOfLinks_LayerCollapse,numOfNodes_SnapOff, ...
        numOfNodes_PoreBodyFilling,numOfNodes_LayerCollapse,Pc_min];
    if network.flowVisualization
        network.IO.visualization(network,'SI',t);
    end
end

network.ImbibitionData = network.ImbibitionData(1:t,:);
rowNames = {'Sw', 'Pc (Pa)', 'Kr_w', 'Kr_nw', 'No. of invaded', ...
    'numOfLinks_SnapOff', 'numOfLinks_PistoneLike', 'numOfLinks_LayerCollapse',...
    'numOfNodes_SnapOff', 'numOfNodes_PoreBodyFilling', 'numOfNodes_LayerCollapse', 'Pc_min'};
network.ImbibitionData_table = array2table(network.ImbibitionData, 'VariableNames',rowNames);
network.sequence(1:network.numberOfLinks,1:9) = invaded_Element(1:network.numberOfLinks,1:9);
end