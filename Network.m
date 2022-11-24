classdef Network < handle & Fluids 
    
    properties
        Nodes
        Links
        xDimension
        yDimension
        zDimension
        
        numberOfNodes        
        numberOfLinks
        numOfInletLinks      
        numOfOutletLinks
        numOfSquareElements
        numOfCircularElements
        numOfTriangularElements
        numOfSquarePores
        numOfCircularPores
        numOfTriangularPores
        maxCoordinationNumber
        averageCoordinationNumber
        numOfIsolatedElements
        
        PSD 
        ThSD
        Porosity
        poreVolume
        networkVolume
        absolutePermeability
        absolutePermeability_m2
        
        totalFlowRate 
        velocity
        capillaryNumber
        pecletNumber
        
        Pc_drain_max        
        DrainageData 
        ImbibitionData   
        sequence
        waterSaturation
        thresholdPressure %control
        
        % Reactive        
        BreakThroughCurve_singlePhase
        
    end    
    methods
        %% Cunstructor function
        function obj = Network(fileName) 
            
            % Opening the files
            link_1_fileID = fopen(strcat(fileName, '_link1.dat'));
            obj.numberOfLinks = str2double(fgetl(link_1_fileID));
            link_2_fileID = fopen(strcat(fileName, '_link2.dat'));
            
            node_2_fileID = fopen(strcat(fileName, '_node2.dat'));            
            node_1_fileID = fopen(strcat(fileName, '_node1.dat'));
            temp = str2num(fgetl(node_1_fileID));
            obj.numberOfNodes = temp(1);
            
            % Network dimension
            obj.xDimension = temp(2);
            obj.yDimension = temp(3);
            obj.zDimension = temp(4);
            
            % Initializing Nodes and Links parameters
            obj.Nodes = cell(obj.numberOfNodes,1);
            obj.Links = cell(obj.numberOfLinks,1);
            
            % 
            for i = 1:obj.numberOfNodes
                node_1_values = str2num(fgetl(node_1_fileID));
                node_2_values = str2num(fgetl(node_2_fileID));
                obj.Nodes{i} = Node(node_1_values(1),... %pore index
                                    node_1_values(2),... % pore x coordinate
                                    node_1_values(3),... % pore y coordinate
                                    node_1_values(4),... % pore z coordinate
                                    node_1_values(5),... %pore connection number
                                    node_1_values(6:end),... % inlet-outlet status and connected link index
                                    node_2_values(2),... % pore volume
                                    node_2_values(3),... % pore radius  
                                    node_2_values(4),... % pore shape factor 
                                    node_2_values(5)); % pore clay volume               
            end        
            
            for i = 1:obj.numberOfLinks
               link_1_values = str2num(fgetl(link_1_fileID));
               link_2_values = str2num(fgetl(link_2_fileID));
               obj.Links{i} = Link(link_1_values(1),... %index 
                                    link_1_values(2),... %pore1Index,... 
                                    link_1_values(3),... %pore2Index,...
                                    link_1_values(4),... %radius,...
                                    link_1_values(5),... %shapeFactor,...
                                    link_1_values(6),... %length,...
                                    link_2_values(4),... %pore1Length,...
                                    link_2_values(5),... %pore2Length,...
                                    link_2_values(6),... %linkLength,...
                                    link_2_values(7),... %volume,...
                                    link_2_values(8)); %clayVolume
                                
                                
            end
            
            %closing the files
            fclose(link_1_fileID); fclose(link_2_fileID);
            fclose(node_1_fileID); fclose(node_2_fileID);    
            
        end
        
        %% Single Phase
        %Network Properties calculation
        function calculateNetworkProperties(obj, inletPressure, outletPressure)
            obj.ThSD = zeros(obj.numberOfLinks,1);
            obj.PSD = zeros(obj.numberOfNodes,1);
            CoordinationNumber = zeros(obj.numberOfNodes,1);
            obj.numOfInletLinks = 0;
            obj.numOfOutletLinks = 0; 
            obj.averageCoordinationNumber = 0;
            obj.numOfIsolatedElements = 0;
            obj.numOfTriangularElements = 0;
            obj.numOfCircularElements = 0;
            obj.numOfSquareElements = 0;
            obj.numOfTriangularPores = 0;
            obj.numOfCircularPores = 0;
            obj.numOfSquarePores = 0;
            nodesVolume = 0;
            linksVolume = 0;

            for ii = 1:obj.numberOfNodes
                obj.Nodes{ii}.calculateElementsProperties
                nodesVolume = nodesVolume + (obj.Nodes{ii}.volume); 
                CoordinationNumber(ii,1) = obj.Nodes{ii}.connectionNumber; 
                obj.PSD(ii,1) = 2 * obj.Nodes{ii}.radius; 
                %Isolated element
                if obj.Nodes{ii}.connectionNumber == 0
                    obj.numOfIsolatedElements = obj.numOfIsolatedElements + 1;
                end
                if strcmp(obj.Nodes{ii}.geometry , 'Circle')== 1
                    obj.numOfCircularPores = obj.numOfCircularPores+1;
                    obj.numOfCircularElements = obj.numOfCircularElements+1;
                elseif strcmp(obj.Nodes{ii}.geometry , 'Triangle')== 1
                    obj.numOfTriangularPores = obj.numOfTriangularPores+1;
                else
                    obj.numOfSquarePores = obj.numOfSquarePores+1;
                end
            end 

            for ii = 1:obj.numberOfLinks                 
                obj.Links{ii}.calculateElementsProperties
                linksVolume = linksVolume + (obj.Links{ii}.volume);
                obj.ThSD (ii,1)= 2 * obj.Links{ii}.radius;
                if obj.Links{ii}.isInlet
                    obj.numOfInletLinks = obj.numOfInletLinks + 1;
                elseif obj.Links{ii}.isOutlet
                    obj.numOfOutletLinks = obj.numOfOutletLinks+1;                 
                end             
                if strcmp(obj.Links{ii}.geometry , 'Circle')== 1
                    obj.numOfCircularElements = obj.numOfCircularElements+1;
                elseif strcmp(obj.Links{ii}.geometry , 'Triangle')== 1
                    obj.numOfTriangularElements = obj.numOfTriangularElements+1;
                else
                    obj.numOfSquareElements = obj.numOfSquareElements+1;
                end
            end 

            obj.averageCoordinationNumber = sum(CoordinationNumber)/obj.numberOfNodes;
            obj.maxCoordinationNumber = max(CoordinationNumber); 
            obj.networkVolume = obj.xDimension * obj.yDimension * obj.zDimension;
            obj.poreVolume = linksVolume + nodesVolume;
            obj.Porosity = obj.poreVolume / (obj.xDimension * obj.yDimension * obj.zDimension);      
            calculateAbsolutePermeability(obj, inletPressure, outletPressure);
        end 
        
        % Pressure distribution calculation of single phase flow       
        function pressureDistribution_singlePhaseFlow (obj, inletPressure, outletPressure)
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            B = zeros(obj.numberOfNodes, 1);
     
            for ii = 1:obj.numberOfLinks
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;

                % if the link is connected to inlet (index of node 1 is -1 which does not exist) 
                if obj.Links{ii}.isInlet
                    obj.Links{ii}.nodeLinkSystemConductanceSinglePhase = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductanceSinglePhase) +...
                        ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductanceSinglePhase)))^-1;
                    
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + obj.Links{ii}.nodeLinkSystemConductanceSinglePhase;
                    B(node2Index) = obj.Links{ii}.nodeLinkSystemConductanceSinglePhase * inletPressure;
                    
                % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                elseif obj.Links{ii}.isOutlet
                     obj.Links{ii}.nodeLinkSystemConductanceSinglePhase = ( (obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductanceSinglePhase) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductanceSinglePhase)))^-1;
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + obj.Links{ii}.nodeLinkSystemConductanceSinglePhase;
                    B(node1Index) = obj.Links{ii}.nodeLinkSystemConductanceSinglePhase * outletPressure;
                    
                %if the link is neither inlet nor outlet    
                else
                    obj.Links{ii}.nodeLinkSystemConductanceSinglePhase = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductanceSinglePhase) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductanceSinglePhase) +...
                        (obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductanceSinglePhase)))^-1;   
                
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + obj.Links{ii}.nodeLinkSystemConductanceSinglePhase;
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + obj.Links{ii}.nodeLinkSystemConductanceSinglePhase;
                    Factor(node1Index, node2Index) = Factor(node1Index, node2Index) - obj.Links{ii}.nodeLinkSystemConductanceSinglePhase;
                    Factor(node2Index, node1Index) = Factor(node2Index, node1Index) - obj.Links{ii}.nodeLinkSystemConductanceSinglePhase;
                   
                end     
            end
            
            % using GMRES method to solve the pressure distribution             
            nodesPressure = gmres(Factor, B,[], 1e-10, obj.numberOfNodes);
%             nodesPressure = pcg(Factor, B, 1e-7, 1000);

            %assign the pressure values to each node
            for ii = 1:obj.numberOfNodes
                if nodesPressure(ii) > inletPressure
                    obj.Nodes{ii}.waterPressure = inletPressure; 
                elseif nodesPressure(ii) < outletPressure
                    obj.Nodes{ii}.waterPressure = outletPressure; 
                else
                    obj.Nodes{ii}.waterPressure = nodesPressure(ii); 
                end                
            end
            
            %assign pressure values to links, since the surface where
            %flowrate is calculated through might pass through the links
            for ii = 1:obj.numberOfLinks
                if obj.Links{ii}.isInlet
                    obj.Links{ii}.waterPressure =...
                        (1+obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure)/2;
                elseif obj.Links{ii}.isOutlet
                    obj.Links{ii}.waterPressure =...
                        obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure/2;                    
                else
                    obj.Links{ii}.waterPressure =...
                        (obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure + ...
                        obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure) / 2;
                end
            end                   
           
        end  
               
        % Flow rate calculation for each phase in the netwrok
        function calculateFlowRate(obj, inletPressure, outletPressure)
            % Fluid = water
            pressureDistribution_singlePhaseFlow(obj, inletPressure,outletPressure);  
            
            % calculate flow rate in Inlet_Links
            obj.totalFlowRate = 0;  
            for ii = 1:obj.numberOfLinks 
                
                node1Index = obj.Links{ii}.pore1Index;
                
                if obj.Links{ii}.isOutlet
                    
                    % calculate the flow rate of the fluid
                    obj.totalFlowRate = obj.totalFlowRate + ...
                        obj.Links{ii}.nodeLinkSystemConductanceSinglePhase * ...
                        (obj.Nodes{node1Index}.waterPressure - outletPressure);  
                end
            end
            
            % calculate velocity through the network 
            obj.velocity = obj.totalFlowRate/(obj.yDimension * obj.zDimension); 
            
            % for quasi-static, capillaryNumber must be less than 10e-4
            obj.capillaryNumber = obj.waterViscosity * obj.velocity/ obj.sig_ow;  
        end 
        
        % AbsolutePermeability
        % Piri: eliminate boundary condition        
        function calculateAbsolutePermeability_box(obj, inletPressure, outletPressure) 
            % Fluid = water
            pressureDistribution_singlePhaseFlow(obj, inletPressure,outletPressure);  
            
            % calculate flow rate in Inlet_Links
            obj.totalFlowRate = 0;  
            sumPA_a = 0;
            sumA_a = 0;
            sumPA_b = 0;
            sumA_b = 0;
            x_coor = zeros(obj.numberOfNodes,1);
            
            for ii = 1:obj.numberOfNodes
                x_coor(ii,1) = obj.Nodes{ii}.x_coordinate;
            end
            x_outlet = max(x_coor);
            x_inlet = min(x_coor);
            intervalx = (x_outlet - x_inlet)/20; % Should check > linkLength to avoid index -1 & 0 in nodes
            a = x_inlet + intervalx;
            b = x_outlet - intervalx;  
                for i = 1:obj.numberOfNodes 
                    if obj.Nodes{i}.x_coordinate == a
                        sumA_a = sumA_a + obj.Nodes{i}.area ;
                        sumPA_a = sumPA_a + obj.Nodes{i}.area * obj.Nodes{i}.waterPressure;
                    elseif obj.Nodes{i}.x_coordinate == b
                        sumA_b = sumA_b + obj.Nodes{i}.area ;
                        sumPA_b = sumPA_b + obj.Nodes{i}.area * obj.Nodes{i}.waterPressure;
                    end
                end
                for i = 1:obj.numberOfLinks
                    if ~obj.Links{i}.isOutlet && ~obj.Links{i}.isInlet
                        if obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate > a && obj.Nodes{obj.Links{i}.pore1Index}.x_coordinate < a || ...
                             obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate < a && obj.Nodes{obj.Links{i}.pore1Index}.x_coordinate > a
                            sumA_a = sumA_a + obj.Links{i}.area ;
                            sumPA_a = sumPA_a + obj.Links{i}.area * obj.Links{i}.waterPressure;
                        elseif obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate > b && obj.Nodes{obj.Links{i}.pore1Index}.x_coordinate < b || ...
                                obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate < b && obj.Nodes{obj.Links{i}.pore1Index}.x_coordinate > b
                            sumA_b = sumA_b + obj.Links{i}.area ;
                            sumPA_b = sumPA_b + obj.Links{i}.area * obj.Links{i}.waterPressure;
                        end
                    elseif obj.Links{i}.isOutlet
                    
                    % calculate the flow rate of the fluid
                    obj.totalFlowRate = obj.totalFlowRate + ...
                        obj.Links{i}.nodeLinkSystemConductanceSinglePhase * ...
                        (obj.Nodes{obj.Links{i}.pore1Index}.waterPressure - outletPressure); 
                    end
                end  
            
                Pa = sumPA_a / sumA_a;
                Pb = sumPA_b / sumA_b; 
            
            % calculate velocity through the network 
            obj.velocity = obj.totalFlowRate/(obj.yDimension * obj.zDimension); 
            
            % for quasi-static, capillaryNumber must be less than 10e-4
            obj.capillaryNumber = obj.waterViscosity * obj.velocity/ obj.sig_ow;        
            unitConvertor = 1/0.987*10^15; % unit conversion from m2 to miliDarcy
            obj.absolutePermeability = unitConvertor * obj.velocity * (b-a) * obj.waterViscosity/ (Pa -Pb );
            format longE
            obj.absolutePermeability_m2 = obj.velocity * (b-a) * obj.waterViscosity/ (Pa -Pb );
        end  
        
        function calculateAbsolutePermeability(obj, inletPressure, outletPressure) 
            calculateFlowRate(obj, inletPressure, outletPressure);                        
            unitConvertor = 1/0.987*10^15; % unit conversion from m2 to miliDarcy
            obj.absolutePermeability = unitConvertor * obj.velocity * obj.xDimension * obj.waterViscosity/ (inletPressure -outletPressure );
            format longE
            obj.absolutePermeability_m2 = obj.velocity * obj.xDimension * obj.waterViscosity/ (inletPressure -outletPressure );
        end   
        
        % Single phase Reactive transport   
        % Pressure distribution calculation in Cylindrical pore_Single-Phase       
        function pressureDistribution_Cylindrical (obj, inletPressure, outletPressure) 
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            B = zeros(obj.numberOfNodes, 1);
     
            for ii = 1:obj.numberOfLinks
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index; 
                
                % Calculate conductance based on Raoof thesis
                obj.Links{ii}.cylindricalConductance = pi * obj.Links{ii}.radius ^ 4 /...
                    (8 * obj.waterViscosity * obj.Links{ii}.length); 

                % if the link is connected to inlet (index of node 1 is -1 which does not exist) 
                if obj.Links{ii}.isInlet
                     
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + ...
                        obj.Links{ii}.cylindricalConductance;
                    B(node2Index) = obj.Links{ii}.cylindricalConductance * inletPressure;
           
                % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                elseif obj.Links{ii}.isOutlet 
                    
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + ...
                        obj.Links{ii}.cylindricalConductance;
                    B(node1Index) = obj.Links{ii}.cylindricalConductance * outletPressure; 
                    
                %if the link is neither inlet nor outlet    
                else 
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + ...
                        obj.Links{ii}.cylindricalConductance;
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + ...
                        obj.Links{ii}.cylindricalConductance;
                    Factor(node1Index, node2Index) = Factor(node1Index, node2Index) - ...
                        obj.Links{ii}.cylindricalConductance;
                    Factor(node2Index, node1Index) = Factor(node2Index, node1Index) - ...
                        obj.Links{ii}.cylindricalConductance;                   
                end     
            end
            
            % using GMRES method to solve the pressure distribution 
            nodesPressure = gmres(Factor, B,[], 1e-7, 1000);
            
            %assign the pressure values to each node 
            for ii = 1:obj.numberOfNodes
                if nodesPressure(ii) > inletPressure
                    obj.Nodes{ii}.waterPressure = inletPressure; 
                elseif nodesPressure(ii) < outletPressure
                    obj.Nodes{ii}.waterPressure = outletPressure; 
                else
                    obj.Nodes{ii}.waterPressure = nodesPressure(ii); 
                end                 
            end
            
            %assign pressure values to links, since the surface where
            %flowrate is calculated through might pass through the links
            for ii = 1:obj.numberOfLinks
                if obj.Links{ii}.isInlet
                    obj.Links{ii}.waterPressure =...
                        (1+obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure)/2;
                elseif obj.Links{ii}.isOutlet
                    obj.Links{ii}.waterPressure =...
                        obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure/2;                    
                else
                    obj.Links{ii}.waterPressure =...
                        (obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure + ...
                        obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure) / 2;
                end
            end                     
            
        end    
        % Single phase Reactive transport_Raoof_2017 
        function calculateReactiveTransport_SinglePhaseDiffusion(obj, inletPressure, outletPressure, soluteConcentration, poreVolumeInjected) 
            
            % calculate pressure distribution
            pressureDistribution_Cylindrical (obj, inletPressure, outletPressure); 
            
            residenceTime_link = zeros(obj.numberOfLinks,1);
            flowRate_link = zeros(obj.numberOfLinks,1);
            obj.totalFlowRate = 0;            
            effectiveDiffusion = 10^ (-9);         
            diffusion_link = zeros(obj.numberOfLinks,1);
            
            % calculate flowrate of links residence time                          
            for ii = 1:obj.numberOfLinks                
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;
                
                if ~obj.Links{ii}.isInlet && ~obj.Links{ii}.isOutlet      
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);  
                    
                elseif obj.Links{ii}.isInlet 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(inletPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);  
                    obj.totalFlowRate = obj.totalFlowRate + flowRate_link(ii);
                else 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        outletPressure); 
                end 
                residenceTime_link(ii) = obj.Links{ii}.volume/flowRate_link(ii);
                diffusion_link (ii) = effectiveDiffusion * obj.Links{ii}.area / obj.Links{ii}.linkLength;
            end
            timeStep = min(nonzeros(residenceTime_link)); 
%             timeStep = timeStep * 10;
            flowRate_node = zeros(obj.numberOfNodes,1); 
            diffusion_node = zeros(obj.numberOfNodes,1);  
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            Known = zeros(obj.numberOfNodes, 1);
            
            % calculation of 3 Unknowns (concentration of nodes & links) in each timeStep 
             
            t = 0; 
            time = 0;
            simulationTime = obj.poreVolume / obj.totalFlowRate/15;
            injectionTime = poreVolumeInjected / obj.totalFlowRate;
            timePlot = zeros(round(simulationTime/timeStep)+1 ,1);
            flux_averagedConcentration = zeros(round(simulationTime/timeStep)+1 ,1);
            obj.BreakThroughCurve_singlePhase = zeros(round(simulationTime/timeStep)+1 ,2);
            soluteConcentration1 = soluteConcentration;
            % Plot & Animation
            h = animatedline;
            h.Color = 'b';
            h.LineStyle = '-';
            h.LineWidth = 2; 
            axis([0 simulationTime 0 1])
            title('Break Through Curve')
            xlabel('Time(s)')            
            ylabel('DimensionlessConcentration(-)') 
            fig = figure; 
            figure('name','BTC')
             
            while time < simulationTime
                
                if time > injectionTime*0.6
                    soluteConcentration = 0;
                end
                
                t = t+1;
                time = time + timeStep; 
                timePlot(t) = time;
                sumOfConcentration = 0;
                sumOfFlowRate = 0;
            % calculate concentration of nodes & links based on eq. 8 & 9
            for i = 1:obj.numberOfNodes  
                    
                for j = 1:obj.Nodes{i}.connectionNumber 
                    
                    connectedLinkIndex = obj.Nodes{i}.connectedLinks(j);                    
                    connectedNodeIndex = obj.Nodes{i}.connectedNodes(j);                 
                    
                    diffusion_node(i) = diffusion_node(i) + diffusion_link(connectedLinkIndex);  
                    Known(i,1) = Known(i,1) + obj.Links{i}.concentration(t)*...
                                timeStep / obj.Nodes{i}.volume *  diffusion_link(connectedLinkIndex);
                            
                    % determine link flowing into this node
                    if connectedNodeIndex ~= 0 && connectedNodeIndex ~= -1                        
                        
                        if obj.Nodes{connectedNodeIndex}.waterPressure > obj.Nodes{i}.waterPressure                             
                            flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex);  
                            Known(i,1) = Known(i,1) + obj.Links{i}.concentration(t)*...
                                timeStep / obj.Nodes{i}.volume *  flowRate_link(connectedLinkIndex);
                        end
                            
                    elseif connectedNodeIndex == -1                                                
                            flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex);
                            Known(i,1) = Known(i,1) + obj.Links{i}.concentration(t)*...
                                timeStep / obj.Nodes{i}.volume *  flowRate_link(connectedLinkIndex);
                    end
                end
                
                Factor(i, i) = 1 + timeStep * (flowRate_node(i)+  diffusion_node(i)) / obj.Nodes{i}.volume  ;
                Known(i,1) = Known(i,1) + obj.Nodes{i}.concentration(t);
            end
            
            nodesConcentration_new = gmres (Factor, Known,[], 1e-10, obj.numberOfNodes);
                        
            % asign new concentration of nodes
            for i = 1:obj.numberOfNodes
                if nodesConcentration_new(i) > soluteConcentration
                obj.Nodes{i}.concentration(t+1) = soluteConcentration;
                else
                obj.Nodes{i}.concentration(t+1) = nodesConcentration_new(i);
                end
            end
            
            % calculate new concentration of links
            for i = 1:obj.numberOfLinks 
                
                node1Index = obj.Links{i}.pore1Index;
                node2Index = obj.Links{i}.pore2Index;
                
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet
                    if obj.Nodes{node1Index}.waterPressure > obj.Nodes{node2Index}.waterPressure
                    obj.Links{i}.concentration(t+1) = (obj.Links{i}.concentration(t)+...
                         timeStep  / obj.Links{i}.volume * (flowRate_link(i)+ ...
                         diffusion_link(i))*obj.Nodes{node1Index}.concentration(t+1)+...
                         timeStep  / obj.Links{i}.volume *diffusion_link(i) * obj.Nodes{node2Index}.concentration(t+1))/...
                         (1 + timeStep / obj.Links{i}.volume * (flowRate_link(i)+diffusion_link(i)));
                    else 
                        obj.Links{i}.concentration(t+1) = (obj.Links{i}.concentration(t)+...
                         timeStep  / obj.Links{i}.volume * (flowRate_link(i)+ ...
                         diffusion_link(i))*obj.Nodes{node2Index}.concentration(t+1)+...
                         timeStep  / obj.Links{i}.volume *diffusion_link(i) * obj.Nodes{node1Index}.concentration(t+1))/...
                         (1 + timeStep / obj.Links{i}.volume * (flowRate_link(i)+diffusion_link(i)));
                    end 
                elseif obj.Links{i}.isInlet
                    obj.Links{i}.concentration(t+1) = (obj.Links{i}.concentration(t)+...
                         timeStep  / obj.Links{i}.volume * (flowRate_link(i)+ ...
                         diffusion_link(i))* soluteConcentration +...
                         timeStep  / obj.Links{i}.volume *diffusion_link(i) * obj.Nodes{node2Index}.concentration(t+1))/...
                         (1 + timeStep / obj.Links{i}.volume * (flowRate_link(i)+diffusion_link(i)));
                elseif obj.Links{i}.isOutlet 
                     obj.Links{i}.concentration(t+1) = (obj.Links{i}.concentration(t)+...
                         timeStep  / obj.Links{i}.volume * (flowRate_link(i)+ ...
                         diffusion_link(i))*obj.Nodes{node1Index}.concentration(t+1))/...
                         (1 + timeStep / obj.Links{i}.volume * (flowRate_link(i)+diffusion_link(i))); 
                    sumOfConcentration = sumOfConcentration + ...
                        obj.Links{i}.concentration(t)*flowRate_link(i);
                    sumOfFlowRate = sumOfFlowRate + flowRate_link(i);
                end 
            end
            % calculate BreakThroughCurve at outlet of network 
            flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration1;
%             nameFile = vtkWriter(obj,'Diff',t);
             fileName = obj.vtkWriter_glyph('Diff',t);
%             movefile nameFile D:\Uni\PUT\Dissertation\PNM_Code\MatlabPNM\Ver.2\seq\Diff
            % Plot & Animation
            addpoints(h,timePlot(t),flux_averagedConcentration(t));  
            % GIF 
%             plot(timePlot(1:t),flux_averagedConcentration(1:t), 'b','LineWidth',2);            
%             axis([0 simulationTime 0 1])
%             title('Break Through Curve')
%             xlabel('Time(s)')            
%             ylabel('DimensionlessConcentration(-)') 
            drawnow 
            frame = getframe(fig);
            im{t} = frame2im(frame);
            end
            close;
            figure;                                 
            axis([0 simulationTime 0 1])
            for idx = 1:length(timePlot)
                subplot(15,15,idx)
                imshow(im{idx});
            end
            filename = 'BTC.gif'; % Specify the output file name
            for idx = 1:length(timePlot)
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
                end
            end
            obj.BreakThroughCurve_singlePhase(:,1) = timePlot;
            obj.BreakThroughCurve_singlePhase(:,2) = flux_averagedConcentration; 
            
            % Plot & Animation 
%             for k = 1:length(timePlot) 
%                 addpoints(h,timePlot(k),flux_averagedConcentration(k));
%                 drawnow
%                 pause(0.05); 
%             end
%             plot(timePlot,flux_averagedConcentration,'*');
        end 
        % Single phase Reactive transport_Raoof_2010  
        function calculateReactiveTransport_SinglePhaseDesorption(obj, inletPressure, outletPressure, soluteConcentration, poreVolumeInjected)
            
            % calculate pressure distribution
            pressureDistribution_Cylindrical (obj, inletPressure, outletPressure); 
            
            % mass transfer coefficient: alpha & distribution coefficient: K_d
            alpha = 45/24/3600*ones(obj.numberOfLinks,1);
            K_d = 10/45*ones(obj.numberOfLinks,1);
            obj.capillaryNumber = 1;
            obj.pecletNumber = 1; 
            
            residenceTime_link = zeros(obj.numberOfLinks,1);
            flowRate_link = zeros(obj.numberOfLinks,1);
            obj.totalFlowRate = 0;            
            
            % calculate flow rate of links residence time based on eq. 1 & 10                         
            for ii = 1:obj.numberOfLinks                
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;
                
                if ~obj.Links{ii}.isInlet && ~obj.Links{ii}.isOutlet      
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);  
                elseif obj.Links{ii}.isInlet 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(inletPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);   
                else 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        outletPressure);                     
                    obj.totalFlowRate = obj.totalFlowRate + flowRate_link(ii);
                end  
                residenceTime_link(ii) = obj.Links{ii}.volume/flowRate_link(ii);
            end
            
            % Set time step to avoid numerical dispersion
            timeStep = min(nonzeros(residenceTime_link));
            timeStep = timeStep * 40;             
            % for perfect mixing must be less than 1: diffusion is dominant
            % rather than advection
            diffusionCoefficent = 1e-9;
            obj.pecletNumber = obj.xDimension * obj.velocity / diffusionCoefficent;
            
            
            flowRate_node = zeros(obj.numberOfNodes,1);
            B = zeros(obj.numberOfLinks,1);            
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            Known = zeros(obj.numberOfNodes, 1);
            
            % calculation of 3 Unknowns (concentration of nodes & 
            % concentration & adsorbedConcentration of links) in each timeStep 
            
            t = 0; 
            time = 0;                
            simulationTime = obj.poreVolume / obj.totalFlowRate*2;
            injectionTime = poreVolumeInjected / obj.totalFlowRate;
            timePlot = zeros(round(simulationTime/timeStep)+1 ,1);
            flux_averagedConcentration = zeros(round(simulationTime/timeStep)+1 ,1);
            obj.BreakThroughCurve_singlePhase = zeros(round(simulationTime/timeStep)+1 ,2);
            soluteConcentration1 = soluteConcentration;
            
            %Plot & Animation
            h = animatedline;
            h.Color = 'b';
            h.LineStyle = '-';
            h.LineWidth = 2; 
%             axis([0 simulationTime 0 1])
            title('Break Through Curve')
            xlabel('Time(s)')            
            ylabel('DimensionlessConcentration(-)') 
            fig = figure;
            
            while time < simulationTime 
                
                if time > injectionTime
                    soluteConcentration = 0;
                end
                t = t+1
                time = time + timeStep;
                timePlot(t) = time;
                sumOfConcentration = 0;
                sumOfFlowRate = 0;
                
            % calculate concentration of nodes: based on eq.7
            for i = 1:obj.numberOfNodes 
                
                I = timeStep / obj.Nodes{i}.volume;
                for j = 1:obj.Nodes{i}.connectionNumber 
                    
                    connectedLinkIndex = obj.Nodes{i}.connectedLinks(j);
                    connectedNodeIndex = obj.Nodes{i}.connectedNodes(j);
                    adsorbedConcentration = obj.Links{connectedLinkIndex}.adsorbedConcentration(t);
                    linksConcentration = obj.Links{connectedLinkIndex}.concentration(t);
                    
                    B(connectedLinkIndex) = 1 + ...
                        (flowRate_link(connectedLinkIndex) * timeStep /...
                        obj.Links{connectedLinkIndex}.volume) + ...
                        (timeStep * alpha(connectedLinkIndex) * K_d(connectedLinkIndex)) - ...
                        (timeStep^2 * (alpha(connectedLinkIndex))^2 * K_d(connectedLinkIndex))/...
                        (1 + timeStep * alpha(connectedLinkIndex));
                    F = 1 / B(connectedLinkIndex) * ...
                        (flowRate_link(connectedLinkIndex))^2 * timeStep / obj.Links{connectedLinkIndex}.volume;
                    G = flowRate_link(connectedLinkIndex) /  B(connectedLinkIndex);
                    H = flowRate_link(connectedLinkIndex) * timeStep * alpha(connectedLinkIndex)/...
                        B(connectedLinkIndex)/(1 + timeStep * alpha(connectedLinkIndex));
                                
                    if connectedNodeIndex ~= 0 && connectedNodeIndex ~= -1
                        
                        % determine link flowing into this node
                        if obj.Nodes{connectedNodeIndex}.waterPressure > obj.Nodes{i}.waterPressure 
                            
                            flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex);       
                            
                            Factor(i, connectedNodeIndex) = -I * F;                    
                            Known(i,1) = Known(i,1) + G * linksConcentration + H * adsorbedConcentration;
                        end
                    elseif connectedNodeIndex == -1
                        
                        flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex); 
                        
                        Known(i,1) = Known(i,1) + ...
                                G * linksConcentration + H * adsorbedConcentration + ...
                                F * soluteConcentration;                    
                    end
                end
                Factor(i, i) = 1+ timeStep * flowRate_node(i) / obj.Nodes{i}.volume  ;
                Known(i,1) = obj.Nodes{i}.concentration(t) + I * Known(i,1);  
            end  
            
            nodesConcentration_new = gmres(Factor, Known,[], 1e-10, obj.numberOfNodes); 
            
            % asign new concentration of nodes
            for i = 1:obj.numberOfNodes
%                 if nodesConcentration_new(i) > soluteConcentration
%                 obj.Nodes{i}.concentration(t+1) = soluteConcentration;
%                 elseif nodesConcentration_new(i) < 0
%                  obj.Nodes{i}.concentration(t+1) = 0;
%                 else
                obj.Nodes{i}.concentration(t+1) = nodesConcentration_new(i);
%                 end  
            end 
           
            % calculate new concentration & adsorbedConcentration of links:
            % based on eq.3&4
            for i = 1:obj.numberOfLinks 
                
                node1Index = obj.Links{i}.pore1Index;
                node2Index = obj.Links{i}.pore2Index;
                
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet
                    if obj.Nodes{node1Index}.waterPressure > obj.Nodes{node2Index}.waterPressure
                        upstreamNode = node1Index;
                    else
                        upstreamNode = node2Index;
                    end
                    obj.Links{i}.concentration(t+1) = 1/ B(i)*(obj.Links{i}.concentration(t) + ...
                        timeStep * alpha(i) * obj.Links{i}.adsorbedConcentration(t) / ...
                        (1 + alpha(i) * timeStep) + ...
                        flowRate_link(i) * timeStep * obj.Nodes{upstreamNode}.concentration(t+1)/...
                        obj.Links{i}.volume);                      
                elseif obj.Links{i}.isInlet
                    obj.Links{i}.concentration(t+1) = 1/ B(i)*(obj.Links{i}.concentration(t)+ ...
                        timeStep * alpha(i) * obj.Links{i}.adsorbedConcentration(t) / ...
                        (1 + alpha(i) * timeStep) + ...
                        flowRate_link(i) * timeStep * soluteConcentration/obj.Links{i}.volume);
                else
                    obj.Links{i}.concentration(t+1) = 1/ B(i)*(obj.Links{i}.concentration(t) + ...
                        timeStep * alpha(i) * obj.Links{i}.adsorbedConcentration(t) / ...
                        (B(i)*(1 + alpha(i) * timeStep)) + ...
                        flowRate_link(i) * timeStep * obj.Nodes{node1Index}.concentration(t+1)/...
                        obj.Links{i}.volume);
                    
                    % calculation for BreakThroughCurve at outlet of network  
                    sumOfConcentration = sumOfConcentration + ...
                        obj.Links{i}.concentration(t+1)*flowRate_link(i);
                    sumOfFlowRate = sumOfFlowRate + flowRate_link(i);
                end
                obj.Links{i}.adsorbedConcentration(t+1) = (alpha(i) * timeStep * K_d(i) * ...
                    obj.Links{i}.concentration(t+1)+ obj.Links{i}.adsorbedConcentration(t))/...
                    (1 + alpha(i) * timeStep);
            end
            % calculate BreakThroughCurve at outlet of network
            flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration1; 
            
            % Plot & Animation
            addpoints(h,timePlot(t),flux_averagedConcentration(t)); 
            drawnow
            % GIF 
            plot(timePlot(1:t),flux_averagedConcentration(1:t), 'b','LineWidth',2);            
%             axis([0 simulationTime 0 1])
            title('Break Through Curve')
            xlabel('Time(s)')            
            ylabel('DimensionlessConcentration(-)')  
            frame = getframe(fig);
            im{t} = frame2im(frame);
            end
            close;
%             figure;                                 
%             axis([0 simulationTime 0 1])
            for idx = 1:length(timePlot)
                subplot(15,15,idx)
                imshow(im{idx});
            end
            filename = 'BTC.gif'; % Specify the output file name
            for idx = 1:length(timePlot)
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
                end
            end
            
%             obj.BreakThroughCurve_singlePhase(:,1) = timePlot;
%             obj.BreakThroughCurve_singlePhase(:,2) = flux_averagedConcentration;
%             plot(timePlot, flux_averagedConcentration,'-'); 
%             title('Break Through Curve')
%             xlabel('Time(s)')            
%             ylabel('DimensionlessConcentration(-)')
%             ylim([0 1]
        end
           
        
        %% Two Phases 
        % Contact Angle
        function contactAngleDistribution(obj) 
            
          %check & sort pores based on radius 
          [~, pIndex] = sort(obj.PSD);
          volFrac = 0.85;  
          n = round(volFrac * obj.numberOfNodes); m = obj.numberOfNodes - n;      
          min = 45 * pi /180; max = 55 * pi /180; largerInterval = (55-45)/(n-1)* pi /180;
          advCnntAng1 = min:largerInterval:max;
          min = 55 * pi /180; max = 75 * pi /180; largerInterval = (75-55)/(m-1)* pi /180;
          advCnntAng2 = min:largerInterval:max;
          
          for i = 1:obj.numberOfNodes 
              if i <= n
                obj.Nodes{pIndex(i)}.newContactAngle = advCnntAng1(i);
              else
                obj.Nodes{pIndex(i)}.newContactAngle = advCnntAng2(i-n);
              end                  
          end   
          for i = 1:obj.numberOfLinks          
              if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet
                obj.Links{i}.newContactAngle = ...
                    (obj.Nodes{obj.Links{i}.pore1Index}.newContactAngle +...
                    obj.Nodes{obj.Links{i}.pore2Index}.newContactAngle)/2;
              elseif obj.Links{i}.isInlet
                  obj.Links{i}.newContactAngle = ...
                    obj.Nodes{obj.Links{i}.pore2Index}.newContactAngle;
              else
                  obj.Links{i}.newContactAngle = ...
                    obj.Nodes{obj.Links{i}.pore1Index}.newContactAngle;
              end
          end            
        end
        
        % Conductance 2-Phase Flow and Water Saturation
        % Piri: eliminate boundary condition 
        function calculateConductance_and_Saturation_Drainage_box (obj,Pc)
            
            Pc = abs(Pc);     
            waterVolume = 0;
            vol = 0;             
            x_coor = zeros(obj.numberOfNodes,1);
            
            for ii = 1:obj.numberOfNodes
                x_coor(ii,1) = obj.Nodes{ii}.x_coordinate;
            end
            x_outlet = max(x_coor);
            x_inlet = min(x_coor);
            intervalx = (x_outlet - x_inlet)/20;
            a = x_inlet + intervalx;
            b = x_outlet - intervalx;  
                for i = 1:obj.numberOfNodes
                    obj.Nodes{i}.calculateConductance_Drainage(Pc);
                    if obj.Nodes{i}.x_coordinate >= a && obj.Nodes{i}.x_coordinate < b
                            waterVolume = waterVolume + (obj.Nodes{i}.waterCrossSectionArea )...
                                / obj.Nodes{i}.area *obj.Nodes{i}.volume + obj.Nodes{i}.clayVolume;
                            vol = vol + obj.Nodes{i}.volume + obj.Nodes{i}.clayVolume;
                    end
                end
                for i = 1:obj.numberOfLinks
                    obj.Links{i}.calculateConductance_Drainage(Pc); 
                    if ~obj.Links{i}.isOutlet
                        if obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate >= a && obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate < b
                            waterVolume = waterVolume + (obj.Links{i}.waterCrossSectionArea )...
                                / obj.Links{i}.area * obj.Links{i}.volume + obj.Links{i}.clayVolume;
                            vol = vol + obj.Links{i}.volume + obj.Links{i}.clayVolume;
                        end
                    end
                end  
            obj.waterSaturation = waterVolume / vol;    
        end
        
        function calculateConductance_and_Saturation_Drainage(obj, Pc)  
            Pc = abs(Pc);     
            waterVolume = 0;
            vol = 0; 
            
            for i = 1:obj.numberOfNodes 
                    obj.Nodes{i}.calculateConductance_Drainage(Pc); 
                    
                % Water Saturation Calculation 
                if ~obj.Nodes{i}.isInlet && ~obj.Nodes{i}.isOutlet 
                    waterVolume = waterVolume + (obj.Nodes{i}.waterCrossSectionArea )...
                        / obj.Nodes{i}.area *obj.Nodes{i}.volume + obj.Nodes{i}.clayVolume;                
                    vol = vol + obj.Nodes{i}.volume + obj.Nodes{i}.clayVolume;
                end
            end
            for i = 1:obj.numberOfLinks 
                    obj.Links{i}.calculateConductance_Drainage(Pc);     
                 
                % Water Saturation Calculation                 
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet 
                     waterVolume = waterVolume + (obj.Links{i}.waterCrossSectionArea )...
                         / obj.Links{i}.area * obj.Links{i}.volume + obj.Links{i}.clayVolume;
                     vol = vol + obj.Links{i}.volume + obj.Links{i}.clayVolume;
                end
            end             
            obj.waterSaturation = waterVolume / vol;          
        end        
        function calculateConductance_and_Saturation_Imb(obj, Pc, NodeL, NodeL_W, LinkL, LinkL_W, cluster_A_nums, cluster_A, cluster_B_nums, cluster_B)  
            Pc = abs(Pc);     
            waterVolume = 0;   
            vol = 0; 
            waterArea = zeros(obj.numberOfLinks,4);
            
            for i = 1:obj.numberOfNodes 
                    if (any(NodeL(i) == cluster_A_nums(:)) ||  any(NodeL(i) == cluster_B_nums(:)))
                        obj.Nodes{i}.calculateConductance_Imbibition(obj, Pc);   
                        
                    else
                        if isnan(obj.Nodes{i}.imbPressureTrapped) && obj.Nodes{i}.occupancy == 'B'
                            obj.Nodes{i}.imbPressureTrapped = max(Pc, obj.Nodes{i}.imbThresholdPressure_SnapOff);  
                        end 
                        if any(obj.Nodes{i}.imbPressureTrapped) 
                            obj.Nodes{i}.calculateConductance_Imbibition(obj, obj.Nodes{i}.imbPressureTrapped); 
                        else%if (any(NodeL_W(i) == cluster_A(:)) ||  any(NodeL_W(i) == cluster_B(:)))
                            obj.Nodes{i}.calculateConductance_Imbibition(obj, Pc);
                        end
                    end           
                % Water Saturation Calculation
                waterArea(i,1)=obj.Nodes{i}.area;
                waterArea(i,2)=obj.Nodes{i}.waterCrossSectionArea;
                if ~obj.Nodes{i}.isInlet && ~obj.Nodes{i}.isOutlet 
                    waterVolume = waterVolume + (obj.Nodes{i}.waterCrossSectionArea )...
                        / obj.Nodes{i}.area *obj.Nodes{i}.volume + obj.Nodes{i}.clayVolume;                
                    vol = vol + obj.Nodes{i}.volume + obj.Nodes{i}.clayVolume;
                end
            end
            
            for i = 1:obj.numberOfLinks 
                    if (any(LinkL(i) == cluster_A_nums(:)) ||  any(LinkL(i) == cluster_B_nums(:)))
                        obj.Links{i}.calculateConductance_Imbibition(obj, Pc);    
                    else
                        if isnan(obj.Links{i}.imbPressureTrapped) && obj.Links{i}.occupancy == 'B'
                            obj.Links{i}.imbPressureTrapped = max(Pc, obj.Links{i}.imbThresholdPressure_SnapOff);    
                        end 
                        if any(obj.Links{i}.imbPressureTrapped) 
                            obj.Links{i}.calculateConductance_Imbibition(obj, obj.Links{i}.imbPressureTrapped); 
                        else%if (any(LinkL_W(i) == cluster_A(:)) ||  any(LinkL_W(i) == cluster_B(:)))
                            obj.Links{i}.calculateConductance_Imbibition(obj, Pc);
                        end
                    end 
                % Water Saturation Calculation                
                waterArea(i,3)=obj.Links{i}.area;
                waterArea(i,4)=obj.Links{i}.waterCrossSectionArea;
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet 
                     waterVolume = waterVolume + (obj.Links{i}.waterCrossSectionArea )...
                         / obj.Links{i}.area * obj.Links{i}.volume + obj.Links{i}.clayVolume;
                     vol = vol + obj.Links{i}.volume + obj.Links{i}.clayVolume;
                end
            end             
            obj.waterSaturation = waterVolume / vol;     
        end
          
        % Clustering
        function [NumberOfClusters, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = Clustering_oil(obj)
            % Arguments of HoshenKopelman function
            NodeS = zeros(obj.numberOfNodes,1);
            LinksOfNode = zeros(obj.numberOfNodes,obj.maxCoordinationNumber);
            NodeNext = zeros(obj.numberOfNodes,obj.maxCoordinationNumber);
            
            for i = 1:obj.numberOfNodes
                NodeS(i,1) = obj.Nodes{i}.occupancy; % nodes with oil(1),nodes with water(0)
                if any(obj.Nodes{i}.oilLayerExist)
                    NodeS(i,1)='B';
                end
                if ~obj.Nodes{i}.isInlet && ~obj.Nodes{i}.isOutlet
                LinksOfNode(i,1:obj.Nodes{i}.connectionNumber) = obj.Nodes{i}.connectedLinks;
                NodeNext(i,1:obj.Nodes{i}.connectionNumber) = obj.Nodes{i}.connectedNodes;
                else
                    a = 1;
                    for j = 1:obj.Nodes{i}.connectionNumber
                        if ~obj.Links{obj.Nodes{i}.connectedLinks(j)}.isInlet && ~obj.Links{obj.Nodes{i}.connectedLinks(j)}.isOutlet                                                        
                            LinksOfNode(i,a) = obj.Nodes{i}.connectedLinks(j);
                            NodeNext(i,a) = obj.Nodes{i}.connectedNodes(j);
                            a = a+1;
                        end
                    end
                end
            end  
            LinkS = zeros(obj.numberOfLinks,1);
            for i =1:obj.numberOfLinks
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet
                    LinkS(i,1) = obj.Links{i}.occupancy; % throats with oil(1), throats with water(0)
                    if any(obj.Links{i}.oilLayerExist)
                        LinkS(i,1)='B';
                    end
                end
            end            
            OFlag = 'B'; %oil clusters has numbers 1:NumberOfClusters %water clusters are 0 
            
            % HoshenKopelman Algorithm for clustering
            [NumberOfClusters, NodeL, LinkL] = modifiedHKNonLattice(NodeS, LinkS,NodeNext, LinksOfNode, OFlag);
            a = 0;
            % Modify number of inlet & outlet Links of Clusters
            for i =1:obj.numberOfLinks
                if obj.Links{i}.isInlet 
                    if any(obj.Links{i}.oilLayerExist) || obj.Links{i}.occupancy == 'B'
                        if obj.Nodes{obj.Links{i}.pore2Index}.occupancy == 'B' || any(obj.Nodes{obj.Links{i}.pore2Index}.oilLayerExist)
                            LinkL(i,1) = NodeL(obj.Links{i}.pore2Index);
                        else
                            a = a + 1;
                            LinkL(i,1) = max(NodeL)+a;
                        end
                    end
                    
                elseif obj.Links{i}.isOutlet
                    if (obj.Links{i}.occupancy == 'B' || any(obj.Links{i}.oilLayerExist))
                        if obj.Nodes{obj.Links{i}.pore1Index}.occupancy == 'B' || any(obj.Nodes{obj.Links{i}.pore1Index}.oilLayerExist)
                            LinkL(i,1) = NodeL(obj.Links{i}.pore1Index);  
                        else
                            a = a + 1;
                            LinkL(i,1) = max(NodeL)+a;
                        end
                    end
                end
            end             
            
            inlet_cluster_indx = zeros(obj.numOfInletLinks,2);
            outlet_cluster_indx = zeros(obj.numOfOutletLinks,2);
            inlet = 1;
            outlet = 1;
            for i = 1:obj.numberOfLinks
                if obj.Links{i}.isInlet
                    inlet_cluster_indx(inlet,1) = obj.Links{i}.index;
                    inlet_cluster_indx(inlet,2) = LinkL(i,1);
                    inlet = inlet +1;
                elseif obj.Links{i}.isOutlet
                    outlet_cluster_indx(outlet,1) = obj.Links{i}.index;
                    outlet_cluster_indx(outlet,2) = LinkL(i,1);
                    outlet = outlet + 1;    
                end
            end                       
            
            a = 0;
            A = zeros(max(obj.numOfOutletLinks , obj.numOfInletLinks),1);
            for i = 1:length(outlet_cluster_indx)
                if outlet_cluster_indx(i,2) ~= 0
                    for j = 1:length(inlet_cluster_indx(:,2))
                        if outlet_cluster_indx(i,2) == inlet_cluster_indx(j,2)
                            if ~any(outlet_cluster_indx(i,2) == A(:,1))
                                a = a+1;
                                A(a,1) = outlet_cluster_indx(i,2);
                                break
                            end
                        end
                    end
                end
            end
            cluster_A_nums = nonzeros(A);    
            
            b = 0;
            B = zeros(obj.numberOfLinks,1);
            for i = 1:length(outlet_cluster_indx)
                if outlet_cluster_indx(i,2) ~= 0
                    if ~any(outlet_cluster_indx(i,2) == cluster_A_nums(:))
                        if ~any(outlet_cluster_indx(i,2) == B(:))
                            b = b +1;
                            B(b,1) = outlet_cluster_indx(i,2);
                        end
                    end
                end
            end
            cluster_B_nums = nonzeros(B); 
        end     
        function [NumberOfClusters, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = Clustering_water(obj)
            % Arguments of HoshenKopelman function
            NodeS = zeros(obj.numberOfNodes,1);
            LinksOfNode = zeros(obj.numberOfNodes,obj.maxCoordinationNumber);
            NodeNext = zeros(obj.numberOfNodes,obj.maxCoordinationNumber); 
            water = 0; 
            for i = 1:obj.numberOfNodes
                NodeS(i,1) = obj.Nodes{i}.occupancy; % nodes with oil(B),nodes with water(A)
                
                if obj.Nodes{i}.occupancy == 'A'
                    water = 1;
                end
                if obj.Nodes{i}.shapeFactor <= 1 / 16 % There is water in the corners
                    NodeS(i,1) = 'A';
                end
                if ~obj.Nodes{i}.isInlet && ~obj.Nodes{i}.isOutlet
                    LinksOfNode(i,1:obj.Nodes{i}.connectionNumber) = obj.Nodes{i}.connectedLinks;
                    NodeNext(i,1:obj.Nodes{i}.connectionNumber) = obj.Nodes{i}.connectedNodes;
                else
                    a = 1;
                    for j = 1:obj.Nodes{i}.connectionNumber
                        if ~obj.Links{obj.Nodes{i}.connectedLinks(j)}.isInlet && ~obj.Links{obj.Nodes{i}.connectedLinks(j)}.isOutlet                                                        
                            LinksOfNode(i,a) = obj.Nodes{i}.connectedLinks(j);
                            NodeNext(i,a) = obj.Nodes{i}.connectedNodes(j);
                            a = a+1;
                        end
                    end
                end
            end  
            LinkS = zeros(obj.numberOfLinks,1);
            for i =1:obj.numberOfLinks
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet 
                     LinkS(i,1) = obj.Links{i}.occupancy; % throats with oil(B), throats with water(A)
                    
                 if obj.Links{i}.occupancy == 'A'
                        water = 1;
                 end
                     if obj.Links{i}.shapeFactor <= 1 / 16 % There is water in the corners
                        LinkS(i,1) = 'A';  
                     end
                end
            end            
            OFlag = 'A'; %oil clusters has numbers B:NumberOfClusters %water clusters are A
            
            if  water == 1
                % HoshenKopelman Algorithm for clustering
                [NumberOfClusters, NodeL, LinkL] = modifiedHKNonLattice(NodeS, LinkS,NodeNext, LinksOfNode, OFlag);
                
                % Modify number of inlet & outlet Links of Clusters
                a = 0;
                for i =1:obj.numberOfLinks
                    if obj.Links{i}.isInlet
                        if  obj.Links{i}.occupancy == 'A' || any(obj.Links{i}.waterCornerExist)
                            if obj.Nodes{obj.Links{i}.pore2Index}.occupancy == 'A' || obj.Nodes{obj.Links{i}.pore2Index}.shapeFactor <= 1 / 16
                                LinkL(i,1) = NodeL(obj.Links{i}.pore2Index);
                            else
                                a = a + 1;
                                LinkL(i,1) = max(NodeL)+a;
                            end
                        end
                        
                    elseif obj.Links{i}.isOutlet
                        if obj.Links{i}.occupancy == 'A' || any(obj.Links{i}.waterCornerExist)
                            if obj.Nodes{obj.Links{i}.pore1Index}.occupancy == 'A' || obj.Nodes{obj.Links{i}.pore1Index}.shapeFactor <= 1 / 16
                                LinkL(i,1) = NodeL(obj.Links{i}.pore1Index);
                            else
                                a = a + 1;
                                LinkL(i,1) = max(NodeL)+a;
                            end
                        end
                    end
                end
                
                inlet_cluster_indx = zeros(obj.numOfInletLinks,2);
                outlet_cluster_indx = zeros(obj.numOfOutletLinks,2);
                inlet = 1;
                outlet = 1;
                for i = 1:obj.numberOfLinks
                    if obj.Links{i}.isInlet
                        inlet_cluster_indx(inlet,1) = obj.Links{i}.index;
                        inlet_cluster_indx(inlet,2) = LinkL(i,1);
                        inlet = inlet +1;
                    elseif obj.Links{i}.isOutlet
                        outlet_cluster_indx(outlet,1) = obj.Links{i}.index;
                        outlet_cluster_indx(outlet,2) = LinkL(i,1);
                        outlet = outlet + 1;
                    end
                end
                
                a = 0;
                A = zeros(max(obj.numOfOutletLinks , obj.numOfInletLinks),1);
                for i = 1:length(outlet_cluster_indx)
                    if outlet_cluster_indx(i,2) ~= 0
                        for j = 1:length(inlet_cluster_indx(:,2))
                            if outlet_cluster_indx(i,2) == inlet_cluster_indx(j,2)
                                if ~any(outlet_cluster_indx(i,2) == A(:,1))
                                    a = a+1;
                                    A(a,1) = outlet_cluster_indx(i,2);
                                    break
                                end
                            end
                        end
                    end
                end
                cluster_A_nums = nonzeros(A);
                
                b = 0;
                B = zeros(obj.numberOfLinks,1);
                for i = 1:length(outlet_cluster_indx)
                    if outlet_cluster_indx(i,2) ~= 0
                        if ~any(outlet_cluster_indx(i,2) == cluster_A_nums(:))
                            if ~any(outlet_cluster_indx(i,2) == B(:))
                                b = b +1;
                                B(b,1) = outlet_cluster_indx(i,2);
                            end
                        end
                    end
                end
                cluster_B_nums = nonzeros(B);
            else
                NumberOfClusters=-1; NodeL=-2*ones(obj.numberOfNodes,1); LinkL=-2*ones(obj.numberOfLinks,1); cluster_A_nums=-1; cluster_B_nums=-1;
            end
        end 
       
        % Pressure distribution calculation in pore_Two-Phases 
        function pressureDistribution_TwoPhases(obj, inletPressure, outletPressure) 
            
            Factor_W = zeros(obj.numberOfNodes, obj.numberOfNodes);
            B_W = zeros(obj.numberOfNodes, 1);
            Factor_O = zeros(obj.numberOfNodes, obj.numberOfNodes);
            B_O = zeros(obj.numberOfNodes, 1);
            
            % calculation of pressure distribution
            for ii = 1:obj.numberOfLinks  
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;
                
                % if the link is connected to inlet (index of node 1 is -1 which does not exist) 
                if obj.Links{ii}.isInlet
                    
                    obj.Links{ii}.nodeLinkSystemConductance_O = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.oilConductance) +...
                        ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.oilConductance)))^-1;     
                    
                    Factor_O(node2Index, node2Index) = Factor_O(node2Index, node2Index) + obj.Links{ii}.nodeLinkSystemConductance_O;
                    B_O(node2Index) = obj.Links{ii}.nodeLinkSystemConductance_O * inletPressure;
                    
                    obj.Links{ii}.nodeLinkSystemConductance_W = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.waterConductance) +...
                        ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.waterConductance)))^-1;   
                    
                    Factor_W(node2Index, node2Index) = Factor_W(node2Index, node2Index) + obj.Links{ii}.nodeLinkSystemConductance_W;
                    B_W(node2Index) = obj.Links{ii}.nodeLinkSystemConductance_W * inletPressure;
                    
                % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                elseif obj.Links{ii}.isOutlet
                    
                    obj.Links{ii}.nodeLinkSystemConductance_O = ( (obj.Links{ii}.linkLength /...
                        obj.Links{ii}.oilConductance) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.oilConductance)))^-1;
                    
                    Factor_O(node1Index, node1Index) = Factor_O(node1Index, node1Index) + obj.Links{ii}.nodeLinkSystemConductance_O;
                    B_O(node1Index) = obj.Links{ii}.nodeLinkSystemConductance_O * outletPressure;  
                    
                    obj.Links{ii}.nodeLinkSystemConductance_W = ( (obj.Links{ii}.linkLength /...
                        obj.Links{ii}.waterConductance) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.waterConductance)))^-1;
                    
                    Factor_W(node1Index, node1Index) = Factor_W(node1Index, node1Index) + obj.Links{ii}.nodeLinkSystemConductance_W;
                    B_W(node1Index) = obj.Links{ii}.nodeLinkSystemConductance_W * outletPressure;   
                    
                %if the link is neither inlet nor outlet    
                else
                    obj.Links{ii}.nodeLinkSystemConductance_W = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.waterConductance) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.waterConductance) +...
                        (obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.waterConductance)))^-1;  
                    
                    Factor_W(node1Index, node1Index) = Factor_W(node1Index, node1Index) + obj.Links{ii}.nodeLinkSystemConductance_W;
                    Factor_W(node2Index, node2Index) = Factor_W(node2Index, node2Index) + obj.Links{ii}.nodeLinkSystemConductance_W;
                    Factor_W(node1Index, node2Index) = Factor_W(node1Index, node2Index) - obj.Links{ii}.nodeLinkSystemConductance_W;
                    Factor_W(node2Index, node1Index) = Factor_W(node2Index, node1Index) - obj.Links{ii}.nodeLinkSystemConductance_W;  
                    
                    obj.Links{ii}.nodeLinkSystemConductance_O = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.oilConductance) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.oilConductance) +...
                        (obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.oilConductance)))^-1;                
                    Factor_O(node1Index, node1Index) = Factor_O(node1Index, node1Index) + obj.Links{ii}.nodeLinkSystemConductance_O;
                    Factor_O(node2Index, node2Index) = Factor_O(node2Index, node2Index) + obj.Links{ii}.nodeLinkSystemConductance_O;
                    Factor_O(node1Index, node2Index) = Factor_O(node1Index, node2Index) - obj.Links{ii}.nodeLinkSystemConductance_O;
                    Factor_O(node2Index, node1Index) = Factor_O(node2Index, node1Index) - obj.Links{ii}.nodeLinkSystemConductance_O;                   
                end     
            end
            
            % using Preconditioned conjugate gradients method to solve the
            % pressure distribution  
            nodesWaterPressure = gmres(Factor_W, B_W,[], 1e-10, obj.numberOfNodes);            
            nodesOilPressure = gmres(Factor_O, B_O,[], 1e-10, obj.numberOfNodes);
%             nodesWaterPressure = pcg(Factor_W, B_W, 1e-7, 1000);
%             nodesOilPressure = pcg(Factor_O, B_O, 1e-7, 1000);
            %assign the pressure values to each node
            for ii = 1:obj.numberOfNodes
                if nodesWaterPressure(ii)> inletPressure
                    obj.Nodes{ii}.waterPressure = inletPressure;
                elseif nodesWaterPressure(ii)< outletPressure
                    obj.Nodes{ii}.waterPressure = outletPressure;
                else
                obj.Nodes{ii}.waterPressure = nodesWaterPressure(ii);        
                end
                if nodesOilPressure(ii)> inletPressure
                    obj.Nodes{ii}.oilPressure = inletPressure;
                elseif nodesOilPressure(ii)< outletPressure
                    obj.Nodes{ii}.oilPressure = outletPressure;
                else
                obj.Nodes{ii}.oilPressure = nodesOilPressure(ii);        
                end
            end
            
            %assign pressure values to links, since the surface where
            %flowrate is calculated through might pass through the links
            for ii = 1:obj.numberOfLinks
                if obj.Links{ii}.isInlet
                    obj.Links{ii}.waterPressure =...
                        (1+obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure)/2;
                    obj.Links{ii}.oilPressure =...
                        (1+obj.Nodes{obj.Links{ii}.pore2Index}.oilPressure)/2;
                elseif obj.Links{ii}.isOutlet
                    obj.Links{ii}.waterPressure =...
                        obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure/2; 
                    obj.Links{ii}.oilPressure =...
                        obj.Nodes{obj.Links{ii}.pore1Index}.oilPressure/2;      
                else
                    obj.Links{ii}.waterPressure =...
                        (obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure + ...
                        obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure) / 2;
                    obj.Links{ii}.oilPressure =...
                        (obj.Nodes{obj.Links{ii}.pore1Index}.oilPressure + ...
                        obj.Nodes{obj.Links{ii}.pore2Index}.oilPressure) / 2;
                end
            end              
        end
        
        % Relative Permeability
        % Piri: eliminate boundary condition 
        function [krw, kro] = calculateRelativePermeability_Drainage_box(obj, outletPressure, LinkL, LinkL_W, cluster_A_nums, cluster_A)  
            
            waterFlowRate = 0;   
            oilFlowRate = 0;
            
            %search through all the links
             
            sumPA_a = 0;
            sumA_a = 0;
            sumPA_b = 0;
            sumA_b = 0;
            sumPA_ao = 0;
            sumA_ao = 0;
            sumPA_bo = 0;
            sumA_bo = 0;
            x_coor = zeros(obj.numberOfNodes,1);
            
            for ii = 1:obj.numberOfNodes
                x_coor(ii,1) = obj.Nodes{ii}.x_coordinate;
            end
            x_outlet = max(x_coor);
            x_inlet = min(x_coor);
            intervalx = (x_outlet - x_inlet)/20; % TODO: Should check > linkLength to avoid index -1 & 0 in nodes
            a = x_inlet + intervalx;
            b = x_outlet - intervalx;  
                for i = 1:obj.numberOfNodes 
                    if obj.Nodes{i}.x_coordinate == a
                        sumA_a = sumA_a + obj.Nodes{i}.waterCrossSectionArea ;
                        sumPA_a = sumPA_a + obj.Nodes{i}.waterCrossSectionArea * obj.Nodes{i}.waterPressure;
                        sumA_ao = sumA_ao + obj.Nodes{i}.oilCrossSectionArea ;
                        sumPA_ao = sumPA_ao + obj.Nodes{i}.oilCrossSectionArea * obj.Nodes{i}.oilPressure;
                    elseif obj.Nodes{i}.x_coordinate == b
                        sumA_b = sumA_b + obj.Nodes{i}.waterCrossSectionArea ;
                        sumPA_b = sumPA_b + obj.Nodes{i}.waterCrossSectionArea * obj.Nodes{i}.waterPressure;
                        sumA_bo = sumA_bo + obj.Nodes{i}.oilCrossSectionArea ;
                        sumPA_bo = sumPA_bo + obj.Nodes{i}.oilCrossSectionArea * obj.Nodes{i}.oilPressure;
                    end
                end
                for i = 1:obj.numberOfLinks
                    if ~obj.Links{i}.isOutlet && ~obj.Links{i}.isInlet
                        if obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate > a && obj.Nodes{obj.Links{i}.pore1Index}.x_coordinate < a || ...
                             obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate < a && obj.Nodes{obj.Links{i}.pore1Index}.x_coordinate > a
                            sumA_a = sumA_a + obj.Links{i}.waterCrossSectionArea ;
                            sumPA_a = sumPA_a + obj.Links{i}.waterCrossSectionArea * obj.Links{i}.waterPressure;
                            sumA_ao = sumA_ao + obj.Links{i}.oilCrossSectionArea ;
                            sumPA_ao = sumPA_ao + obj.Links{i}.oilCrossSectionArea * obj.Links{i}.oilPressure;
                        elseif obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate > b && obj.Nodes{obj.Links{i}.pore1Index}.x_coordinate < b || ...
                                obj.Nodes{obj.Links{i}.pore2Index}.x_coordinate < b && obj.Nodes{obj.Links{i}.pore1Index}.x_coordinate > b
                            sumA_b = sumA_b + obj.Links{i}.waterCrossSectionArea ;
                            sumPA_b = sumPA_b + obj.Links{i}.waterCrossSectionArea * obj.Links{i}.waterPressure;
                            sumA_bo = sumA_bo + obj.Links{i}.oilCrossSectionArea ;
                            sumPA_bo = sumPA_bo + obj.Links{i}.oilCrossSectionArea * obj.Links{i}.oilPressure;
                        end
                    elseif obj.Links{i}.isOutlet
                    
                        % calculate the flow rate of the fluid
%                         if any(LinkL_W(ii) == cluster_A(:))  
                        waterFlowRate = waterFlowRate + ...
                            obj.Links{i}.nodeLinkSystemConductance_W * ...
                            (obj.Nodes{obj.Links{i}.pore1Index}.waterPressure - outletPressure);  
%                         end
                        
%                         if any(LinkL(ii) == cluster_A_nums(:)) 
                        % calculate the flow rate of the fluid
                        oilFlowRate = oilFlowRate + ...
                            obj.Links{i}.nodeLinkSystemConductance_O * ...
                            (obj.Nodes{obj.Links{i}.pore1Index}.oilPressure - outletPressure);  
%                         end
                    end
                end  
            
                Pa = sumPA_a / sumA_a;
                Pb = sumPA_b / sumA_b;
                Pao = sumPA_ao / sumA_ao;
                Pbo = sumPA_bo / sumA_bo;   
            
            % calculate velocity through the network 
            watervelocity = waterFlowRate/(obj.yDimension * obj.zDimension);             
            oilvelocity =oilFlowRate/(obj.yDimension * obj.zDimension);       
            unitConvertor = 1/0.987*10^15; % unit conversion from m2 to miliDarcy
            kw = unitConvertor * watervelocity * (b-a) * obj.waterViscosity/ (Pa -Pb );  
            ko = unitConvertor * oilvelocity * (b-a) * obj.oilViscosity/ (Pao -Pbo );                        
            
            krw = kw/obj.absolutePermeability;
            kro = ko/obj.absolutePermeability;
            if krw > 1
                krw = 1;
            elseif krw <0
                krw = 0;
            end 
            if kro > 1
                kro = 1;
            elseif kro <0
                kro = 0;
            end
        end  
        
        function [krw, kro] = calculateRelativePermeability_Drainage(obj, outletPressure, LinkL, LinkL_W, cluster_A_nums, cluster_A)
              
            waterFlowRate = 0;   
            oilFlowRate = 0;
            
            %search through all the links
            for ii = 1:obj.numberOfLinks 
                  
                node1Index = obj.Links{ii}.pore1Index;
                if obj.Links{ii}.isOutlet   
                     
                        % calculate the flow rate of the fluid
%                         if any(LinkL_W(ii) == cluster_A(:))  
                        waterFlowRate = waterFlowRate + ...
                            obj.Links{ii}.nodeLinkSystemConductance_W * ...
                            (obj.Nodes{node1Index}.waterPressure - outletPressure);  
%                         end
                        
%                         if any(LinkL(ii) == cluster_A_nums(:)) 
                        % calculate the flow rate of the fluid
                        oilFlowRate = oilFlowRate + ...
                            obj.Links{ii}.nodeLinkSystemConductance_O * ...
                            (obj.Nodes{node1Index}.oilPressure - outletPressure);  
%                         end
                end 
            end                
            krw = waterFlowRate/obj.totalFlowRate;
            if krw > 1
                krw = 1;
            elseif krw <0
                krw = 0;
            end
            kro = oilFlowRate * obj.oilViscosity/obj.totalFlowRate / obj.waterViscosity;
            if kro > 1
                kro = 1;
            elseif kro <0
                kro = 0;
            end
        end
        function [krw, kro] = calculateRelativePermeability_Imbibition (obj, outletPressure, LinkL, LinkL_W, cluster_A_nums, cluster_A)
              
            waterFlowRate = 0;   
            oilFlowRate = 0; 
            %search through all the links
            for ii = 1:obj.numberOfLinks 
                
                node1Index = obj.Links{ii}.pore1Index;                
                if obj.Links{ii}.isOutlet  
                        % calculate the flow rate of the fluid
                        if any(LinkL_W(ii) == cluster_A(:))  
                        waterFlowRate = waterFlowRate + ...
                            abs(obj.Links{ii}.nodeLinkSystemConductance_W * ...
                            (outletPressure - obj.Nodes{node1Index}.waterPressure));
                        end
                        % calculate the flow rate of the fluid
                        if any(LinkL(ii) == cluster_A_nums(:))  
                        oilFlowRate = oilFlowRate + ...
                            abs(obj.Links{ii}.nodeLinkSystemConductance_O * ...
                            (outletPressure - obj.Nodes{node1Index}.oilPressure)); 
                        end
                end

            end
            krw = waterFlowRate/obj.totalFlowRate;
            if krw > 1
                krw = 1;
            elseif krw <0
                krw = 0;
            end 
            kro = oilFlowRate * obj.oilViscosity/obj.totalFlowRate / obj.waterViscosity;  
            if kro > 1
                kro = 1;
            elseif kro <0
                kro = 0;
            end
        end
        
        % Primary Drainage  
        function PrimaryDrainage(obj, inletPressure, outletPressure)              
             % determining the capillary pressure level interval
             Pc_threshold = zeros(2*obj.numberOfLinks,1);  
             Pc_threshold_n = zeros(obj.numberOfLinks,1); 
             for i = 1:obj.numberOfLinks   
                 obj.Links{i}.waterPressure = 0;
                 obj.Links{i}.calculateThresholdPressurePistonLike_drainage();
                 Pc_threshold(i,1) = obj.Links{i}.drainThresholdPressure_PistonLike;
             end 
             for i = 1:obj.numberOfNodes                
                  obj.Nodes{i}.waterPressure = 0;               
                  obj.Nodes{i}.calculateThresholdPressurePistonLike_drainage();
                  Pc_threshold(i+obj.numberOfLinks) = obj.Nodes{i}.drainThresholdPressure_PistonLike; 
             end  
             max_Pc = max(Pc_threshold); 
             min_Pc = min(nonzeros(Pc_threshold)); 
             Pc_interval = (max_Pc-min_Pc)/10;  
             
             Pc = min_Pc; 
             lastMaxP = Pc;
             t = 1; 
             invaded = 0;
             obj.Pc_drain_max = max_Pc;
             obj.DrainageData = zeros(10,5); 
             % initializing clusters
             LinkL = zeros(obj.numberOfLinks); 
             cluster_A_nums =[];
             [~, NodeL_W, LinkL_W,cluster_A,cluster_B] = Clustering_water(obj); 
             
             % Cycle of increasing Pressure
             while Pc <= obj.Pc_drain_max *1.0001 
              press = 1; 
             % Find new inlet-Links with threshold pressure < Pc             
             for i = 1:obj.numberOfLinks                  
                  node1Index = obj.Links{i}.pore1Index;
                  node2Index = obj.Links{i}.pore2Index;
                  if obj.Links{i}.isInlet && obj.Links{i}.occupancy == 'A'
                         if (any(LinkL_W(i) == cluster_A(:))) && Pc_threshold(i) <= Pc 
                             obj.Links{i}.occupancy = 'B';                              
                             invaded = invaded + 1;   
                             if  obj.Nodes{node2Index}.occupancy == 'A' &&obj.Nodes{node2Index}.drainThresholdPressure_PistonLike <=Pc  &&...
                                     (any(NodeL_W(node2Index) == cluster_A(:)) ||  any(NodeL_W(node2Index) == cluster_B(:)))
                                 
                                 obj.Nodes{node2Index}.occupancy = 'B';
                                 invaded = invaded + 1; 
                                 for j=1:obj.Nodes{node2Index}.connectionNumber
                                     if obj.Nodes{node2Index}.connectedLinks(j)~=i
                                         Pc_threshold_n(obj.Nodes{node2Index}.connectedLinks(j),1)= Pc_threshold(obj.Nodes{node2Index}.connectedLinks(j));
                                     end
                                 end
                             end
                         end
                  elseif obj.Links{i}.isOutlet && obj.Links{i}.occupancy == 'A'                     
                           if obj.Nodes{node1Index}.occupancy == 'B'  
                               Pc_threshold_n(i,1)= Pc_threshold(i); 
                           end
                  elseif ~obj.Links{i}.isOutlet && ~obj.Links{i}.isInlet && obj.Links{i}.occupancy == 'A'                    
                      if obj.Nodes{node1Index}.occupancy == 'B' || obj.Nodes{node2Index}.occupancy == 'B'
                          if (any(LinkL_W(i) == cluster_A(:)) ||  any(LinkL_W(i) == cluster_B(:)))
                              Pc_threshold_n(i,1)= Pc_threshold(i);  
                          end
                      end  
                  end 
             end
             deltaS = 0;
             deltaV = 0;  
             if (any(nonzeros(Pc_threshold_n)) && min(nonzeros(Pc_threshold_n))<= Pc) %|| Pc == max_Pc
                 pressure = 1;
             else
                 pressure = 0;
             end 
             lastMaxP = Pc;
             % Add Links which have Pc_threshold < Pc in each steps and also have oil-saturated neighbour Node 
             while  pressure == 1 && deltaS <= 0.1 
                    
                 %check & sort Links based on Pc_Threshold
                 [~, ix] = sort(Pc_threshold_n(1:end), 1);
                 i = ix(obj.numberOfLinks - length(nonzeros(Pc_threshold_n))+1);
                 
                 if lastMaxP < Pc_threshold_n(i)
                     lastMaxP = Pc_threshold_n(i);
                 end
                 Pc_threshold_n(i) = 0;   
                 node1Index = obj.Links{i}.pore1Index;
                 node2Index = obj.Links{i}.pore2Index;
                 
                 if obj.Links{i}.isOutlet && obj.Links{i}.occupancy == 'A'                     
                           if obj.Nodes{node1Index}.occupancy == 'B'  
                               
                               obj.Links{i}.occupancy = 'B'; 
                               invaded = invaded + 1;
                               deltaV = deltaV + obj.Links{i}.volume ;  
                           end
                  elseif ~obj.Links{i}.isOutlet && ~obj.Links{i}.isInlet && obj.Links{i}.occupancy == 'A'                        
                      if obj.Nodes{node1Index}.occupancy == 'B' || obj.Nodes{node2Index}.occupancy == 'B'  &&...
                              (any(LinkL_W(i) == cluster_A(:)) ||  any(LinkL_W(i) == cluster_B(:)))
                          
                          obj.Links{i}.occupancy = 'B';
                          invaded = invaded + 1;
                          deltaV = deltaV + obj.Links{i}.volume ;  
                          
                          if  obj.Nodes{node2Index}.occupancy == 'A' && obj.Nodes{node2Index}.drainThresholdPressure_PistonLike <=Pc  && ...                                  
                                   (any(NodeL_W(node2Index) == cluster_A(:)) ||  any(NodeL_W(node2Index) == cluster_B(:)))
                               
                              obj.Nodes{node2Index}.occupancy = 'B'; 
                              invaded = invaded + 1;
                              deltaV = deltaV + obj.Nodes{node2Index}.volume ;  
                              for j=1:obj.Nodes{node2Index}.connectionNumber
                                 if obj.Nodes{node2Index}.connectedLinks(j)~=i
                                     Pc_threshold_n(obj.Nodes{node2Index}.connectedLinks(j),1)= Pc_threshold(obj.Nodes{node2Index}.connectedLinks(j));
                                 end
                             end
                          end   
                          
                          if obj.Nodes{node1Index}.occupancy == 'A' && obj.Nodes{node1Index}.drainThresholdPressure_PistonLike <=Pc && ...
                                   (any(NodeL_W(node1Index) == cluster_A(:)) ||  any(NodeL_W(node1Index) == cluster_B(:)))
                               
                              obj.Nodes{node1Index}.occupancy = 'B';
                              deltaV = deltaV + obj.Nodes{node1Index}.volume ; 
                              invaded = invaded + 1; 
                              for j=1:obj.Nodes{node1Index}.connectionNumber
                                 if obj.Nodes{node1Index}.connectedLinks(j)~=i
                                     Pc_threshold_n(obj.Nodes{node1Index}.connectedLinks(j),1)= Pc_threshold(obj.Nodes{node1Index}.connectedLinks(j));
                                 end
                             end
                          end
                      end                      
                 end                 
                 deltaS = deltaV /obj.poreVolume ; 
                   if deltaS > 0.1
                     press = 0;
                   end
                 if any (nonzeros(Pc_threshold_n)) && min(nonzeros(Pc_threshold_n))<= Pc
                     pressure = 1;
                 else
                     pressure = 0;
                 end
                 if invaded == 1 || invaded == 2
                     break
                 end
             end
               if Pc == max_Pc
                   lastMaxP = max_Pc;
               end
               
             % Updating element saturations and conductances
             calculateConductance_and_Saturation_Drainage(obj, lastMaxP); 
             pressureDistribution_TwoPhases(obj, inletPressure, outletPressure); 
             
             % Relative Permeability Calculation              
             [Krw , Kro] = calculateRelativePermeability_Drainage (obj, outletPressure, LinkL, LinkL_W, cluster_A_nums, cluster_A);             
             obj.DrainageData(t,:) = [obj.waterSaturation, lastMaxP, Krw, Kro, invaded]; 
%              fileName = obj.vtkWriter('PD',t);
             fileName = obj.vtkWriter_glyph('PD',t);
             if invaded ~= 0                 
                 [~, ~, LinkL,cluster_A_nums,~] = Clustering_oil(obj);  
             end
             [~, NodeL_W, LinkL_W,cluster_A,cluster_B] = Clustering_water(obj);
             
              % Pc Step Calculation 
              if press ~= 0 
                 Pc = Pc + Pc_interval;   
              end
             t = t + 1;     
             end               
             
        end   
         
        % Secondary Imbibition 
        function ScoendaryImbibition(obj, inletPressure, outletPressure)    
            
            %counter for invaded elements
            numOfLinks_SnapOff = 0;
            numOfLinks_PistoneLike = 0;
            numOfLinks_LayerCollapse = 0;
            numOfNodes_SnapOff = 0;
            numOfNodes_PoreBodyFilling = 0;
            numOfNodes_LayerCollapse = 0;
            obj.thresholdPressure = zeros(obj.numberOfLinks, 14);               
            Pc_imb = obj.Pc_drain_max; 
            Pc_min = Pc_imb;
            
            % Calculating throat Snap-Off & Pistone-Like displacement & layer collapse            
            for i = 1:obj.numberOfLinks  
                
%                 newContactAngle = obj.Links{i}.advancingContactAngle - obj.Links{1}.wettabilityAlteration;
                if obj.Links{i}.occupancy == 'B' % if the throat is oil filled   
%                     obj.Links{i}.advancingContactAngle = obj.Links{i}.newContactAngle;
                    obj.Links{i}.calculateThresholdPressurePistonLike_Imbibition (obj.Pc_drain_max);
                    obj.Links{i}.calculateThresholdPressureSnapOff (obj.Pc_drain_max); 
%                     obj.Links{i}.calculateThresholdPressureLayerCollapse (obj.Pc_drain_max);
                    if obj.Links{i}.isInlet
                        obj.thresholdPressure(i,1) = -1; 
                    elseif obj.Links{i}.isOutlet
                        obj.thresholdPressure(i,1) = 1; 
                    end
                    obj.thresholdPressure(i,2) = obj.Links{i}.imbThresholdPressure_PistonLike;
                    obj.thresholdPressure(i,3) = obj.Links{i}.imbThresholdPressure_SnapOff; 
%                     obj.thresholdPressure(i,4:7) = obj.Links{i}.imbThresholdPressure_LayerCollapse;
%                     obj.Links{i}.calculateThresholdPressurePistonLike_Imbibition_R (obj.Pc_drain_max); 
%                     obj.thresholdPressure(i,3) = obj.Links{i}.imbThresholdPressure_PistonLike;
%                     obj.Links{i}.calculateThresholdPressurePistonLike_Imbibition_P (obj.Pc_drain_max); 
%                     obj.thresholdPressure(i,4) = obj.Links{i}.imbThresholdPressure_PistonLike;
%                     obj.Links{i}.calculateThresholdPressurePistonLike_Imbibition_2003 (obj.Pc_drain_max); 
%                     obj.thresholdPressure(i,5) = obj.Links{i}.imbThresholdPressure_PistonLike;
                end
            end
            
            % Calculating Pore Snap-Off & Pore-Body Filling displacement & layer collapse            
            for i = 1:obj.numberOfNodes                
%                 obj.Nodes{i}.advancingContactAngle = obj.Nodes{i}.newContactAngle;
                if obj.Nodes{i}.occupancy == 'B' % if the throat is oil filled
                    obj.Nodes{i}.calculateThresholdPressurePoreBodyFilling (obj);
                    obj.Nodes{i}.calculateThresholdPressurePistonLike_Imbibition (obj.Pc_drain_max);
                    obj.Nodes{i}.calculateThresholdPressureSnapOff (obj.Pc_drain_max);
                    %                     obj.Nodes{i}.calculateThresholdPressureLayerCollapse (obj.Pc_drain_max);
                    if obj.Nodes{i}.isInlet
                        obj.thresholdPressure(i,8) = -1;
                    elseif obj.Nodes{i}.isOutlet
                        obj.thresholdPressure(i,8) = 1;
                    end
                    obj.thresholdPressure(i,9) = obj.Nodes{i}.imbThresholdPressure_PistonLike;
                    obj.thresholdPressure(i,10) = obj.Nodes{i}.imbThresholdPressure_SnapOff;
                    %                     obj.thresholdPressure(i,11:14) = obj.Nodes{i}.imbThresholdPressure_LayerCollapse;
                    %                     obj.Nodes{i}.calculateThresholdPressurePistonLike_Imbibition_R (obj.Pc_drain_max);
                    %                     obj.thresholdPressure(i,10) = obj.Nodes{i}.imbThresholdPressure_PistonLike;
                    %                     obj.Nodes{i}.calculateThresholdPressurePistonLike_Imbibition_P (obj.Pc_drain_max);
                    %                     obj.thresholdPressure(i,11) = obj.Nodes{i}.imbThresholdPressure_PistonLike;
                    %                     obj.Nodes{i}.calculateThresholdPressurePistonLike_Imbibition_2003 (obj.Pc_drain_max);
                    %                     obj.thresholdPressure(i,12) = obj.Nodes{i}.imbThresholdPressure_PistonLike;
                end
            end
            
            Pc_interval = Pc_imb /10;
            t = 0;
            obj.ImbibitionData = zeros(100,12);
            invaded_Element = zeros(2*(obj.numberOfLinks+obj.numberOfNodes), 11);
            e = 0;
            obj.sequence = zeros(2*obj.numberOfLinks, 11);
            percList = -1000000*ones(obj.numberOfNodes+obj.numberOfLinks,1);
            poreVolumeInjected = 0;
            PVInjected = 0;
            [~, NodeL, LinkL, cluster_A_nums, cluster_B_nums] = Clustering_oil(obj);
            inv = false;
            
            % Invasion & Percolation List
            for i = 1:obj.numberOfLinks
                
                node1Index = obj.Links{i}.pore1Index;
                node2Index = obj.Links{i}.pore2Index;
                
                if (any(LinkL(i) == cluster_A_nums(:)) || any(LinkL(i) == cluster_B_nums(:)))
                    
                    if obj.Links{i}.isInlet
                        
                        if any(obj.Links{i}.imbThresholdPressure_PistonLike)
                            
                            percList(i) = obj.Links{i}.imbThresholdPressure_PistonLike ;
                        end
                    elseif obj.Links{i}.isOutlet
                        
                        if  obj.Nodes{node1Index}.occupancy == 'A' && any(obj.Links{i}.imbThresholdPressure_PistonLike)
                            
                            percList(i) = obj.Links{i}.imbThresholdPressure_PistonLike ;
                        elseif  obj.Nodes{node1Index}.occupancy == 'B' && ...
                                any(obj.Links{i}.imbThresholdPressure_SnapOff) % if the throat is non circular
                            
                            percList(i) = obj.Links{i}.imbThresholdPressure_SnapOff;
                        end
                    else
                        if  (obj.Nodes{node1Index}.occupancy == 'A' && obj.Nodes{node2Index}.occupancy == 'B') || ...
                                (obj.Nodes{node1Index}.occupancy == 'B' && obj.Nodes{node2Index}.occupancy == 'A') && ...
                                any(obj.Links{i}.imbThresholdPressure_PistonLike)
                            
                            
                            percList(i) = obj.Links{i}.imbThresholdPressure_PistonLike ;
                            
                        elseif obj.Nodes{node1Index}.occupancy == 'B' &&...
                                obj.Nodes{node2Index}.occupancy == 'B' &&...
                                any(obj.Links{i}.imbThresholdPressure_SnapOff)% if the throat is non circular
                            
                            percList(i) = obj.Links{i}.imbThresholdPressure_SnapOff;
                        end
                    end
                end
            end
            a = obj.numberOfLinks;
            for i = 1:obj.numberOfNodes
                if (any(NodeL(i) == cluster_A_nums(:)) || any(NodeL(i) == cluster_B_nums(:)))
                    
                    filledThroats = 0;
                    for j = 1:obj.Nodes{i}.connectionNumber
                        if (obj.Links{obj.Nodes{i}.connectedLinks(j)}.occupancy == 'A')
                            filledThroats = filledThroats + 1;
                        end
                    end
                    
                    if filledThroats ~= 0 &&  any(obj.Nodes{i}.imbThresholdPressure_PoreBodyFilling)
                        
                        percList(a+i) = obj.Nodes{i}.imbThresholdPressure_PoreBodyFilling;
                        
                    elseif filledThroats == 0 && any(obj.Nodes{i}.imbThresholdPressure_SnapOff)% if the node is non circular
                        
                        percList(a+i) = obj.Nodes{i}.imbThresholdPressure_SnapOff; % snap off threshold pressure
                    end
                    
                end
            end
            
            while (~isempty(cluster_A_nums) || ~isempty(cluster_B_nums)) && Pc_imb >-99999
                
                press = 1;
                deltaS = 0;
                
                % Percolation Section
                if (max(percList)) >= Pc_imb
                    pressure = 1;
                else
                    pressure = 0;
                end
                
                while pressure == 1 && (max(percList)) >= Pc_imb  && deltaS <= 0.2
                    
                    inv = true;
                    % Descending sorting of threshold pressures
                    [PcTh, ix] = max(percList(1:end));
                    if PcTh ~= -1000000
                        indexElement = ix(1);
                        if Pc_min > percList(indexElement)
                            Pc_min = percList(indexElement);
                        end
                        
                        [~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = Clustering_oil(obj);
                        [~, NodeL_W, LinkL_W,cluster_A,cluster_B] = Clustering_water(obj);
                        
                        % if the first element is a throat
                        if indexElement <= obj.numberOfLinks
                            
                            linkIndex = indexElement;
                            node1Index = obj.Links{linkIndex}.pore1Index;
                            node2Index = obj.Links{linkIndex}.pore2Index;
                            if any(LinkL(linkIndex) == cluster_A_nums(:)) || any(LinkL(linkIndex) == cluster_B_nums(:))
                                
                                if obj.Links{linkIndex}.isInlet
                                    
                                    if obj.Links{linkIndex}.imbThresholdPressure_PistonLike >= Pc_imb
                                        
                                        obj.Links{linkIndex}.occupancy = 'A';
                                        obj.Links{linkIndex}.oilLayerExistance()
                                        obj.Links{linkIndex}.calculateThresholdPressureLayerCollapse(obj.Pc_drain_max);
                                        poreVolumeInjected = poreVolumeInjected + obj.Links{linkIndex}.volume;
                                        numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;
                                        obj.Links{linkIndex}.isInvaded = true;
                                        e = e+1;
                                        invaded_Element(e,1:4) = [linkIndex, percList(linkIndex), ...
                                            obj.Links{linkIndex}.imbThresholdPressure_PistonLike, obj.Links{linkIndex}.imbThresholdPressure_SnapOff];
                                        
                                        if  obj.Nodes{node2Index}.occupancy == 'B'
                                            % Updating pore body filling of the pore
                                            obj.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                            if any(obj.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling)   
                                                percList(obj.numberOfLinks+node2Index) = obj.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling;
                                                 
                                            end
                                        end
                                    end
                                    
                                elseif obj.Links{linkIndex}.isOutlet
                                    if obj.Links{linkIndex}.imbThresholdPressure_PistonLike >= Pc_imb && obj.Nodes{node1Index}.occupancy == 'A'
                                        
                                        obj.Links{linkIndex}.occupancy = 'A';
                                        obj.Links{linkIndex}.oilLayerExistance()
                                        obj.Links{linkIndex}.calculateThresholdPressureLayerCollapse(obj.Pc_drain_max);
                                        
                                        poreVolumeInjected = poreVolumeInjected + obj.Links{linkIndex}.volume;
                                        numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;
                                        obj.Links{linkIndex}.isInvaded = true;
                                        e = e+1;
                                        invaded_Element(e,1:4) = [linkIndex, percList(linkIndex), ...
                                            obj.Links{linkIndex}.imbThresholdPressure_PistonLike, obj.Links{linkIndex}.imbThresholdPressure_SnapOff];
                                        
                                    elseif obj.Links{linkIndex}.imbThresholdPressure_SnapOff >= Pc_imb && obj.Nodes{node1Index}.occupancy == 'B'
                                        
                                        obj.Links{linkIndex}.occupancy = 'A';
                                        obj.Links{linkIndex}.oilLayerExistance()
                                        obj.Links{linkIndex}.calculateThresholdPressureLayerCollapse(obj.Pc_drain_max);
                                        numOfLinks_SnapOff = numOfLinks_SnapOff + 1;
                                        poreVolumeInjected = poreVolumeInjected + obj.Links{linkIndex}.volume;
                                        obj.Links{linkIndex}.isInvaded = true;
                                        e = e+1;
                                        invaded_Element(e,1:4) = [linkIndex, percList(linkIndex), ...
                                            obj.Links{linkIndex}.imbThresholdPressure_PistonLike, obj.Links{linkIndex}.imbThresholdPressure_SnapOff];
                                        percList(obj.numberOfLinks+node1Index) = -1000000;
                                        
                                        % Updating pore body filling of the pore
                                        obj.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                        if any(obj.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling)  
                                            
                                            percList(obj.numberOfLinks+node1Index) = obj.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling;
                                        end
                                    end
                                    
                                elseif ~obj.Links{linkIndex}.isOutlet && ~obj.Links{linkIndex}.isInlet
                                    
                                    if obj.Nodes{node1Index}.occupancy == 'A' || obj.Nodes{node2Index}.occupancy == 'A'
                                        if obj.Links{linkIndex}.imbThresholdPressure_PistonLike >= Pc_imb
                                            
                                            obj.Links{linkIndex}.occupancy = 'A';
                                            obj.Links{linkIndex}.oilLayerExistance()
                                            obj.Links{linkIndex}.calculateThresholdPressureLayerCollapse(obj.Pc_drain_max);
                                            poreVolumeInjected = poreVolumeInjected + obj.Links{linkIndex}.volume;
                                            numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;
                                            obj.Links{linkIndex}.isInvaded = true;
                                            e = e+1;
                                            invaded_Element(e,1:4) = [linkIndex, percList(linkIndex), ...
                                                obj.Links{linkIndex}.imbThresholdPressure_PistonLike, obj.Links{linkIndex}.imbThresholdPressure_SnapOff];
                                            
                                            if obj.Nodes{node1Index}.occupancy == 'B'
                                                
                                                % Updating pore body filling of the pore
                                                obj.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                                if any(obj.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling)   
                                                    percList(obj.numberOfLinks+node1Index) = obj.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling;
                                                end
                                            end
                                            if  obj.Nodes{node2Index}.occupancy == 'B'
                                                
                                                % Updating pore body filling of the pore
                                                obj.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                                if any(obj.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling)  
                                                    percList(obj.numberOfLinks+node2Index) = obj.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling;
                                                end
                                            end
                                        end
                                    elseif obj.Links{linkIndex}.imbThresholdPressure_SnapOff >= Pc_imb
                                        
                                        obj.Links{linkIndex}.occupancy = 'A';
                                        obj.Links{linkIndex}.oilLayerExistance()
                                        obj.Links{linkIndex}.calculateThresholdPressureLayerCollapse(obj.Pc_drain_max);
                                        numOfLinks_SnapOff = numOfLinks_SnapOff + 1;
                                        
                                        % Updating pore body filling of the pore
                                        obj.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                        if any(obj.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling)  
                                            percList(obj.numberOfLinks+node1Index) = obj.Nodes{node1Index}.imbThresholdPressure_PoreBodyFilling;
                                        end
                                        
                                        % Updating pore body filling of the pore
                                        obj.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                        if any(obj.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling) 
                                            percList(obj.numberOfLinks+node2Index) = obj.Nodes{node2Index}.imbThresholdPressure_PoreBodyFilling;
                                        end 
                                    end
                                end
                            end
                            
                            % if the first element is a pore
                        else
                            nodeIndex = indexElement-obj.numberOfLinks;
                            if any(NodeL(nodeIndex) == cluster_A_nums(:)) || any(NodeL(nodeIndex) == cluster_B_nums(:))
                                
                                filledThroats = 0;
                                for j = 1:obj.Nodes{nodeIndex}.connectionNumber
                                    
                                    if obj.Links{obj.Nodes{nodeIndex}.connectedLinks(j)}.occupancy == 'A'
                                        
                                        filledThroats = filledThroats + 1;
                                    end
                                end
                                
                                if filledThroats ~= 0 &&  any(obj.Nodes{nodeIndex}.imbThresholdPressure_PoreBodyFilling) && ...
                                        obj.Nodes{nodeIndex}.imbThresholdPressure_PoreBodyFilling >= Pc_imb
                                    
                                    obj.Nodes{nodeIndex}.occupancy = 'A'; % make the pore water type
                                    obj.Nodes{nodeIndex}.oilLayerExistance()
                                    obj.Nodes{nodeIndex}.calculateThresholdPressureLayerCollapse(obj.Pc_drain_max);
                                    poreVolumeInjected = poreVolumeInjected + obj.Nodes{nodeIndex}.volume;
                                    numOfNodes_PoreBodyFilling = numOfNodes_PoreBodyFilling + 1;
                                    obj.Nodes{nodeIndex}.isInvaded = true;
                                    e = e+1;
                                    invaded_Element(e,5:8) = [nodeIndex, percList(obj.numberOfLinks+nodeIndex), ...
                                        obj.Nodes{nodeIndex}.imbThresholdPressure_PoreBodyFilling, obj.Nodes{nodeIndex}.imbThresholdPressure_SnapOff];
                                    
                                    for j=1:obj.Nodes{nodeIndex}.connectionNumber
                                        percList(obj.Nodes{nodeIndex}.connectedLinks(j)) = -1000000;
                                        if obj.Links{obj.Nodes{nodeIndex}.connectedLinks(j)}.occupancy == 'B'
                                            percList(obj.Nodes{nodeIndex}.connectedLinks(j))=...
                                                obj.Links{obj.Nodes{nodeIndex}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                                        end
                                    end
                                    
                                elseif filledThroats == 0 && any(obj.Nodes{nodeIndex}.imbThresholdPressure_SnapOff) && ...% if the node is non circular
                                        obj.Nodes{nodeIndex}.imbThresholdPressure_SnapOff >= Pc_imb
                                    
                                    obj.Nodes{nodeIndex}.occupancy = 'A'; % make the pore water type
                                    obj.Nodes{nodeIndex}.oilLayerExistance()
                                    obj.Nodes{nodeIndex}.calculateThresholdPressureLayerCollapse(obj.Pc_drain_max);
                                    poreVolumeInjected = poreVolumeInjected + obj.Nodes{nodeIndex}.volume;
                                    numOfNodes_SnapOff = numOfNodes_SnapOff + 1;
                                    obj.Nodes{nodeIndex}.isInvaded = true;
                                    e = e+1;
                                    invaded_Element(e,5:8) = [nodeIndex, percList(obj.numberOfLinks+nodeIndex), ...
                                        obj.Nodes{nodeIndex}.imbThresholdPressure_PoreBodyFilling, obj.Nodes{nodeIndex}.imbThresholdPressure_SnapOff];
                                    
                                    for j=1:obj.Nodes{nodeIndex}.connectionNumber
                                        percList(obj.Nodes{nodeIndex}.connectedLinks(j)) = -1000000;
                                        if obj.Links{obj.Nodes{nodeIndex}.connectedLinks(j)}.occupancy == 'B'
                                            percList(obj.Nodes{nodeIndex}.connectedLinks(j))=...
                                                obj.Links{obj.Nodes{nodeIndex}.connectedLinks(j)}.imbThresholdPressure_PistonLike;
                                        end
                                    end
                                end
                            end
                        end
                        
                        percList(indexElement) = -1000000;
                        deltaS = poreVolumeInjected /obj.poreVolume ;
                        if deltaS > 0.2
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
                    for ii = 1:obj.numberOfNodes
                        
                        if any(obj.Nodes{ii}.oilLayerExist) && any(obj.Nodes{ii}.imbThresholdPressure_LayerCollapse(1,:))...
                                && (any(NodeL(ii) == cluster_A_nums(:)) || any(NodeL(ii) == cluster_B_nums(:)) )
                            
                            % Updating Pc of layer collapse
                            % Cheking layer collapse
                            for jj = 1:4
                                if ~isnan(obj.Nodes{ii}.imbThresholdPressure_LayerCollapse(1,j)) && ...
                                        obj.Nodes{ii}.imbThresholdPressure_LayerCollapse(1,j) > Pc_imb
                                    
                                    obj.Nodes{ii}.oilLayerExist(1,j) = nan;
                                    
                                    numOfNodes_LayerCollapse = numOfNodes_LayerCollapse + 1;
                                end
                            end
                        end
                    end
                    for ii = 1:obj.numberOfLinks
                        
                        if any(obj.Links{ii}.oilLayerExist) && any(obj.Links{ii}.imbThresholdPressure_LayerCollapse(1,:))...
                                && (any(LinkL(ii) == cluster_A_nums(:))|| any(LinkL(ii) == cluster_B_nums(:)) )
                            
                            % Updating Pc of layer collapse
                            % Cheking layer collapse
                            for jj = 1:4
                                
                                if ~isnan(obj.Links{ii}.imbThresholdPressure_LayerCollapse(1,j)) && ...
                                        obj.Links{ii}.imbThresholdPressure_LayerCollapse(1,j) > Pc_imb
                                    
                                    obj.Links{ii}.oilLayerExist(1,j) = nan;
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
                    calculateConductance_and_Saturation_Imb(obj, Pc_imb,NodeL, NodeL_W, LinkL, LinkL_W, cluster_A_nums, cluster_A, cluster_B_nums, cluster_B);
                    pressureDistribution_TwoPhases(obj, inletPressure, outletPressure);                      

                    PVInjected = poreVolumeInjected + PVInjected ;
%                    obj.calculateReactiveTransport_TwoPhaseDesorption(inletPressure, outletPressure, 1, poreVolumeInjected, newContactAngle)
                    poreVolumeInjected = 0;
                    [Krw_imb, Kro_imb] =...
                        calculateRelativePermeability_Imbibition(obj, outletPressure, LinkL, LinkL_W, cluster_A_nums, cluster_A);
                    obj.ImbibitionData(t,:) = ...
                        [obj.waterSaturation,Pc_imb,Krw_imb, Kro_imb,invaded, ...
                        numOfLinks_SnapOff,numOfLinks_PistoneLike, ...
                        numOfLinks_LayerCollapse,numOfNodes_SnapOff, ...
                        numOfNodes_PoreBodyFilling,numOfNodes_LayerCollapse,Pc_min];
                    
                    fileName = obj.vtkWriter_glyph('SI',t);
                    [~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = Clustering_oil(obj);
                    [~, NodeL_W, LinkL_W,cluster_A,cluster_B] = Clustering_water(obj);
                end
                inv = false;
                if press ~= 0
                    Pc_imb = Pc_imb - Pc_interval;
                end
                
            end
            
            obj.ImbibitionData = obj.ImbibitionData(1:t,:);
            obj.sequence(1:obj.numberOfLinks,1:9) = invaded_Element(1:obj.numberOfLinks,1:9);
        end 
        
        % Two-phase Reactive transportBased on Raoof paper 2013
        function calculateReactiveTransport_TwoPhaseDesorption(obj, inletPressure, outletPressure, soluteConcentration, poreVolumeInjected)
               
            % mass transfer coefficient: alpha & distribution coefficient: K_d
            alpha_link = 0.001*ones(obj.numberOfLinks,1);
            Kd_link = 0.0005*ones(obj.numberOfLinks,1);
            alpha_node = 0.001*ones(obj.numberOfNodes,1);
            Kd_node = 0.0005*ones(obj.numberOfNodes,1);
            
            obj.capillaryNumber = 1;
            
%             while obj.capillaryNumber > 10^(-7) % for Capillary dominant flow
%             inletPressure = inletPressure/2;
            
            % calculate pressure distribution
            pressureDistribution_TwoPhases(obj, inletPressure, outletPressure, 1000); 
            
            residenceTime_link = zeros(obj.numberOfLinks,1);
            residenceTime_node = zeros(obj.numberOfNodes,1);
            waterVolume_link = zeros(obj.numberOfLinks,1);  
            waterVolume_node = zeros(obj.numberOfNodes,1);
            flowRate_link = zeros(obj.numberOfLinks,1);            
            flowRate_node = zeros(obj.numberOfNodes,1);
            
            obj.totalFlowRate = 0;
            
            % calculate flowrate of links residence time                          
            for ii = 1:obj.numberOfLinks                
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;
                
                if ~obj.Links{ii}.isInlet && ~obj.Links{ii}.isOutlet      
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.nodeLinkSystemConductance_W * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);  
                    
                elseif obj.Links{ii}.isInlet 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.nodeLinkSystemConductance_W * ...
                        abs(inletPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);   
                else 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.nodeLinkSystemConductance_W * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        outletPressure);                    
                    obj.totalFlowRate = obj.totalFlowRate + flowRate_link(ii);
                end 
                waterVolume_link(ii) = obj.Links{ii}.waterCrossSectionArea * ...
                    obj.Links{ii}.linkLength;
                residenceTime_link(ii) = waterVolume_link(ii) / flowRate_link(ii);
            end 
            
            obj.velocity = obj.totalFlowRate * obj.xDimension / obj.poreVolume;
            obj.capillaryNumber = obj.waterViscosity * obj.velocity/ obj.sig_ow;            
%             end
            for i = 1:obj.numberOfNodes  
                for j = 1:obj.Nodes{i}.connectionNumber 
                    
                    connectedLinkIndex = obj.Nodes{i}.connectedLinks(j);
                    connectedNodeIndex = obj.Nodes{i}.connectedNodes(j); 
                    
                    if connectedNodeIndex ~= 0 && connectedNodeIndex ~= -1 
                        % determine link flowing into this node
                        if obj.Nodes{connectedNodeIndex}.waterPressure > obj.Nodes{i}.waterPressure  
                            flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex);   
                        end
                    elseif connectedNodeIndex == -1 
                        flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex);  
                    end
                end
                waterVolume_node(i) = obj.Nodes{i}.waterCrossSectionArea /obj.Nodes{i}.area * ...
                    obj.Nodes{i}.volume;
                residenceTime_node(i) = waterVolume_node(i) / flowRate_node(i);
            end
            
            %Set TimeStep
            linkTime = min(nonzeros(residenceTime_link));
            nodeTime = min(nonzeros(residenceTime_node));
            timeStep = min(linkTime,nodeTime)/2;            
            
            B = zeros(obj.numberOfLinks,1);            
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            Known = zeros(obj.numberOfNodes, 1);
            
            % calculation of 6 Unknowns (concentration & adsorbedConcentrations of nodes & 
            % concentration & adsorbedConcentrations of links) in each timeStep 
            
            t = 0; 
            time = 0;
            timePlot = zeros(round(poreVolumeInjected/timeStep)+1 ,1);
            flux_averagedConcentration = zeros(round(poreVolumeInjected/timeStep)+1 ,1);
            
            while time < poreVolumeInjected
                t = t+1
                time = time + timeStep; 
                timePlot(t) = time;
                sumOfConcentration = 0;
                sumOfFlowRate = 0;
                
            % calculate concentration of nodes: based on eq.10
            for i = 1:obj.numberOfNodes 
                
                adsorbedConcentration_SolidFluid_node = ...
                    obj.Nodes{i}.adsorbedConcentration_SolidFluid(t);
                adsorbedConcentration_FluidFluid_node = ...
                    obj.Nodes{i}.adsorbedConcentration_FluidFluid(t);
                J = alpha_node(i) * timeStep / (1+alpha_node(i) * timeStep);
                E = 1 + timeStep * flowRate_node(i) / waterVolume_node(i) + ...
                    2 * alpha_node(i) * Kd_node(i) * timeStep - ...
                    2 * (alpha_node(i) * timeStep) ^ 2 * Kd_node(i) / ...
                    (1 + timeStep * flowRate_node(i));
                I = timeStep / waterVolume_node(i);
                
                for j = 1:obj.Nodes{i}.connectionNumber 
                    
                    connectedLinkIndex = obj.Nodes{i}.connectedLinks(j);
                    connectedNodeIndex = obj.Nodes{i}.connectedNodes(j);
                    
                    adsorbedConcentration_SolidFluid_link = ...
                        obj.Links{connectedLinkIndex}.adsorbedConcentration_SolidFluid(t);
                    adsorbedConcentration_FluidFluid_link = ...
                        obj.Links{connectedLinkIndex}.adsorbedConcentration_FluidFluid(t);
                    linksConcentration = obj.Links{connectedLinkIndex}.concentration(t);
                    
                    B(connectedLinkIndex) = 1 + ...
                        (flowRate_link(connectedLinkIndex) * timeStep /...
                        waterVolume_link(connectedLinkIndex)) + ...
                        2*(timeStep * alpha_link(connectedLinkIndex) * Kd_link(connectedLinkIndex)) - ...
                        2*(timeStep^2 * (alpha_link(connectedLinkIndex))^2 * Kd_link(connectedLinkIndex))/...
                        (1 + timeStep * alpha_link(connectedLinkIndex));
                    F = 1 / B(connectedLinkIndex) * ...
                        (flowRate_link(connectedLinkIndex))^2 * timeStep / waterVolume_link(connectedLinkIndex);
                    G = flowRate_link(connectedLinkIndex) /  B(connectedLinkIndex);
                    H = flowRate_link(connectedLinkIndex) * timeStep * alpha_link(connectedLinkIndex)/...
                        B(connectedLinkIndex)/(1 + timeStep * alpha_link(connectedLinkIndex));
                                
                    if connectedNodeIndex ~= 0 && connectedNodeIndex ~= -1
                        
                        % determine link flowing into this node
                        if obj.Nodes{connectedNodeIndex}.waterPressure > obj.Nodes{i}.waterPressure        
                            
                            Factor(i, connectedNodeIndex) = -I * F;                    
                            Known(i,1) = Known(i,1) + G * linksConcentration +...
                                H * (adsorbedConcentration_SolidFluid_link + ...
                                adsorbedConcentration_FluidFluid_link);
                        end
                    elseif connectedNodeIndex == -1 
                        
                        Known(i,1) = Known(i,1) + G * linksConcentration + ...
                            H * (adsorbedConcentration_SolidFluid_link + ...
                            adsorbedConcentration_FluidFluid_link) + ...
                            F * soluteConcentration;                    
                    end
                end
                Factor(i, i) = E;
                Known(i,1) = obj.Nodes{i}.concentration(t) + I * Known(i,1) + ...
                    J * (adsorbedConcentration_SolidFluid_node + ...
                    adsorbedConcentration_FluidFluid_node);                    
            end
            
            nodesConcentration_new = pcg(Factor, Known, 1e-10, obj.numberOfNodes);  
            
            % asign new concentration of nodes & 
            % calculate new adsorbedConcentrations of nodes:
            for i = 1:obj.numberOfNodes
                if nodesConcentration_new(i) > soluteConcentration
                obj.Nodes{i}.concentration(t+1) = soluteConcentration;
                else
                obj.Nodes{i}.concentration(t+1) = nodesConcentration_new(i);
                end 
                
                obj.Nodes{i}.adsorbedConcentration_SolidFluid(t+1) = ...
                    (alpha_node(i) * timeStep * Kd_node(i) * ...
                    obj.Nodes{i}.concentration(t+1)+ ...
                    obj.Nodes{i}.adsorbedConcentration_SolidFluid(t))/...
                    (1 + alpha_node(i) * timeStep);
                
                obj.Nodes{i}.adsorbedConcentration_FluidFluid(t+1) = ...
                    (alpha_node(i) * timeStep * Kd_node(i) * ...
                    obj.Nodes{i}.concentration(t+1)+ ...
                    obj.Nodes{i}.adsorbedConcentration_FluidFluid(t))/...
                    (1 + alpha_node(i) * timeStep);
                
            end 
            
            % calculate new concentration & adsorbedConcentrations of links: 
            for i = 1:obj.numberOfLinks 
                
                node1Index = obj.Links{i}.pore1Index;
                node2Index = obj.Links{i}.pore2Index;
                
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet
                    if obj.Nodes{node1Index}.waterPressure > obj.Nodes{node2Index}.waterPressure
                        upstreamNode = node1Index;
                    else
                        upstreamNode = node2Index;
                    end
                    obj.Links{i}.concentration(t+1) = obj.Links{i}.concentration(t)/ B(i) + ...
                        timeStep * alpha_link(i) * (obj.Links{i}.adsorbedConcentration_SolidFluid(t)+...
                        obj.Links{i}.adsorbedConcentration_FluidFluid(t))/ ...
                        (B(i)*(1 + alpha_link(i) * timeStep)) + ...
                        flowRate_link(i) * timeStep * obj.Nodes{upstreamNode}.concentration(t+1);                      
                elseif obj.Links{i}.isInlet
                    obj.Links{i}.concentration(t+1) = obj.Links{i}.concentration(t)/ B(i) + ...
                        timeStep * alpha_link(i) * (obj.Links{i}.adsorbedConcentration_SolidFluid(t)+...
                        obj.Links{i}.adsorbedConcentration_FluidFluid(t)) / ...
                        (B(i)*(1 + alpha_link(i) * timeStep)) + ...
                        flowRate_link(i) * timeStep * soluteConcentration;
                else
                    obj.Links{i}.concentration(t+1) = obj.Links{i}.concentration(t)/ B(i) + ...
                        timeStep * alpha_link(i) * (obj.Links{i}.adsorbedConcentration_SolidFluid(t)+...
                        obj.Links{i}.adsorbedConcentration_FluidFluid(t)) / ...
                        (B(i)*(1 + alpha_link(i) * timeStep)) + ...
                        flowRate_link(i) * timeStep * obj.Nodes{node1Index}.concentration(t+1);
                    % calculation for BreakThroughCurve at outlet of network
                    sumOfConcentration = sumOfConcentration + ...
                        obj.Nodes{node1Index}.concentration(t)*flowRate_node(node1Index);
                    sumOfFlowRate = sumOfFlowRate + flowRate_node(node1Index);
                end
                obj.Links{i}.adsorbedConcentration_SolidFluid(t+1) = ...
                    (alpha_link(i) * timeStep * Kd_link(i) * ...
                    obj.Links{i}.concentration(t+1)+ ...
                    obj.Links{i}.adsorbedConcentration_SolidFluid(t))/...
                    (1 + alpha_link(i) * timeStep);
                
                obj.Links{i}.adsorbedConcentration_FluidFluid(t+1) = ...
                    (alpha_link(i) * timeStep * Kd_link(i) * ...
                    obj.Links{i}.concentration(t+1)+ ...
                    obj.Links{i}.adsorbedConcentration_FluidFluid(t))/...
                    (1 + alpha_link(i) * timeStep);
            end
            
            % calculate BreakThroughCurve at outlet of network
            flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration;
            
            end 
            plot(flux_averagedConcentration,timePlot,'*');
        end
        
        %% vtk file generation
        function fileName = vtkWriter(obj,process,i)
            num = num2str(i, '%0.3d');
            fileName = strcat(process,num,'.vtk');
            vtkFileID = fopen(fileName,'w');
            if vtkFileID == -1
                error('Cannot open file for writing.');
            end
            title = fileName;
            inletLinks = obj.numberOfLinks - obj.numOfInletLinks - obj.numOfOutletLinks;
            points = 9 * obj.numOfSquarePores + 7 * (obj.numOfCircularPores + obj.numOfTriangularPores);
            surfaces = 6 * obj.numOfSquarePores + 8 * obj.numOfCircularPores + 5 * obj.numOfTriangularPores;
            subSurfaces = 30 * obj.numOfSquarePores + 32 * obj.numOfCircularPores + 23 * obj.numOfTriangularPores; 
            fprintf ( vtkFileID, '# vtk DataFile Version 2.0\n' );
            fprintf ( vtkFileID, '%s\n', title );
            fprintf ( vtkFileID, 'ASCII\n' );
            fprintf ( vtkFileID, '\n' );
            fprintf ( vtkFileID, 'DATASET POLYDATA\n' );
            
            fprintf ( vtkFileID, 'POINTS %d double\n', points );
            % pore coordination
            for i = 1:obj.numberOfNodes
                fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate );
            end
            % corner points of pores
            for i = 1:obj.numberOfNodes
                    r = obj.Nodes{i}.radius;
                if strcmp(obj.Nodes{i}.geometry , 'Square')== 1 
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-r, obj.Nodes{i}.y_coordinate-r, obj.Nodes{i}.z_coordinate-r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+r, obj.Nodes{i}.y_coordinate-r, obj.Nodes{i}.z_coordinate-r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+r, obj.Nodes{i}.y_coordinate+r, obj.Nodes{i}.z_coordinate-r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-r, obj.Nodes{i}.y_coordinate+r, obj.Nodes{i}.z_coordinate-r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-r, obj.Nodes{i}.y_coordinate-r, obj.Nodes{i}.z_coordinate+r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+r, obj.Nodes{i}.y_coordinate-r, obj.Nodes{i}.z_coordinate+r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+r, obj.Nodes{i}.y_coordinate+r, obj.Nodes{i}.z_coordinate+r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-r, obj.Nodes{i}.y_coordinate+r, obj.Nodes{i}.z_coordinate+r );
                elseif strcmp(obj.Nodes{i}.geometry , 'Circle')== 1
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-r, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate);
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate, obj.Nodes{i}.y_coordinate+r, obj.Nodes{i}.z_coordinate );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate-r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate, obj.Nodes{i}.y_coordinate-r, obj.Nodes{i}.z_coordinate );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate+r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+r, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate);
                else
                    y1 = r;%/tan(obj.Nodes{i}.halfAngle1);
                    y2 = r;%/tan(obj.Nodes{i}.halfAngle2);
                    z1 = r;%/sin(obj.Nodes{i}.halfAngle3);
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-r, obj.Nodes{i}.y_coordinate+y1, obj.Nodes{i}.z_coordinate-r);
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-r, obj.Nodes{i}.y_coordinate-y2, obj.Nodes{i}.z_coordinate-r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-r, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate+z1 );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+r, obj.Nodes{i}.y_coordinate+y1, obj.Nodes{i}.z_coordinate-r);
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+r, obj.Nodes{i}.y_coordinate-y2, obj.Nodes{i}.z_coordinate-r );
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+r, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate+z1 );                
                end
            end  
            
            fprintf ( vtkFileID, 'LINES %d %d\n', inletLinks, 3*inletLinks);      
            for i = 1:obj.numberOfLinks
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet
                    pore1Index = obj.Links{i}.pore2Index-1;
                    pore2Index = obj.Links{i}.pore1Index-1;
                    fprintf( vtkFileID,'%d %d %d %d\n', 2, pore2Index,pore1Index);  
                end
            end 
            
            fprintf ( vtkFileID, 'POLYGONS %d %d\n', surfaces, subSurfaces);
            previouspoints = obj.numberOfNodes;
            for i = 1:obj.numberOfNodes 
                
                if strcmp(obj.Nodes{i}.geometry , 'Square')== 1
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 0+previouspoints,1+previouspoints,2+previouspoints,3+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 4+previouspoints,5+previouspoints,6+previouspoints,7+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 0+previouspoints,1+previouspoints,5+previouspoints,4+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 2+previouspoints,3+previouspoints,7+previouspoints,6+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 0+previouspoints,4+previouspoints,7+previouspoints,3+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 1+previouspoints,2+previouspoints,6+previouspoints,5+previouspoints);                    
                    previouspoints = previouspoints + 8;
                elseif strcmp(obj.Nodes{i}.geometry , 'Circle')== 1
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 0+previouspoints,1+previouspoints,2+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 0+previouspoints,2+previouspoints,3+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 0+previouspoints,3+previouspoints,4+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 0+previouspoints,4+previouspoints,1+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 5+previouspoints,1+previouspoints,2+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 5+previouspoints,2+previouspoints,3+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 5+previouspoints,3+previouspoints,4+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 5+previouspoints,4+previouspoints,1+previouspoints);                  
                    previouspoints = previouspoints + 6;
                else
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 0+previouspoints,1+previouspoints,2+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 0+previouspoints,1+previouspoints,4+previouspoints,3+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 1+previouspoints,2+previouspoints,5+previouspoints,4+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d %d\n', 4, 2+previouspoints,0+previouspoints,3+previouspoints,5+previouspoints);
                    fprintf( vtkFileID,'%d %d %d %d\n', 3, 3+previouspoints,4+previouspoints,5+previouspoints);                  
                    previouspoints = previouspoints + 6;
                end
            end       
            
            fprintf ( vtkFileID, 'CELL_DATA  %d \n',1*( surfaces+inletLinks));
            
            % Initializing
            if strcmp(process , 'init')== 1 
            fprintf ( vtkFileID, 'SCALARS pressure float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            % Pressure distribution
            for i = 1:obj.numberOfLinks 
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet  
                        fprintf( vtkFileID,'%3.5f \n', obj.Links{i}.waterPressure);  
                end
            end 
            for i = 1:obj.numberOfNodes
                for j = 1:obj.Nodes{i}.surfaces
                fprintf( vtkFileID,'%3.5f \n', obj.Nodes{i}.waterPressure);
                end 
            end 
            
            % Diffusion
            elseif strcmp(process , 'Diff')== 1 
            fprintf ( vtkFileID, 'SCALARS concentration float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes
                for j = 1:obj.Nodes{i}.surfaces
                fprintf( vtkFileID,'%3.5f \n', obj.Nodes{i}.concentration);
                end 
            end 
            
            % Primary Drainage
            elseif strcmp(process , 'PD')== 1
            fprintf ( vtkFileID, 'SCALARS occupancy float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfLinks
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet
                fprintf( vtkFileID,'%3.5f \n', obj.Links{i}.occupancy); 
                end
            end 
            for i = 1:obj.numberOfNodes
                for j = 1:obj.Nodes{i}.surfaces
                fprintf( vtkFileID,'%3.5f \n', obj.Nodes{i}.occupancy);
                end 
            end 
%             fprintf ( vtkFileID, 'SCALARS waterSaturation float  %d \n', 1);
%             fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
%             for i = 1:obj.numberOfLinks
%                 fprintf( vtkFileID,'%3.5f \n', obj.Links{i}.waterSaturation); 
%             end 
%             for i = 1:obj.numberOfNodes
%                 for j = 1:obj.Nodes{i}.surfaces
%                 fprintf( vtkFileID,'%3.5f \n', obj.Nodes{i}.waterSaturation);
%                 end 
%             end 
            end
                
        end  
        
        function fileName = vtkWriter_glyph(obj,process,i)
            num = num2str(i, '%0.3d');
            fileName = strcat(process,'_glyph',num,'.vtk');
            vtkFileID = fopen(fileName,'w');
            if vtkFileID == -1
                error('Cannot open file for writing.');
            end
            title = fileName;
            imaginaryPoints = obj.numOfInletLinks + obj.numOfOutletLinks; 
            points = obj.numberOfNodes+imaginaryPoints; 
            fprintf ( vtkFileID, '# vtk DataFile Version 2.0\n' );
            fprintf ( vtkFileID, '%s\n', title );
            fprintf ( vtkFileID, 'ASCII\n' );
            fprintf ( vtkFileID, '\n' );
            fprintf ( vtkFileID, 'DATASET POLYDATA\n' );
            
            fprintf ( vtkFileID, 'POINTS %d double\n', points );
            for i = 1:obj.numberOfNodes 
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate); 
            end    
            
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate-3*obj.Nodes{i}.radius, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate);
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate+3*obj.Nodes{i}.radius, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate);
                end
            end    
            
            fprintf ( vtkFileID, 'LINES %d %d\n', obj.numberOfLinks, 3*obj.numberOfLinks);   
            imagine = -1;
            for i = 1:obj.numberOfLinks
                    pore1Index = obj.Links{i}.pore1Index-1;
                    pore2Index = obj.Links{i}.pore2Index-1;
                if ~obj.Links{i}.isInlet && ~obj.Links{i}.isOutlet
                    fprintf( vtkFileID,'%d %d %d %d\n', 2, pore2Index,pore1Index); 
                elseif obj.Links{i}.isInlet 
                    imagine = imagine + 1;
                    fprintf( vtkFileID,'%d %d %d %d\n', 2, pore2Index,obj.numberOfNodes+imagine); 
                else
                    imagine = imagine + 1;
                    fprintf( vtkFileID,'%d %d %d %d\n', 2, obj.numberOfNodes+imagine,pore1Index); 
                end
            end 
            
            % Initializing
            if strcmp(process , 'init')== 1 
            fprintf ( vtkFileID, 'POINT_DATA  %d \n', obj.numberOfNodes+imaginaryPoints);
            fprintf ( vtkFileID, 'SCALARS radius float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes 
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius); 
            end 
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius*0.001);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius*0.001); 
                end
            end   
            fprintf ( vtkFileID, 'SCALARS pressure float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes 
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.waterPressure); 
            end 
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', 1);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 0); 
                end
            end   
            
            fprintf ( vtkFileID, 'CELL_DATA  %d \n',obj.numberOfLinks);
            fprintf ( vtkFileID, 'SCALARS pressure float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n'); 
            for i = 1:obj.numberOfLinks  
                fprintf( vtkFileID,'%3.8f \n', obj.Links{i}.waterPressure);
            end  
            
            % Diffusion
            elseif strcmp(process , 'Diff')== 1 
            fprintf ( vtkFileID, 'POINT_DATA  %d \n', obj.numberOfNodes+imaginaryPoints);
            fprintf ( vtkFileID, 'SCALARS radius float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes 
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius); 
            end 
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius*0.001);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius*0.001); 
                end
            end  
            fprintf ( vtkFileID, 'SCALARS pressure float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes 
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.waterPressure); 
            end 
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', 1);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 0); 
                end
            end  
            fprintf ( vtkFileID, 'SCALARS concentration float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes 
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.concentration); 
            end
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', 1);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 0); 
                end
            end  
            
            fprintf ( vtkFileID, 'CELL_DATA  %d \n',obj.numberOfLinks );
            fprintf ( vtkFileID, 'SCALARS radius float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n'); 
            for i = 1:obj.numberOfLinks  
                fprintf( vtkFileID,'%3.8f \n', obj.Links{i}.radius);
            end  
            fprintf ( vtkFileID, 'SCALARS pressure float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n'); 
            for i = 1:obj.numberOfLinks  
                fprintf( vtkFileID,'%3.8f \n', obj.Links{i}.waterPressure);
            end    
            fprintf ( vtkFileID, 'SCALARS concentration float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n'); 
            for i = 1:obj.numberOfLinks  
                fprintf( vtkFileID,'%3.8f \n', obj.Links{i}.concentration);
            end   
            % Primary Drainage
            elseif strcmp(process , 'PD')== 1 || strcmp(process , 'SI')== 1 
            fprintf ( vtkFileID, 'POINT_DATA  %d \n', obj.numberOfNodes+imaginaryPoints);
            fprintf ( vtkFileID, 'SCALARS radius float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes 
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius); 
            end 
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius*0.001);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.radius*0.001); 
                end
            end  
            fprintf ( vtkFileID, 'SCALARS pressure float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes 
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.waterPressure); 
            end 
            for i = 1:obj.numberOfNodes
                if strcmp(process , 'PD')== 1
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', 0);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 1); 
                end
                elseif strcmp(process , 'SI')== 1 
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', 1);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 0); 
                end
                end
            end  
            fprintf ( vtkFileID, 'SCALARS waterSaturation float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n');
            for i = 1:obj.numberOfNodes 
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.waterSaturation); 
            end
            for i = 1:obj.numberOfNodes
                if strcmp(process , 'PD')== 1
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', 0);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 1); 
                end
                elseif strcmp(process , 'SI')== 1 
                if obj.Nodes{i}.isInlet 
                    fprintf( vtkFileID,'%3.8f \n', 1);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 0); 
                end
                end
            end  
            
            fprintf ( vtkFileID, 'CELL_DATA  %d \n',obj.numberOfLinks);
            fprintf ( vtkFileID, 'SCALARS radius float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n'); 
            for i = 1:obj.numberOfLinks   
                fprintf( vtkFileID,'%3.8f \n', obj.Links{i}.radius);
            end  
            fprintf ( vtkFileID, 'SCALARS pressure float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n'); 
            for i = 1:obj.numberOfLinks   
                fprintf( vtkFileID,'%3.8f \n', obj.Links{i}.waterPressure); 
            end    
            fprintf ( vtkFileID, 'SCALARS waterSaturation float  %d \n', 1);
            fprintf ( vtkFileID, 'LOOKUP_TABLE default\n'); 
            for i = 1:obj.numberOfLinks  
                fprintf( vtkFileID,'%3.8f \n', obj.Links{i}.waterSaturation);
            end   
            end
                
        end 
    end
end

