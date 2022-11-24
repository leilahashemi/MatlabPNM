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
        function calculateReactiveTransport_SinglePhaseDiffusion_R(obj, inletPressure, outletPressure, soluteConcentration, poreVolumeInjected) 
            
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
            timeStep = timeStep * 10;
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
            
              obj.vtkWriter_glyph('Diff',t);
              
%             movefile nameFile D:\Uni\PUT\Dissertation\PNM_Code\MatlabPNM\Ver.2\seq\Diff
            % Plot & Animation
            addpoints(h,timePlot(t),flux_averagedConcentration(t));  
            % GIF 
            plot(timePlot(1:t),flux_averagedConcentration(1:t), 'b','LineWidth',2);            
            axis([0 simulationTime 0 1])
            title('Break Through Curve')
            xlabel('Time(s)')            
            ylabel('DimensionlessConcentration(-)') 
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
        
        % Single phase Reactive transport_Edgar         
        function calculateReactiveTransport_SinglePhaseDiffusion(obj, inletPressure, outletPressure, soluteConcentration, simulationVolume, poreVolumeInjected) 
            
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
%                 obj.Links{ii}.flow = flowRate_link(ii)*10^14;
                residenceTime_link(ii) = obj.Links{ii}.volume/flowRate_link(ii);
                diffusion_link (ii) = effectiveDiffusion * obj.Links{ii}.area / obj.Links{ii}.linkLength;
            end
            obj.vtkWriter_glyph('init',0);
            
            timeStep = min(nonzeros(residenceTime_link)); 
%             timeStep = timeStep /2; 
            
            % calculation of 3 Unknowns (concentration of nodes & links) in each timeStep 
             
            t = 0; 
            time = 0;
            simulationTime = simulationVolume / obj.totalFlowRate;
            injectionTime = poreVolumeInjected / obj.totalFlowRate;
            
            timePlot = zeros(round(simulationTime/timeStep)+1 ,1);
            flux_averagedConcentration = zeros(round(simulationTime/timeStep)+1 ,1);
            obj.BreakThroughCurve_singlePhase = zeros(round(simulationTime/timeStep)+1 ,2);
            soluteConcentration1 = soluteConcentration;
                  
            while time < simulationTime
                
                if time > injectionTime
                    soluteConcentration = 0;
                end
                
                t = t+1;
                time = time + timeStep;
                timePlot(t) = time;
                sumOfConcentration = 0;
                sumOfFlowRate = 0;
                 
                diffusion_node = zeros(obj.numberOfNodes,1);
                Factor = zeros(obj.numberOfNodes + obj.numberOfLinks, obj.numberOfNodes + obj.numberOfLinks);
                Known = zeros(obj.numberOfNodes + obj.numberOfLinks, 1); 
                % calculate concentration of nodes & links based on eq. 8 & 9
                for i = 1:obj.numberOfNodes + obj.numberOfLinks
                    if i <= obj.numberOfNodes
                        for j = 1:obj.Nodes{i}.connectionNumber 
                            
                            connectedLinkIndex = obj.Nodes{i}.connectedLinks(j);
                            connectedNodeIndex = obj.Nodes{i}.connectedNodes(j);
                            jj = obj.numberOfNodes + connectedLinkIndex;
                            
                            diffusion_node(i) = diffusion_node(i) + diffusion_link(connectedLinkIndex);                             
                            Factor(i, jj) = -timeStep / obj.Nodes{i}.volume * diffusion_link(connectedLinkIndex);
                            
                            % determine link flowing into this node
                            if connectedNodeIndex ~= 0 && connectedNodeIndex ~= -1
                                
                                if obj.Nodes{connectedNodeIndex}.waterPressure > obj.Nodes{i}.waterPressure
                                    Factor(i, i) = Factor(i, i) + timeStep / obj.Nodes{i}.volume * flowRate_link(connectedLinkIndex);
                                    Factor(i, jj) = Factor(i, jj) - timeStep / obj.Nodes{i}.volume * flowRate_link(connectedLinkIndex); 
                                end                                
                                Factor(i, i) = Factor(i, i) + timeStep / obj.Nodes{i}.volume * diffusion_link(connectedLinkIndex);
                            elseif connectedNodeIndex == -1
                                Factor(i, i) = Factor(i, i) + ...
                                    timeStep / obj.Nodes{i}.volume * (flowRate_link(connectedLinkIndex) + diffusion_link(connectedLinkIndex));
                                Factor(i, jj) = Factor(i, jj) - timeStep / obj.Nodes{i}.volume * flowRate_link(connectedLinkIndex);  
                            end
                        end
                        
                        Factor(i, i) = 1 + Factor(i, i);
                        Known(i,1) = obj.Nodes{i}.concentration(t);
                    else
                        jj = i - obj.numberOfNodes;
                        node1Index = obj.Links{jj}.pore1Index;
                        node2Index = obj.Links{jj}.pore2Index;                        
                        Known(i,1) = obj.Links{jj}.concentration(t); 
                        Factor(i, i) = 1 + timeStep  / obj.Links{jj}.volume * (flowRate_link(jj)+ diffusion_link(jj)) ;
                        
                        if ~obj.Links{jj}.isInlet && ~obj.Links{jj}.isOutlet
                            if obj.Nodes{node1Index}.waterPressure > obj.Nodes{node2Index}.waterPressure
                                Factor(i, node1Index) = - timeStep  / obj.Links{jj}.volume * ...
                                    (flowRate_link(jj)+ diffusion_link(jj)) ;
                                Factor(i, node2Index) = - timeStep  / obj.Links{jj}.volume * ...
                                    (diffusion_link(jj)) ;
                            elseif obj.Nodes{node2Index}.waterPressure > obj.Nodes{node1Index}.waterPressure
                                Factor(i, node2Index) = - timeStep  / obj.Links{jj}.volume * ...
                                    (flowRate_link(jj)+ diffusion_link(jj)) ;
                                Factor(i, node1Index) = - timeStep  / obj.Links{jj}.volume * ...
                                    (diffusion_link(jj)) ;
                            end
                        elseif obj.Links{jj}.isInlet
                            Factor(i, node2Index) = - timeStep  / obj.Links{jj}.volume * ...
                                (diffusion_link(jj)) ;
                            Known(i,1) = Known(i,1) + soluteConcentration * timeStep  / obj.Links{jj}.volume * ...
                                    (flowRate_link(jj)+ diffusion_link(jj));
                        elseif obj.Links{jj}.isOutlet 
                            Factor(i, i) = 0 ; 
                            Known(i,1) = Known(i,1) - obj.Links{jj}.concentration(t);
                        end
                    end
                end
                
                nodesConcentration_new = gmres (Factor, Known,[], 1e-10, obj.numberOfNodes + obj.numberOfLinks);
                
                % asign new concentration of nodes
                for i = 1:obj.numberOfNodes
                    if nodesConcentration_new(i) > soluteConcentration1
                        obj.Nodes{i}.concentration(t+1) = soluteConcentration1;
                    else
                        obj.Nodes{i}.concentration(t+1) = nodesConcentration_new(i);
                    end
                end
                
                % calculate new concentration of links
                for i = 1:obj.numberOfLinks  
                    jj = i + obj.numberOfNodes;
                    if nodesConcentration_new(jj) > soluteConcentration1
                        obj.Links{i}.concentration(t+1) = soluteConcentration1;
                    else
                        obj.Links{i}.concentration(t+1) = nodesConcentration_new(jj);
                    end
                    if obj.Links{i}.isOutlet 
                        obj.Links{i}.concentration(t+1) = obj.Nodes{obj.Links{ii}.pore1Index}.concentration(t+1);
                        sumOfConcentration = sumOfConcentration + ...
                            obj.Links{i}.concentration(t)*flowRate_link(i);
                        sumOfFlowRate = sumOfFlowRate + flowRate_link(i);
                    end
                end
                % calculate BreakThroughCurve at outlet of network
                flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration1;
                
            obj.BreakThroughCurve_singlePhase(t,1) = timePlot(t);
            obj.BreakThroughCurve_singlePhase(t,2) = flux_averagedConcentration(t);
%             obj.vtkWriter_glyph('Diff',t); 
            end
             
%             obj.BreakThroughCurve_singlePhase(:,1) = timePlot;
%             obj.BreakThroughCurve_singlePhase(:,2) = flux_averagedConcentration;
            
            % Plot   
            plot(timePlot,flux_averagedConcentration,'*');
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
                t = t+1;
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
         
        %% vtk file generation 
        function fileName = vtkWriter_glyph(obj,process,ii)
            num = num2str(ii, '%0.3d');
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
                    fprintf( vtkFileID,'%3.8f \n', 0);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 0); 
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
                    fprintf( vtkFileID,'%3.8f \n',0);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 0); 
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
                fprintf( vtkFileID,'%3.8f \n', obj.Nodes{i}.concentration(ii)); 
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
                fprintf( vtkFileID,'%3.8f \n', obj.Links{i}.concentration(ii));
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
                    fprintf( vtkFileID,'%3.8f \n', 0);                     
                elseif obj.Nodes{i}.isOutlet
                    fprintf( vtkFileID,'%3.8f \n', 0); 
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

