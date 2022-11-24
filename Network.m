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
        maxCoordinationNumber  
        
        Porosity
        poreVolume
        networkVolume
        absolutePermeability
        ThSD
        PSD
        
        totalFlowRate
        totalFlowRate_o
        totalFlowRate_w
        velocity
        capillaryNumber
        pecletNumber 
        Pc_drain_max
        
        DrainageData
        ImbibitionData 
        BreakThroughCurve_singlePhase
        BreakThroughCurve_twoPhase
    end    
    methods
        %% Cunstructor function
        function obj = Network(fileName) 
            
            % Opening the files
            link_1_fileID = fopen(strcat(fileName, '_link1.dat'));
            obj.numberOfLinks = str2num(fgetl(link_1_fileID));
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
        
        %% Network Properties calculation
        function calculateNetworkProperties(obj, inletPressure, outletPressure)
            obj.ThSD = zeros(obj.numberOfLinks,1);
            obj.PSD = zeros(obj.numberOfNodes,1);
            obj.numOfInletLinks = 0;
            obj.numOfOutletLinks = 0; 
            nodesVolume = 0;
            linksVolume = 0;

            for ii = 1:obj.numberOfNodes
                nodesVolume = nodesVolume + (obj.Nodes{ii}.volume); 
                obj.maxCoordinationNumber(ii,1) = obj.Nodes{ii}.connectionNumber; 
                obj.PSD(ii,1) = 2 * obj.Nodes{ii}.radius;
%                 %Isolated element
%                 if obj.Nodes{ii}.connectionNumber == 0
%                     ii
%                 end
            end 
            obj.maxCoordinationNumber = max(obj.maxCoordinationNumber); 
%             hist(obj.PSD);
            
            for ii = 1:obj.numberOfLinks 
                linksVolume = linksVolume + (obj.Links{ii}.volume);
                obj.ThSD (ii,1)= 2 * obj.Links{ii}.radius;
                if obj.Links{ii}.isInlet
                    obj.numOfInletLinks = obj.numOfInletLinks + 1;
                elseif obj.Links{ii}.isOutlet
                    obj.numOfOutletLinks = obj.numOfOutletLinks+1;                 
                end
            end 
%             hist(obj.ThSD);
            obj.networkVolume = obj.xDimension * obj.yDimension * obj.zDimension;
            obj.poreVolume = linksVolume + nodesVolume;
            obj.Porosity = obj.poreVolume * 100 / (obj.xDimension * obj.yDimension * obj.zDimension);      
            calculateAbsolutePermeability(obj, inletPressure, outletPressure)
        end 
        
        %% Pressure distribution calculation of single phase flow       
        function pressureDistribution_singlePhaseFlow (obj, inletPressure, outletPressure)
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            B = zeros(obj.numberOfNodes, 1);
     
            for ii = 1:obj.numberOfLinks
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;

                % if the link is connected to inlet (index of node 1 is -1 which does not exist) 
                if obj.Links{ii}.isInlet
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance) +...
                        ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance)))^-1;
                    
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + nodeLinkSystemConductance;
                    B(node2Index) = nodeLinkSystemConductance * inletPressure;
                    
                % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                elseif obj.Links{ii}.isOutlet
                     nodeLinkSystemConductance = ( (obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance)))^-1;
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + nodeLinkSystemConductance;
                    B(node1Index) = nodeLinkSystemConductance * outletPressure;
                    
                %if the link is neither inlet nor outlet    
                else
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance) +...
                        (obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance)))^-1;   
                
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + nodeLinkSystemConductance;
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + nodeLinkSystemConductance;
                    Factor(node1Index, node2Index) = Factor(node1Index, node2Index) - nodeLinkSystemConductance;
                    Factor(node2Index, node1Index) = Factor(node2Index, node1Index) - nodeLinkSystemConductance;
                   
                end     
            end
            
            % using GMRES method to solve the pressure distribution             
            nodesPressure = gmres(Factor, B,[], 1e-10, 1000);
%             nodesPressure = pcg(Factor, B, 1e-7, 1000);
            %assign the pressure values to each node
             x_coor = zeros(obj.numberOfNodes,1); 
            for ii = 1:obj.numberOfNodes
                if nodesPressure(ii) > inletPressure
                    obj.Nodes{ii}.waterPressure = inletPressure; 
                elseif nodesPressure(ii) < outletPressure
                    obj.Nodes{ii}.waterPressure = outletPressure; 
                else
                    obj.Nodes{ii}.waterPressure = nodesPressure(ii); 
                end                
                x_coor(ii,1) = obj.Nodes{ii}.x_coordinate;
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
           
            % Plot Pressure
            % a , b are 2 surfaces perpendicular to x-direction with
            % distance equals to intervalx
            x_outlet = max(x_coor);
            x_inlet = min(x_coor);
            n = 100;
            intervalx = (x_outlet - x_inlet)/n;
            a = x_inlet;
            b = x_inlet + intervalx;
            x = zeros(n,1);
            press_x = zeros(100,1);
            for i = 1:n
                area = 0;
                for ii = 1:obj.numberOfNodes
                    if obj.Nodes{ii}.x_coordinate >= a && obj.Nodes{ii}.x_coordinate < b
%                        press_x(i) = press_x(i) + obj.Nodes{ii}.waterPressure*obj.Nodes{ii}.area;
%                        area= area+obj.Nodes{ii}.area;                       
                       press_x(i) = press_x(i) + obj.Nodes{ii}.waterPressure;
                       area= area+1;                       
                    end
                end
                for ii = 1:obj.numberOfLinks    
                    if ~obj.Links{ii}.isOutlet
                        if obj.Nodes{obj.Links{ii}.pore2Index}.x_coordinate >= a && obj.Nodes{obj.Links{ii}.pore2Index}.x_coordinate < b
%                             press_x(i) = press_x(i) + obj.Links{ii}.waterPressure*obj.Links{ii}.area;
%                             area= area+obj.Links{ii}.area;                            
                            press_x(i) = press_x(i) + obj.Links{ii}.waterPressure;
                            area = area + 1;
                        end
                    end
                end
                press_x(i)=press_x(i)/area;
                x(i) = x(i) + i*intervalx;                
                a = a + intervalx;
                b = b + intervalx;
            end                  
%             plot(x, press_x, '*')
%             title('Pressure drop in x-direction')
%             xlabel('X(m)')
%             xlim([x_inlet x_outlet])
%             ylabel('Pressure(Pa)') 
        end 
        function pressureDistribution_singlePhaseFlow_o (obj, inletPressure, outletPressure)
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            B = zeros(obj.numberOfNodes, 1);
            cond_factor = obj.waterViscosity / obj.oilViscosity;
     
            for ii = 1:obj.numberOfLinks
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;

                % if the link is connected to inlet (index of node 1 is -1 which does not exist) 
                if obj.Links{ii}.isInlet
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance*cond_factor) +...
                        ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance*cond_factor)))^-1;
                    
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + nodeLinkSystemConductance;
                    B(node2Index) = nodeLinkSystemConductance * inletPressure;
                    
                % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                elseif obj.Links{ii}.isOutlet
                     nodeLinkSystemConductance = ( (obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance*cond_factor) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance*cond_factor)))^-1;
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + nodeLinkSystemConductance;
                    B(node1Index) = nodeLinkSystemConductance * outletPressure;
                    
                %if the link is neither inlet nor outlet    
                else
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance*cond_factor) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance*cond_factor) +...
                        (obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance*cond_factor)))^-1;   
                
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + nodeLinkSystemConductance;
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + nodeLinkSystemConductance;
                    Factor(node1Index, node2Index) = Factor(node1Index, node2Index) - nodeLinkSystemConductance;
                    Factor(node2Index, node1Index) = Factor(node2Index, node1Index) - nodeLinkSystemConductance;
                   
                end     
            end
            
            % using GMRES method to solve the pressure distribution 
            nodesPressure = gmres(Factor, B,[], 1e-7, 1000);
%             nodesPressure = pcg(Factor, B, 1e-7, 1000);
            %assign the pressure values to each node
             x_coor = zeros(obj.numberOfNodes,1); 
            for ii = 1:obj.numberOfNodes
                if nodesPressure(ii) > inletPressure
                    obj.Nodes{ii}.oilPressure = inletPressure; 
                elseif nodesPressure(ii) < outletPressure
                    obj.Nodes{ii}.oilPressure = outletPressure; 
                else
                    obj.Nodes{ii}.oilPressure = nodesPressure(ii); 
                end                
                x_coor(ii,1) = obj.Nodes{ii}.x_coordinate;
            end
            
            %assign pressure values to links, since the surface where
            %flowrate is calculated through might pass through the links
            for ii = 1:obj.numberOfLinks
                if obj.Links{ii}.isInlet
                    obj.Links{ii}.oilPressure =...
                        (1+obj.Nodes{obj.Links{ii}.pore2Index}.oilPressure)/2;
                elseif obj.Links{ii}.isOutlet
                    obj.Links{ii}.oilPressure =...
                        obj.Nodes{obj.Links{ii}.pore1Index}.oilPressure/2;                    
                else
                    obj.Links{ii}.oilPressure =...
                        (obj.Nodes{obj.Links{ii}.pore1Index}.oilPressure + ...
                        obj.Nodes{obj.Links{ii}.pore2Index}.oilPressure) / 2;
                end
            end                     
           
            % Plot Pressure
            % a , b are 2 surfaces perpendicular to x-direction with
            % distance equals to intervalx
            x_outlet = max(x_coor);
            x_inlet = min(x_coor);
            n = 100;
            intervalx = (x_outlet - x_inlet)/n;
            a = x_inlet;
            b = x_inlet + intervalx;
            x = zeros(n,1);
            press_x = zeros(100,1);
            for i = 1:n
                area = 0;
                for ii = 1:obj.numberOfNodes
                    if obj.Nodes{ii}.x_coordinate >= a && obj.Nodes{ii}.x_coordinate < b
%                        press_x(i) = press_x(i) + obj.Nodes{ii}.waterPressure*obj.Nodes{ii}.area;
%                        area= area+obj.Nodes{ii}.area;                       
                       press_x(i) = press_x(i) + obj.Nodes{ii}.oilPressure;
                       area= area+1;                       
                    end
                end
                for ii = 1:obj.numberOfLinks    
                    if ~obj.Links{ii}.isOutlet
                        if obj.Nodes{obj.Links{ii}.pore2Index}.x_coordinate >= a && obj.Nodes{obj.Links{ii}.pore2Index}.x_coordinate < b
%                             press_x(i) = press_x(i) + obj.Links{ii}.waterPressure*obj.Links{ii}.area;
%                             area= area+obj.Links{ii}.area;                            
                            press_x(i) = press_x(i) + obj.Links{ii}.oilPressure;
                            area = area + 1;
                        end
                    end
                end
                press_x(i)=press_x(i)/area;
                x(i) = x(i) + i*intervalx;                
                a = a + intervalx;
                b = b + intervalx;
            end                  
%             plot(x, press_x, '*')
%             title('Oil_Pressure drop in x-direction')
%             xlabel('X(m)')
%             xlim([x_inlet x_outlet])
%             ylabel('Pressure(Pa)') 
        end 
       
        %% Flow rate calculation for each phase in the netwrok
        function calculateFlowRate(obj, inletPressure, outletPressure)
            % Fluid = water
            pressureDistribution_singlePhaseFlow(obj, inletPressure,outletPressure); 
            pressureDistribution_singlePhaseFlow_o (obj, inletPressure, outletPressure)
%             pressureDistribution_Cylindrical(obj, inletPressure,outletPressure); 
            obj.totalFlowRate = 0;
            obj.totalFlowRate_o = 0;
            cond_factor = obj.waterViscosity / obj.oilViscosity;
            % calculate flow rate in Inlet_Links
            for ii = 1:obj.numberOfLinks 
                
%                     node2Index = obj.Links{ii}.pore2Index;
                    node1Index = obj.Links{ii}.pore1Index;
                if obj.Links{ii}.isInlet 
                    
%                     %calculate the conductivity of the linkNode system
%                     nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
%                         obj.Links{ii}.conductance) +...
%                         ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance)))^-1;
%                     nodeLinkSystemConductance_o = ((obj.Links{ii}.linkLength /...
%                         obj.Links{ii}.conductance*cond_factor) +...
%                         ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance)))^-1;
%                     % calculate the flow rate of the fluid
%                     obj.totalFlowRate = obj.totalFlowRate + ...
%                         abs(nodeLinkSystemConductance * ...
%                         (inletPressure - obj.Nodes{node2Index}.waterPressure)); 
%                     obj.totalFlowRate_o = obj.totalFlowRate_o + ...
%                         abs(nodeLinkSystemConductance_o * ...
%                         (inletPressure - obj.Nodes{node2Index}.oilPressure));
                elseif obj.Links{ii}.isOutlet
                    %calculate the conductivity of the linkNode system
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance)))^-1;
                    nodeLinkSystemConductance_o = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance*cond_factor) +...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance)))^-1;
                    % calculate the flow rate of the fluid
                    obj.totalFlowRate = obj.totalFlowRate + ...
                        abs(nodeLinkSystemConductance * ...
                        (outletPressure - obj.Nodes{node1Index}.waterPressure)); 
                    obj.totalFlowRate_o = obj.totalFlowRate_o + ...
                        abs(nodeLinkSystemConductance_o * ...
                        (outletPressure - obj.Nodes{node1Index}.oilPressure));
                end
            end
            
            % calculate velocity through the network 
            obj.velocity = obj.totalFlowRate/(obj.yDimension * obj.zDimension); 
            
            % for quasi-static, capillaryNumber must be less than 10e-4
            obj.capillaryNumber = obj.waterViscosity * obj.velocity/ obj.sig_ow;  
        end
        
        %% AbsolutePermeability
        function calculateAbsolutePermeability(obj, inletPressure, outletPressure)
            %AbsolutePermeability calculates the absolute permeability of
            %the network
            calculateFlowRate(obj, inletPressure, outletPressure);
            
            % for pressure difference in the formula the corresponding
            % pressure drop between the vertical surfaces should be
            % calculated (based on Piri B1 formula)
            
            unitConvertor = 1/0.987*10^15; % unit conversion from m2 to miliDarcy
            obj.absolutePermeability = unitConvertor * obj.velocity * obj.xDimension * obj.waterViscosity/ (inletPressure -outletPressure );
        end  
        
        %% Pressure distribution calculation in Cylindrical pore_Single-Phase       
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
%             nodesPressure = pcg(Factor, B, 1e-7, 1000);
                        %assign the pressure values to each node
             x_coor = zeros(obj.numberOfNodes,1); 
            for ii = 1:obj.numberOfNodes
                if nodesPressure(ii) > inletPressure
                    obj.Nodes{ii}.waterPressure = inletPressure; 
                elseif nodesPressure(ii) < outletPressure
                    obj.Nodes{ii}.waterPressure = outletPressure; 
                else
                    obj.Nodes{ii}.waterPressure = nodesPressure(ii); 
                end                
                x_coor(ii,1) = obj.Nodes{ii}.x_coordinate;
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
           
            % Plot Pressure
            % a , b are 2 surfaces perpendicular to x-direction with
            % distance equals to intervalx
            x_outlet = max(x_coor);
            x_inlet = min(x_coor);
            n = 100;
            intervalx = (x_outlet - x_inlet)/n;
            a = x_inlet;
            b = x_inlet + intervalx;
            x = zeros(n,1);
            press_x = zeros(100,1);
            for i = 1:n
                area = 0;
                for ii = 1:obj.numberOfNodes
                    if obj.Nodes{ii}.x_coordinate >= a && obj.Nodes{ii}.x_coordinate < b
%                        press_x(i) = press_x(i) + obj.Nodes{ii}.waterPressure*obj.Nodes{ii}.area;
%                        area= area+obj.Nodes{ii}.area;                       
                       press_x(i) = press_x(i) + obj.Nodes{ii}.waterPressure;
                       area= area+1;                       
                    end
                end
                for ii = 1:obj.numberOfLinks    
                    if ~obj.Links{ii}.isOutlet
                        if obj.Nodes{obj.Links{ii}.pore2Index}.x_coordinate >= a && obj.Nodes{obj.Links{ii}.pore2Index}.x_coordinate < b
%                             press_x(i) = press_x(i) + obj.Links{ii}.waterPressure*obj.Links{ii}.area;
%                             area= area+obj.Links{ii}.area;                            
                            press_x(i) = press_x(i) + obj.Links{ii}.waterPressure;
                            area = area + 1;
                        end
                    end
                end
                press_x(i)=press_x(i)/area;
                x(i) = x(i) + i*intervalx;                
                a = a + intervalx;
                b = b + intervalx;
            end                  
%             plot(x, press_x, '*')
%             title('Pressure drop in x-direction')
%             xlabel('X(m)')
%             xlim([x_inlet x_outlet])
%             ylabel('Pressure(Pa)') 
        end   
        
        %% Conductance Calculation 2-Phase Flow and Water Saturation
        function Sw = calculateConductance_and_Saturation(obj, Pc, Cycle)  
                 
            waterVolume = 0;   
            vol = 0;
            waterArea = zeros(obj.numberOfLinks,4);
            
            for i = 1:obj.numberOfNodes
                if strcmp(Cycle , 'Drainage')== 1
                    [obj.Nodes{i}.waterCrossSectionArea, obj.Nodes{i}.waterConductance,...
                        obj.Nodes{i}.oilCrossSectionArea, obj.Nodes{i}.oilConductance] =...
                        obj.Nodes{i}.calculateConductance_Drainage(Pc, obj.Pc_drain_max);     
               
                else % imbibition 
                    [obj.Nodes{i}.waterCrossSectionArea, obj.Nodes{i}.waterConductance,...
                        obj.Nodes{i}.oilCrossSectionArea, obj.Nodes{i}.oilConductance] =...
                        obj.Nodes{i}.calculateConductance_Imbibition(obj, Pc);                       
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
                if strcmp(Cycle , 'Drainage')== 1                 
                   [obj.Links{i}.waterCrossSectionArea, obj.Links{i}.waterConductance,...
                        obj.Links{i}.oilCrossSectionArea, obj.Links{i}.oilConductance] =...
                        obj.Links{i}.calculateConductance_Drainage(Pc, obj.Pc_drain_max);     
                else % imbibition
                    [obj.Links{i}.waterCrossSectionArea, obj.Links{i}.waterConductance,...
                        obj.Links{i}.oilCrossSectionArea, obj.Links{i}.oilConductance] =...
                        obj.Links{i}.calculateConductance_Imbibition(obj, Pc);      
                    
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
            Sw = waterVolume / vol;     
        end
          
        %% Pressure distribution calculation in pore_Two-Phases 
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
            nodesWaterPressure = gmres(Factor_W, B_W,[], 1e-7, 1000);            
            nodesOilPressure = gmres(Factor_O, B_O,[], 1e-7, 1000);
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
          
        %% Relative Permeability_Drain
        function [krw, kro] = calculateRelativePermeability(obj, inletPressure, LinkL, cluster_A_nums)
              
            waterFlowRate = 0;   
            oilFlowRate = 0;
            
            %search through all the links
            for ii = 1:obj.numberOfLinks 
                  
                node2Index = obj.Links{ii}.pore2Index;
                if obj.Links{ii}.isInlet   
                     
                        % calculate the flow rate of the fluid
                        if any(LinkL(ii) == cluster_A_nums(:))  
                        waterFlowRate = waterFlowRate + ...
                            abs(obj.Links{ii}.nodeLinkSystemConductance_W * ...
                            (inletPressure - obj.Nodes{node2Index}.waterPressure));  
                        end
                        % calculate the flow rate of the fluid
                        oilFlowRate = oilFlowRate + ...
                            abs(obj.Links{ii}.nodeLinkSystemConductance_O * ...
                            (inletPressure - obj.Nodes{node2Index}.oilPressure));  
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
        
        %% Relative Permeability_imb
        function [krw, kro] = calculateRelativePermeability_Imb (obj, outletPressure, LinkL, LinkL_W, cluster_A_nums, cluster_A)
              
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
        
        %% Primary Drainage  
        function PrimaryDrainage(obj, inletPressure, outletPressure)              
             %% determining the capillary pressure level interval
             Pc_threshold = zeros(2*obj.numberOfLinks,1);  
             Pc_threshold_n = zeros(obj.numberOfLinks,1); 
             for i = 1:obj.numberOfLinks                
                 obj.Links{i}.thresholdPressure = obj.Links{i}.calculateThresholdPressurePistonLike();
                 Pc_threshold(i) = obj.Links{i}.thresholdPressure;
             end
             for i = 1:obj.numberOfNodes                
                  obj.Nodes{i}.thresholdPressure = obj.Nodes{i}.calculateThresholdPressurePistonLike();
                  Pc_threshold(i+obj.numberOfLinks) = obj.Nodes{i}.thresholdPressure;
             end 
             
             % Pc_interval
             max_Pc = max(Pc_threshold); 
             min_Pc = min(nonzeros(Pc_threshold)); 
             Pc_interval = (max_Pc-min_Pc)/10;   
             Pc = min_Pc; 
             lastMaxP = Pc;
             t = 2; 
             invaded = 0;
             obj.Pc_drain_max = max_Pc;
             obj.DrainageData = zeros(10,5);
             obj.DrainageData(1,:) = [1, 0, 1, 0, 0]; 
             
             %% Cycle of increasing Pressure
             while Pc <= obj.Pc_drain_max    
             [~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = Clustering_water(obj);    
              press = 1;
             deltaV = 0;
             %% Find new inlet-Links with threshold pressure < Pc             
             for i = 1:obj.numberOfLinks                  
                  node1Index = obj.Links{i}.pore1Index;
                  node2Index = obj.Links{i}.pore2Index;
                  if obj.Links{i}.isInlet && obj.Links{i}.occupancy == 'A'
                         if (any(LinkL(i) == cluster_A_nums(:)) ||  any(LinkL(i) == cluster_B_nums(:)))
                         Pc_threshold_n(i,1)= Pc_threshold(i);
                         deltaV = deltaV + obj.Links{i}.volume; 
                         end
                  elseif obj.Links{i}.isOutlet && obj.Links{i}.occupancy == 'A'                     
                           if obj.Nodes{node1Index}.occupancy == 'B'  
                               Pc_threshold_n(i,1)= Pc_threshold(i);
                               deltaV = deltaV + obj.Links{i}.volume ; 
                           end
                  elseif ~obj.Links{i}.isOutlet && ~obj.Links{i}.isInlet && obj.Links{i}.occupancy == 'A'                    
                      if obj.Nodes{node1Index}.occupancy == 'B' || obj.Nodes{node2Index}.occupancy == 'B'
                          if (any(LinkL(i) == cluster_A_nums(:)) ||  any(LinkL(i) == cluster_B_nums(:)))
                              Pc_threshold_n(i,1)= Pc_threshold(i); 
                              deltaV = deltaV + obj.Links{i}.volume;  
                          end
                      end  
                  end
%                 deltaS = deltaV /obj.poreVolume ;
%                  if deltaS > 0.05                       
%                      break;
%                  end
             end
             deltaS = 0;
             deltaV = 0;  
             if (any (nonzeros(Pc_threshold_n)) && min(nonzeros(Pc_threshold_n))<= Pc) || Pc == max_Pc
                 pressure = 1;
             else
                 pressure = 0;
             end 
             %% Add Links which have Pc_threshold < Pc in each steps and also have oil-saturated neighbour Node 
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
                 
                 if obj.Links{i}.isInlet && obj.Links{i}.occupancy == 'A' &&...
                                   (any(LinkL(i) == cluster_A_nums(:)) ||  any(LinkL(i) == cluster_B_nums(:)))
                         obj.Links{i}.occupancy = 'B'; 
                         invaded = invaded + 1;
                         deltaV = deltaV + obj.Links{i}.volume ;  
                         if  obj.Nodes{node2Index}.occupancy == 'A' &&obj.Nodes{node2Index}.thresholdPressure <=Pc  &&...
                                   (any(NodeL(node2Index) == cluster_A_nums(:)) ||  any(NodeL(node2Index) == cluster_B_nums(:)))
                             
                             obj.Nodes{node2Index}.occupancy = 'B'; 
                             invaded = invaded + 1;
                             deltaV = deltaV + obj.Nodes{node2Index}.volume ;  
                             for j=1:obj.Nodes{node2Index}.connectionNumber
                                 if obj.Nodes{node2Index}.connectedLinks(j)~=i
                                     Pc_threshold_n(obj.Nodes{node2Index}.connectedLinks(j),1)= Pc_threshold(obj.Nodes{node2Index}.connectedLinks(j));
                                 end
                             end
                         end
                  elseif obj.Links{i}.isOutlet && obj.Links{i}.occupancy == 'A'                     
                           if obj.Nodes{node1Index}.occupancy == 'B'  
                               
                               obj.Links{i}.occupancy = 'B'; 
                               invaded = invaded + 1;
                               deltaV = deltaV + obj.Links{i}.volume ;  
                           end
                  elseif ~obj.Links{i}.isOutlet && ~obj.Links{i}.isInlet && obj.Links{i}.occupancy == 'A'                        
                      if obj.Nodes{node1Index}.occupancy == 'B' || obj.Nodes{node2Index}.occupancy == 'B'  &&...
                              (any(LinkL(i) == cluster_A_nums(:)) ||  any(LinkL(i) == cluster_B_nums(:)))
                          
                          obj.Links{i}.occupancy = 'B';
                          invaded = invaded + 1;
                          deltaV = deltaV + obj.Links{i}.volume ;  
                          
                          if  obj.Nodes{node2Index}.occupancy == 'A' && obj.Nodes{node2Index}.thresholdPressure <=Pc  && ...                                  
                                   (any(NodeL(node2Index) == cluster_A_nums(:)) ||  any(NodeL(node2Index) == cluster_B_nums(:)))
                               
                              obj.Nodes{node2Index}.occupancy = 'B'; 
                              invaded = invaded + 1;
                              deltaV = deltaV + obj.Nodes{node2Index}.volume ;  
                              for j=1:obj.Nodes{node2Index}.connectionNumber
                                 if obj.Nodes{node2Index}.connectedLinks(j)~=i
                                     Pc_threshold_n(obj.Nodes{node2Index}.connectedLinks(j),1)= Pc_threshold(obj.Nodes{node2Index}.connectedLinks(j));
                                 end
                             end
                          end   
                          
                          if obj.Nodes{node1Index}.occupancy == 'A' && obj.Nodes{node1Index}.thresholdPressure <=Pc && ...
                                   (any(NodeL(node1Index) == cluster_A_nums(:)) ||  any(NodeL(node1Index) == cluster_B_nums(:)))
                               
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
             end
               if Pc == max_Pc
                   lastMaxP = max_Pc;
               end
             % Updating element saturations and conductances
             Sw = calculateConductance_and_Saturation(obj, lastMaxP, 'Drainage'); 
             pressureDistribution_TwoPhases(obj, inletPressure, outletPressure); 
             
             % Relative Permeability Calculation 
             [Krw , Kro] = calculateRelativePermeability (obj, inletPressure, LinkL, cluster_A_nums);            
             obj.DrainageData(t,:) = [Sw, lastMaxP, Krw, Kro, invaded]; 
             
              % Pc Step Calculation 
              if press ~= 0 
                 Pc = Pc + Pc_interval;   
              end
             t = t + 1;     
             end 
              %% Ploting Pc & Kr
            S = obj.DrainageData(:,1);
            P = obj.DrainageData(:,2);
            kw = obj.DrainageData(:,3);
            ko = obj.DrainageData(:,4);   
            figure('name','Primary Dranage Cappilary Pressure & Relative Permeability Curves',...
                'units','normalized','outerposition',[0 0 1 1])
            subplot(1,2,1);
            grid on
            plot(S,P,'--r')  
            title('Drainage & Imbibition Cappilary Pressure Curves')
            xlabel('Sw')
            xlim([0 1])
            ylabel('Pc (Pa)')

            subplot(1,2,2);
            plot(S,kw,'--r',S,ko,'--b')
            title('Drainage Relative Permeability Curves')
            xlabel('Sw')
            xlim([0 1])
            ylabel('Reative Permeability')
            ylim([0 1])
            legend('Water Relative Permeability','Oil Relative Permeability','Location','North') 

        end   
        
        %% Clustering
        function [NumberOfClusters, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = Clustering(obj)
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
            
            % Modify number of inlet & outlet Links of Clusters
            for i =1:obj.numberOfLinks
                if obj.Links{i}.isInlet 
                    if any(obj.Links{i}.oilLayerExist) || obj.Links{i}.occupancy == 'B'
                        if obj.Nodes{obj.Links{i}.pore2Index}.occupancy == 'B' || any(obj.Nodes{obj.Links{i}.pore2Index}.oilLayerExist)
                            LinkL(i,1) = NodeL(obj.Links{i}.pore2Index);
                        end
                    end
                    
                elseif obj.Links{i}.isOutlet
                    if (obj.Links{i}.occupancy == 'B' || any(obj.Links{i}.oilLayerExist))
                        if obj.Nodes{obj.Links{i}.pore1Index}.occupancy == 'B' || any(obj.Nodes{obj.Links{i}.pore1Index}.oilLayerExist)
                            LinkL(i,1) = NodeL(obj.Links{i}.pore1Index);                
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
%             numA = 0;
            for i = 1:obj.numberOfNodes
                
                NodeS(i,1) = obj.Nodes{i}.occupancy; % nodes with oil(B),nodes with water(A)
                if obj.Nodes{i}.shapeFactor <= 1 / 16 % There is water in the corners
                        NodeS(i,1) = 'A';                              
%                         numA=1;
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
                     if obj.Links{i}.shapeFactor <= 1 / 16 % There is water in the corners
                        LinkS(i,1) = 'A';  
                     end
                end
            end            
            OFlag = 'A'; %oil clusters has numbers B:NumberOfClusters %water clusters are A
            
%             if numA == 1
            % HoshenKopelman Algorithm for clustering
            [NumberOfClusters, NodeL, LinkL] = modifiedHKNonLattice(NodeS, LinkS,NodeNext, LinksOfNode, OFlag);
            
            % Modify number of inlet & outlet Links of Clusters
            for i =1:obj.numberOfLinks
                if obj.Links{i}.isInlet 
                    if  obj.Links{i}.occupancy == 'A' || obj.Links{i}.shapeFactor <= 1 / 16
                        if obj.Nodes{obj.Links{i}.pore2Index}.occupancy == 'A' || obj.Nodes{obj.Links{i}.pore2Index}.shapeFactor <= 1 / 16
                            LinkL(i,1) = NodeL(obj.Links{i}.pore2Index);
                        end
                    end
                    
                elseif obj.Links{i}.isOutlet
                    if obj.Links{i}.occupancy == 'A' || obj.Links{i}.shapeFactor <= 1 / 16
                        if obj.Nodes{obj.Links{i}.pore1Index}.occupancy == 'A' || obj.Nodes{obj.Links{i}.pore1Index}.shapeFactor <= 1 / 16
                            LinkL(i,1) = NodeL(obj.Links{i}.pore1Index);                
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
%             else
%                  NumberOfClusters =0;     
%                  NodeL=[];
%                  LinkL=[];
%                  cluster_A_nums=[];
%                  cluster_B_nums=[];
%             end
        end 
       
        %% Secondary Imbibition 
        function ScoendaryImbibition_Spon(obj, inletPressure, outletPressure)    
            
            %counter for invaded elements
            numOfLinks_SnapOff = 0;
            numOfLinks_PistoneLike = 0;
            numOfLinks_LayerCollapse = 0;
            numOfNodes_SnapOff = 0;
            numOfNodes_PoreBodyFilling = 0;
            numOfNodes_LayerCollapse = 0;
            thresholdPressure = zeros(obj.numberOfLinks, 14);               
             Pc_imb = obj.Pc_drain_max; 
             
            %% Calculating throat Snap-Off & Pistone-Like displacement & layer collapse            
            for i = 1:obj.numberOfLinks 
                obj.Links{i}.recedingContactAngle = 5*pi/18; 
                if obj.Links{i}.occupancy == 'B' % if the throat is oil filled                           
                    [obj.Links{i}.ThresholdPressure_PistonLike, obj.Links{i}.oilLayerExist] =...
                        obj.Links{i}.calculateThresholdPressurePistonLike_Imbibition (obj.Pc_drain_max);
                    obj.Links{i}.ThresholdPressure_SnapOff = ...
                        obj.Links{i}.calculateThresholdPressureSnapOff (obj.Pc_drain_max); 
                    obj.Links{i}.ThresholdPressure_LayerCollapse = obj.Links{i}.calculateThresholdPressureLayerCollapse (obj.Pc_drain_max);
                    if obj.Links{i}.isInlet
                        thresholdPressure(i,1) = -1;
                    elseif obj.Links{i}.isOutlet
                        thresholdPressure(i,1) = 1;
                    end
                    thresholdPressure(i,2) = obj.Links{i}.ThresholdPressure_SnapOff; 
                    thresholdPressure(i,3) = obj.Links{i}.ThresholdPressure_PistonLike;
                    thresholdPressure(i,4:7) = obj.Links{i}.ThresholdPressure_LayerCollapse;
                end
            end
            
            %% Calculating Pore Snap-Off & Pore-Body Filling displacement & layer collapse            
            for i = 1:obj.numberOfNodes                 
                obj.Nodes{i}.recedingContactAngle = 5*pi/18;
                if obj.Nodes{i}.occupancy == 'B' % if the throat is oil filled
                    
                    obj.Nodes{i}.ThresholdPressure_SnapOff = ...
                        obj.Nodes{i}.calculateThresholdPressureSnapOff (obj.Pc_drain_max); 
                    [obj.Nodes{i}.ThresholdPressure_PoreBodyFilling, obj.Nodes{i}.oilLayerExist] = ...
                        obj.Nodes{i}.calculateThresholdPressurePoreBodyFilling (obj);                    
                    obj.Nodes{i}.ThresholdPressure_LayerCollapse = obj.Nodes{i}.calculateThresholdPressureLayerCollapse (obj.Pc_drain_max);
                    if obj.Nodes{i}.isInlet
                        thresholdPressure(i,8) = -1;
                    elseif obj.Nodes{i}.isOutlet
                        thresholdPressure(i,8) = 1;
                    end
                    thresholdPressure(i,9) = obj.Nodes{i}.ThresholdPressure_SnapOff;
                    thresholdPressure(i,10) = obj.Nodes{i}.ThresholdPressure_PoreBodyFilling;                    
                    thresholdPressure(i,11:14) = obj.Nodes{i}.ThresholdPressure_LayerCollapse;
                end  
            end  
             
             poreVolumeInjected = 0;
             PVInjected = 0;
             Pc_interval = Pc_imb /1000; 
             t = 1;       
             obj.ImbibitionData = zeros(100,11);
             obj.ImbibitionData(1,:) = ...
                 [obj.DrainageData(end,1),obj.Pc_drain_max,obj.DrainageData(end,3),obj.DrainageData(end,4),0,0,0,0,0,0,0];
             
             percList = -100000*ones(2*obj.numberOfNodes+2*obj.numberOfLinks,4);
            
            [~,~, ~, cluster_A_nums, cluster_B_nums] = Clustering(obj);  
             
            while (~isempty(cluster_A_nums) || ~isempty(cluster_B_nums)) && Pc_imb > -99999

             [~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = Clustering(obj);             
             [~, NodeL_W, LinkL_W,cluster_A,cluster_B] = Clustering_water(obj);
             %% Percolation List           
                a = 1; 
                for i = 1:obj.numberOfLinks         
                    
                    node1Index = obj.Links{i}.pore1Index;
                    node2Index = obj.Links{i}.pore2Index;
                    
                    if (any(LinkL(i) == cluster_A_nums(:)) || any(LinkL(i) == cluster_B_nums(:))) 
                        
                         if obj.Links{i}.isInlet            
                             
                             if any(obj.Links{i}.ThresholdPressure_PistonLike)&& obj.Links{i}.ThresholdPressure_PistonLike<= Pc_imb
                                
                                percList(a,1) = i;
                                percList(a,2) = 2;
                                percList(a,3) = 2;
                                percList(a,4) = obj.Links{i}.ThresholdPressure_PistonLike ;
                                a = a + 1;
                             end
                            
                         elseif obj.Links{i}.isOutlet
                             
                             if (any(NodeL_W(node1Index) == cluster_A(:)) || any(NodeL_W(node1Index) == cluster_B(:))) && ...
                                     any(obj.Links{i}.ThresholdPressure_PistonLike)&& ...
                                     obj.Links{i}.ThresholdPressure_PistonLike<= Pc_imb && obj.Nodes{node1Index}.occupancy == 'A'
                                 
                                     percList(a,1) = i;
                                     percList(a,2) = 2;
                                     percList(a,3) = 2;
                                     percList(a,4) = obj.Links{i}.ThresholdPressure_PistonLike ;
                                     a = a + 1;
                             elseif  any(obj.Links{i}.ThresholdPressure_SnapOff) ... % if the throat is non circular
                                     && obj.Links{i}.ThresholdPressure_SnapOff<= Pc_imb 
                                     percList(a,1) = i;
                                     percList(a,2) = 2;
                                     percList(a,3) = 1;
                                     percList(a,4) = obj.Links{i}.ThresholdPressure_SnapOff;
                                     a = a + 1;
                             end        
                            
                         else
                             if  ((any(NodeL_W(node1Index) == cluster_A(:)) || any(NodeL_W(node1Index) == cluster_B(:))) ||...
                                 (any(NodeL_W(node2Index) == cluster_A(:)) || any(NodeL_W(node2Index) == cluster_B(:))) ) && ...
                                 (obj.Nodes{node1Index}.occupancy == 'A' || obj.Nodes{node2Index}.occupancy == 'A') && ...
                                     any(obj.Links{i}.ThresholdPressure_PistonLike) && obj.Links{i}.ThresholdPressure_PistonLike<= Pc_imb
                                     
                                     percList(a,1) = i;
                                     percList(a,2) = 2;
                                     percList(a,3) = 2;
                                     percList(a,4) = obj.Links{i}.ThresholdPressure_PistonLike ;
                                     a = a + 1;
                             elseif any(obj.Links{i}.ThresholdPressure_SnapOff)...% if the throat is non circular
                                     && obj.Links{i}.ThresholdPressure_SnapOff<= Pc_imb
                                 percList(a,1) = i;
                                 percList(a,2) = 2;
                                 percList(a,3) = 1;
                                 percList(a,4) = obj.Links{i}.ThresholdPressure_SnapOff;
                                 a = a + 1;
                             end
                         end
                    end
                end
                for i = 1:obj.numberOfNodes
                    
                    if (any(NodeL(i) == cluster_A_nums(:)) || any(NodeL(i) == cluster_B_nums(:)))  
                        
                        filledThroats = 0;
                        T = 0;
                        for j = 1:obj.Nodes{i}.connectionNumber                            
                            if (obj.Links{obj.Nodes{i}.connectedLinks(j)}.occupancy == 'A')                                 
                                filledThroats = filledThroats + 1;
                                if (any(LinkL_W(obj.Nodes{i}.connectedLinks(j))) == cluster_A(:)) ||...
                                        (any(LinkL_W(obj.Nodes{i}.connectedLinks(j))) == cluster_B(:))
                                    T = 1;
                                end
                            end
                        end                                         
                        
                        if filledThroats ~= 0 
                            
                            [obj.Nodes{i}.ThresholdPressure_PoreBodyFilling, obj.Nodes{i}.oilLayerExist] = ...
                                obj.Nodes{i}.calculateThresholdPressurePoreBodyFilling (obj); 
                            if any(obj.Nodes{i}.ThresholdPressure_PoreBodyFilling) && ...
                                obj.Nodes{i}.ThresholdPressure_PoreBodyFilling <= Pc_imb && T == 1
                         
                                percList(a,1) = i;
                                percList(a,2) = 1;
                                percList(a,3) = 2;
                                percList(a,4) = obj.Nodes{i}.ThresholdPressure_PoreBodyFilling; % pore body filling threshold pressure
                                a = a + 1;
                            end
                        elseif any(obj.Nodes{i}.ThresholdPressure_SnapOff)...% if the node is non circular
                                && obj.Nodes{i}.ThresholdPressure_PoreBodyFilling <= Pc_imb 
                            
                            percList(a,1) = i;
                            percList(a,2) = 1;
                            percList(a,3) = 1;
                            percList(a,4) = obj.Nodes{i}.ThresholdPressure_SnapOff; % snap off threshold pressure
                            a = a + 1;
                        end
                        
                    end
                end
                
                percList = sortrows(percList,-4);
                if Pc_imb > percList(1,4) && percList(1,4) ~= -100000
                    Pc_imb = percList(1,4);
                end
                t = t+1; 
                new =a;
                %% Percolation Section                
                while (~isempty(cluster_A_nums) || ~isempty(cluster_B_nums))  && Pc_imb >-99999 && new > 0
                    
                    % Descending sorting of threshold pressures
                    percList = sortrows(percList,-4);                      
                    new = new - 1;         
                     
                        %% if the first element is a throat
                        if percList(1,2) == 2 
                            
                            linkIndex = percList(1,1);
                            node1Index = obj.Links{linkIndex}.pore1Index;
                            node2Index = obj.Links{linkIndex}.pore2Index;
                             
                            if any(LinkL(linkIndex) == cluster_A_nums(:)) || any(LinkL(linkIndex) == cluster_B_nums(:))
                                
                                if percList(1,3) == 2                                                 
                                    
                                    if obj.Links{linkIndex}.isInlet 
                                          
                                        obj.Links{linkIndex}.occupancy = 'A';  
                                        
                                        poreVolumeInjected = poreVolumeInjected + obj.Links{linkIndex}.volume; 
                                        numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;  
                                        obj.Links{linkIndex}.isInvaded = true;
                                        
                                         if  obj.Nodes{node2Index}.occupancy == 'B'
                                            
                                            % Updating pore body filling of the pore
                                            [obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling, obj.Nodes{node2Index}.oilLayerExist] = ...
                                                obj.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                            if any(obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling) && ...
                                                    Pc_imb <= obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling   
                                                new = new+1;
                                                percList(end-new,:) = [node2Index, 1, 2, obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling];
                                            end
                                        end
                                        
                                    elseif obj.Links{linkIndex}.isOutlet
                                        
                                        if  obj.Nodes{node1Index}.occupancy == 'A' 
                                            
                                            obj.Links{linkIndex}.occupancy = 'A';   
                                            
                                            poreVolumeInjected = poreVolumeInjected + obj.Links{linkIndex}.volume; 
                                            numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;   
                                            obj.Links{linkIndex}.isInvaded = true;
                                            
                                            if obj.Nodes{node1Index}.occupancy == 'B'
                                                
                                                % Updating pore body filling of the pore
                                                [obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling, obj.Nodes{node1Index}.oilLayerExist] = ...
                                                    obj.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                                                                                
                                                if any(obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling) && ...
                                                    Pc_imb <= obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling
                                                
                                                    new = new+1;
                                                    percList(end-new,:) = [node1Index, 1, 2, obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling];  
                                                end
                                            end
                                        end
                                        
                                    elseif ~obj.Links{linkIndex}.isOutlet && ~obj.Links{linkIndex}.isInlet
                                        
                                        if ((obj.Nodes{node1Index}.occupancy == 'A'  ) ||...
                                                (obj.Nodes{node2Index}.occupancy == 'A' ))
                                             
                                            obj.Links{linkIndex}.occupancy = 'A'; 
                                                
                                            poreVolumeInjected = poreVolumeInjected + obj.Links{linkIndex}.volume;
                                            numOfLinks_PistoneLike = numOfLinks_PistoneLike + 1;
                                            obj.Links{linkIndex}.isInvaded = true; 
                                                                                         
                                            if obj.Nodes{node1Index}.occupancy == 'B'
                                                
                                                % Updating pore body filling of the pore
                                                [obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling, obj.Nodes{node1Index}.oilLayerExist] = ...
                                                    obj.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                                
                                                if any(obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling) && ...
                                                    Pc_imb <= obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling
                                                    
                                                    new = new+1;
                                                    percList(end-new,:) = [node1Index, 1, 2, obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling];                                    
                                                end
                                            end
                                            
                                            if  obj.Nodes{node2Index}.occupancy == 'B'
                                                
                                                % Updating pore body filling of the pore
                                                [obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling, obj.Nodes{node2Index}.oilLayerExist] = ...
                                                    obj.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                                
                                                if any(obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling) && ...
                                                    Pc_imb <= obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling 
                                                    
                                                    new = new+1;
                                                    percList(end-new,:) = [node2Index, 1, 2, obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling];
                                                end
                                            end
                                        end
                                    end
                                    
                                elseif percList(1,3)==1
                                    
                                    obj.Links{linkIndex}.occupancy = 'A'; 
                                    
                                    poreVolumeInjected = poreVolumeInjected + obj.Links{linkIndex}.volume; 
                                    numOfLinks_SnapOff = numOfLinks_SnapOff + 1;  
                                    obj.Links{linkIndex}.isInvaded = true; 
                                                                         
                                    if obj.Links{linkIndex}.isOutlet  
                                         
                                        if  obj.Nodes{node1Index}.occupancy == 'B'
                                            
                                            % Updating pore body filling of the pore
                                            [obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling, obj.Nodes{node1Index}.oilLayerExist] = ...
                                                obj.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                            
                                            if any(obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling) && ...
                                                    Pc_imb <= obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling                                                
                                                
                                                new = new+1;
                                                percList(end-new,:) = [node1Index, 1, 2, obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling];  
                                            end
                                        end
                                        
                                    elseif ~obj.Links{linkIndex}.isOutlet && ~obj.Links{linkIndex}.isInlet
                                        
                                        if obj.Nodes{node1Index}.occupancy == 'B'
                                                                                    
                                            % Updating pore body filling of the pore
                                            [obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling, obj.Nodes{node1Index}.oilLayerExist] = ...
                                                obj.Nodes{node1Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                                                                    
                                            if any(obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling) && ...
                                                    Pc_imb <= obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling
                                                
                                                new = new+1;                                          
                                                percList(end-new,:) = [node1Index, 1, 2, obj.Nodes{node1Index}.ThresholdPressure_PoreBodyFilling];                                    
                                            end
                                        end
                                        
                                        if obj.Nodes{node2Index}.occupancy == 'B'
                                            
                                            % Updating pore body filling of the pore
                                            [obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling, obj.Nodes{node2Index}.oilLayerExist] = ...
                                                obj.Nodes{node2Index}.calculateThresholdPressurePoreBodyFilling (obj);
                                            
                                            if any(obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling) && ...
                                                    Pc_imb <= obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling
                                                
                                                new = new+1;
                                                percList(end-new,:) = [node2Index, 1, 2, obj.Nodes{node2Index}.ThresholdPressure_PoreBodyFilling];
                                            end
                                        end
                                    end
                                end
                            end
                            
                        %% if the first element is a pore
                        elseif percList(1,2) == 1 
                            
                            nodeIndex = percList(1,1);
                            
                            if any(NodeL(nodeIndex) == cluster_A_nums(:)) || any(NodeL(nodeIndex) == cluster_B_nums(:))
                                
                                if percList(1,3) == 2  
                                      
                                    obj.Nodes{nodeIndex}.occupancy = 'A'; % make the pore water type 
                                    
                                    poreVolumeInjected = poreVolumeInjected + obj.Nodes{nodeIndex}.volume;
                                    numOfNodes_PoreBodyFilling = numOfNodes_PoreBodyFilling + 1;   
                                    obj.Nodes{nodeIndex}.isInvaded = true; 
                                    
                                    for ii = 1:obj.Nodes{nodeIndex}.connectionNumber 
                                        
                                        connectedLink = obj.Nodes{nodeIndex}.connectedLinks(ii);
                                        
                                        if obj.Links{connectedLink}.occupancy == 'B'&& ...
                                                (any(obj.Links{connectedLink}.ThresholdPressure_PistonLike)) && ...
                                                obj.Links{connectedLink}.ThresholdPressure_PistonLike >= Pc_imb
                                                
                                            new = new +1;
                                            percList(end-new,:) = [connectedLink, 2 ,2 , obj.Links{connectedLink}.ThresholdPressure_PistonLike];
                                        end
                                    end
                                elseif percList(1,3) == 1
                                    
                                    obj.Nodes{nodeIndex}.occupancy = 'A';
                                    
                                    obj.Nodes{nodeIndex}.oilLayerExist = nan;
                                    poreVolumeInjected = poreVolumeInjected + obj.Nodes{nodeIndex}.volume;                                                   
                                    numOfNodes_SnapOff = numOfNodes_SnapOff + 1;  
                                    obj.Nodes{nodeIndex}.isInvaded = true; 
                                                                                                                                            
                                    for ii = 1:obj.Nodes{nodeIndex}.connectionNumber                                  
                                        
                                        connectedLink = obj.Nodes{nodeIndex}.connectedLinks(ii);
                                        
                                        if obj.Links{connectedLink}.occupancy == 'B' && ...
                                                (any(obj.Links{connectedLink}.ThresholdPressure_PistonLike)) && ...
                                                obj.Links{connectedLink}.ThresholdPressure_PistonLike >= Pc_imb
                                            
                                            new = new +1;
                                            percList(end-new,:) = [connectedLink, 2 ,2 , obj.Links{connectedLink}.ThresholdPressure_PistonLike];
                                        end
                                    end
                                end
                            end
                        end
                        
                        percList(1,:) = -100000;    
                        [~, NodeL, LinkL,cluster_A_nums,cluster_B_nums] = Clustering(obj);  
                end
                                
                invaded = numOfLinks_SnapOff + numOfLinks_PistoneLike + ...
                    numOfNodes_SnapOff + numOfNodes_PoreBodyFilling + numOfNodes_LayerCollapse;
                
                %% Updating saturations and conductances 
                Sw_imb = calculateConductance_and_Saturation(obj, Pc_imb, 'Imbibition');    
                pressureDistribution_TwoPhases(obj, inletPressure, outletPressure); 
                
                PVInjected = poreVolumeInjected + PVInjected ; 
                [Krw_imb, Kro_imb] =...
                    calculateRelativePermeability_Imb(obj, outletPressure, LinkL, LinkL_W, cluster_A_nums, cluster_A);
                obj.ImbibitionData(t,:) = ...
                    [Sw_imb,Pc_imb,Krw_imb, Kro_imb,invaded, ...
                    numOfLinks_SnapOff,numOfLinks_PistoneLike, ...
                    numOfLinks_LayerCollapse,numOfNodes_SnapOff, ...
                    numOfNodes_PoreBodyFilling,numOfNodes_LayerCollapse];
                
                Pc_imb = Pc_imb - Pc_interval;                         
                [~, ~, ~,cluster_A_nums,cluster_B_nums] = Clustering(obj); 
            end
             
             %% Ploting Pc & Kr
            S = obj.DrainageData(:,1);
            P = obj.DrainageData(:,2);
            kw = obj.DrainageData(:,3);
            ko = obj.DrainageData(:,4); 
            SI = obj.ImbibitionData(1:t,1);
            PI = obj.ImbibitionData(1:t,2);
            kwi = obj.ImbibitionData(1:t,3);
            koi = obj.ImbibitionData(1:t,4);  
            figure('name','Primary Dranage & Secondary Imbibition Cappilary Pressure & Relative Permeability Curves',...
                'units','normalized','outerposition',[0 0 1 1])
            subplot(2,2,[1 3]);
            grid on
            plot(S,P,'--r') 
            hold on
            plot(SI,PI,'--b')
            title('Drainage & Imbibition Cappilary Pressure Curves')
            xlabel('Sw')
            xlim([0 1])
            ylabel('Pc (Pa)')

            subplot(2,2,2);
            plot(S,kw,'--r',S,ko,'--b')
            title('Drainage Relative Permeability Curves')
            xlabel('Sw')
            xlim([0 1])
            ylabel('Reative Permeability')
            ylim([0 1])
            legend('Water Relative Permeability','Oil Relative Permeability','Location','North')
            subplot(2,2,4);
            plot(SI,kwi,'--r',SI,koi,'--b')
            title('Imbibition Relative Permeability Curves')
            xlabel('Sw')
            xlim([0 1])
            ylabel('Reative Permeability')
            ylim([0 1])
            legend('Water Relative Permeability','Oil Relative Permeability','Location','North')


        end   
        %% vtk file generation
        function vtkOutput(obj)
            vtkFileID = fopen('output.vtk','w');
            if vtkFileID == -1
                error('Cannot open file for writing.');
            end
            title = 'output';
            fprintf ( vtkFileID, '# vtk DataFile Version 2.0\n' );
            fprintf ( vtkFileID, '%s\n', title );
            fprintf ( vtkFileID, 'ASCII\n' );
            fprintf ( vtkFileID, '\n' );
            fprintf ( vtkFileID, 'DATASET POLYDATA\n' );
            fprintf ( vtkFileID, 'POINTS %d double\n', obj.numberOfNodes );
            for i = 1:obj.numberOfNodes
                fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate );
            end
            
        end
    end
end

