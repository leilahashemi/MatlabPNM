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
            nodesPressure = gmres(Factor, B,[], 1e-7, 1000);
            
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
        
        %% Flow rate calculation for each phase in the netwrok
        function calculateFlowRate(obj, inletPressure, outletPressure)
            % Fluid = water
%             pressureDistribution_singlePhaseFlow(obj, inletPressure,outletPressure); 
            pressureDistribution_Cylindrical(obj, inletPressure,outletPressure); 
            obj.totalFlowRate = 0;
            
            % calculate flow rate in Inlet_Links
            for ii = 1:obj.numberOfLinks 
                
                    node2Index = obj.Links{ii}.pore2Index;
                    
                if obj.Links{ii}.isInlet 
                    
                    %calculate the conductivity of the linkNode system
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance) +...
                        ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance)))^-1;
                    
                    % calculate the flow rate of the fluid
                    obj.totalFlowRate = obj.totalFlowRate + ...
                        abs(nodeLinkSystemConductance * ...
                        (inletPressure - obj.Nodes{node2Index}.waterPressure));  
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
            
            unitConvertor = 1.01325E+15; % unit conversion from m2 to miliDarcy
            obj.absolutePermeability = unitConvertor * obj.velocity * obj.xDimension * obj.waterViscosity;
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
        
        %% Single phase Reactive transport_Raoof_2010 
        function calculateReactiveTransport_SinglePhase_circle(obj, inletPressure, outletPressure, soluteConcentration, poreVolumeInjected)
            % Flow & ReactiveTransport Based on Raoof paper 2010
            
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
            timeStep = timeStep * 10;             
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
            PV_time = poreVolumeInjected/obj.totalFlowRate;  
            timePlot = zeros(round(PV_time/timeStep)+1,1);
            obj.BreakThroughCurve_singlePhase = zeros(round(PV_time/timeStep)+1,2);
            flux_averagedConcentration = zeros(round(PV_time/timeStep)+1,1);
            
            while t<50 %time < PV_time  
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
            
            nodesConcentration_new = gmres(Factor, Known,[], 1e-10, 1000); 
            
            % asign new concentration of nodes
            for i = 1:obj.numberOfNodes
                if nodesConcentration_new(i) > soluteConcentration
                obj.Nodes{i}.concentration(t+1) = soluteConcentration;
                elseif nodesConcentration_new(i) < 0
                 obj.Nodes{i}.concentration(t+1) = 0;
                else
                obj.Nodes{i}.concentration(t+1) = nodesConcentration_new(i);
                end 
                % calculation for BreakThroughCurve at outlet of network
%                     sumOfConcentration = sumOfConcentration + ...
%                         obj.Nodes{i}.concentration(t)*flowRate_node(i);
%                     sumOfFlowRate = sumOfFlowRate + flowRate_node(i);
            end
            % calculate BreakThroughCurve at outlet of network
%             flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration;
           
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
            flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration;
            end 
            
            obj.BreakThroughCurve_singlePhase(:,1) = timePlot;
            obj.BreakThroughCurve_singlePhase(:,2) = flux_averagedConcentration;
            plot(timePlot, flux_averagedConcentration,'*'); 
            title('Break Through Curve')
            xlabel('Time(s)')            
            ylabel('DimensionlessConcentration(-)')
%             ylim([0 1]
        end
        function calculateReactiveTransport_SinglePhase(obj, inletPressure, outletPressure, soluteConcentration, poreVolumeInjected)
            % Flow & ReactiveTransport Based on Raoof paper 2010
            
            % mass transfer coefficient: alpha & distribution coefficient: K_d
            alpha = 45/24/3600*ones(obj.numberOfLinks,1);
            K_d = 45/45*ones(obj.numberOfLinks,1);
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
                    flowRate_link(ii) = obj.Links{ii}.conductance * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);  
                elseif obj.Links{ii}.isInlet 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.conductance * ...
                        abs(inletPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);   
                else 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.conductance * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        outletPressure);                     
                    obj.totalFlowRate = obj.totalFlowRate + flowRate_link(ii);
                end  
                residenceTime_link(ii) = obj.Links{ii}.volume/flowRate_link(ii);
            end
            
            % Set time step to avoid numerical dispersion
            timeStep = min(nonzeros(residenceTime_link));
            timeStep = timeStep*500 ;             
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
%             PV_time = poreVolumeInjected/obj.totalFlowRate;  
            timePlot = zeros(50,1);
            obj.BreakThroughCurve_singlePhase = zeros(50,2);
            flux_averagedConcentration = zeros(50,1);
%             
%             timePlot = zeros(round(PV_time/timeStep)+1,1);
%             obj.BreakThroughCurve_singlePhase = zeros(round(PV_time/timeStep)+1,2);
%             flux_averagedConcentration = zeros(round(PV_time/timeStep)+1,1)
            
            while t<50 %time < PV_time  
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
            
            nodesConcentration_new = gmres(Factor, Known,[], 1e-10, 1000); 
            
            % asign new concentration of nodes
            for i = 1:obj.numberOfNodes
                if nodesConcentration_new(i) > soluteConcentration
                obj.Nodes{i}.concentration(t+1) = soluteConcentration;
                elseif nodesConcentration_new(i) < 0
                 obj.Nodes{i}.concentration(t+1) = 0;
                else
                obj.Nodes{i}.concentration(t+1) = nodesConcentration_new(i);
                end 
                % calculation for BreakThroughCurve at outlet of network
%                     sumOfConcentration = sumOfConcentration + ...
%                         obj.Nodes{i}.concentration(t)*flowRate_node(i);
%                     sumOfFlowRate = sumOfFlowRate + flowRate_node(i);
            end
            % calculate BreakThroughCurve at outlet of network
%             flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration;
           
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
            flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration;
            end 
            
            obj.BreakThroughCurve_singlePhase(:,1) = timePlot;
            obj.BreakThroughCurve_singlePhase(:,2) = flux_averagedConcentration;
            plot(timePlot, flux_averagedConcentration,'*'); 
            title('Break Through Curve')
            xlabel('Time(s)')            
            ylabel('DimensionlessConcentration(-)')
%             ylim([0 1]
        end
        
        %% Single phase Reactive transport_Edgar
        function calculateReactiveTransport_Edgar(obj, inletPressure, outletPressure, soluteConcentration, poreVolumeInjected)
            % Flow & ReactiveTransport Based on Edgar 2019 
            
            % calculate pressure distribution
            pressureDistribution_Cylindrical (obj, inletPressure, outletPressure); 
            
            residenceTime_link = zeros(obj.numberOfLinks,1);
            flowRate_link = zeros(obj.numberOfLinks,1);
            obj.totalFlowRate = 0;            
            flowRate_node_In = zeros(obj.numberOfNodes,1);            
                        
            % calculate flowrate of links residence time based on eq. 1 & 10                         
            for ii = 1:obj.numberOfLinks        
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;
                
                if ~obj.Links{ii}.isInlet && ~obj.Links{ii}.isOutlet      
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure); 
                    
                    if obj.Nodes{node1Index}.waterPressure >= obj.Nodes{node2Index}.waterPressure
                        flowRate_node_In (node2Index) = flowRate_node_In (node2Index) + flowRate_link (ii);
                    else
                        flowRate_node_In (node1Index) = flowRate_node_In (node1Index) + flowRate_link (ii);
                    end
                        
                elseif obj.Links{ii}.isInlet 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(inletPressure - ...                         
                        obj.Nodes{node2Index}.waterPressure);  
                    obj.totalFlowRate = obj.totalFlowRate + flowRate_link(ii);                    
                    flowRate_node_In (node2Index) = flowRate_node_In (node2Index) + flowRate_link (ii);
                else 
                    
                    % calculate the flow rate of the fluid
                    flowRate_link(ii) = obj.Links{ii}.cylindricalConductance * ...
                        abs(obj.Nodes{node1Index}.waterPressure - ...                         
                        outletPressure); 
                end 
                residenceTime_link(ii) = obj.Links{ii}.volume/flowRate_link(ii);
            end
            timeStep = min(nonzeros(residenceTime_link));  
            timeStep =  timeStep / 10;
            
            
            flowRate_node = zeros(obj.numberOfNodes,1);
            diff_node = zeros(obj.numberOfNodes,1);
            dif = 10^(-9);
            Nt = obj.numberOfLinks + obj.numberOfNodes; 
            Factor = zeros(Nt, Nt);
            Known = zeros(Nt, 1);
            
            % calculation of 3 Unknowns (concentration of nodes & 
            % concentration & adsorbedConcentration of links) in each timeStep 
            
            t = 0; 
            time = 0;
            timePlot = zeros(round(poreVolumeInjected/timeStep)+1 ,1);
            flux_averagedConcentration = zeros(round(poreVolumeInjected/timeStep)+1 ,1);
            
            while t<50  %time < poreVolumeInjected 
                t = t+1;
                time = time + timeStep; 
                timePlot(t) = time;
                sumOfConcentration = 0;
                sumOfFlowRate = 0;
            % calculate concentration of nodes: based on eq.7
            for i = 1:Nt 
                if i <= obj.numberOfNodes 
                for j = 1:obj.Nodes{i}.connectionNumber 
                    
                    connectedLinkIndex = obj.Nodes{i}.connectedLinks(j);                     
                    connectedNodeIndex = obj.Nodes{i}.connectedNodes(j); 
                    l = obj.numberOfNodes + connectedLinkIndex;
                                        
                    if connectedNodeIndex ~= 0 && connectedNodeIndex ~= -1
                        flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex); 
                        diff_node(i) = diff_node(i) + dif * obj.Links{connectedLinkIndex}.area/obj.Links{connectedLinkIndex}.linkLength;
                        Factor(i, l) = -1 * timeStep / obj.Nodes{i}.volume * ...
                            (dif * obj.Links{connectedLinkIndex}.area/obj.Links{connectedLinkIndex}.linkLength);
                        if obj.Nodes{connectedNodeIndex}.waterPressure > obj.Nodes{i}.waterPressure
                            Factor(i, l) = Factor(i, l)+ ...
                                -1 * timeStep / obj.Nodes{i}.volume *flowRate_link(connectedLinkIndex);
                        end
                    elseif connectedNodeIndex == -1
                        flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex); 
                    diff_node(i) = diff_node(i) + dif * obj.Links{connectedLinkIndex}.area/obj.Links{connectedLinkIndex}.linkLength;
                    Factor(i, l) = -1 * timeStep / obj.Nodes{i}.volume * ...
                       (dif * obj.Links{connectedLinkIndex}.area/obj.Links{connectedLinkIndex}.linkLength+flowRate_link(connectedLinkIndex));
                    end
                end
                
                Known(i,1) = obj.Nodes{i}.concentration(t);
                Factor(i, i) = 1 - timeStep * (-1 * flowRate_node_In(i)-diff_node(i)) / obj.Nodes{i}.volume  ;  
                
                else      
                    j = i-obj.numberOfNodes;
                    Known(i,1) = obj.Links{j}.concentration(t); 
                    if ~obj.Links{j}.isInlet && ~obj.Links{j}.isOutlet
                        Factor(i, obj.Links{j}.pore1Index) = -1 * timeStep / obj.Links{j}.volume * ...
                            (flowRate_link(j) + dif * obj.Links{j}.area/obj.Links{j}.linkLength);
                        Factor(i, obj.Links{j}.pore2Index) = -1 * timeStep / obj.Links{j}.volume * ...
                            (flowRate_link(j) + dif * obj.Links{j}.area/obj.Links{j}.linkLength);
                    elseif obj.Links{j}.isInlet
                        Factor(i, obj.Links{j}.pore2Index) = -1 * timeStep / obj.Links{j}.volume * ...
                            (flowRate_link(j) + dif * obj.Links{j}.area/obj.Links{j}.linkLength);
                         Known(i,1) =  Known(i,1) + soluteConcentration *(timeStep / obj.Links{j}.volume * ...
                            (flowRate_link(j) + dif * obj.Links{j}.area/obj.Links{j}.linkLength));
                    else
                        Factor(i, obj.Links{j}.pore1Index) = 1 + timeStep / obj.Links{j}.volume * flowRate_link(j);
                        Known(i,1) =  Known(i,1) +  obj.Nodes{obj.Links{j}.pore1Index}.concentration(t);
                    end
                    Factor(i, i) = 1 + timeStep / obj.Links{j}.volume *2*( flowRate_link(j)-...
                        dif * obj.Links{j}.area/obj.Links{j}.linkLength); 
                end 
            end
             
            nodesConcentration_new = gmres (Factor, Known,[], 1e-10, 2000);
%              nodesConcentration_new = Factor\ Known;
            % asign new concentration of nodes
            for i = 1:obj.numberOfNodes 
                if nodesConcentration_new(i) > soluteConcentration
                obj.Nodes{i}.concentration(t+1) = soluteConcentration;
                elseif nodesConcentration_new(i)< 0
                  obj.Nodes{i}.concentration(t+1) = 0;
                else
                obj.Nodes{i}.concentration(t+1) = nodesConcentration_new(i);
                end 
            end
           
            % calculate new concentration & adsorbedConcentration of links:
            % based on eq.3&4
             for j = 1:obj.numberOfLinks
                 i = j + obj.numberOfNodes;
                 
                if nodesConcentration_new(i) > soluteConcentration
                obj.Links{j}.concentration(t+1) = soluteConcentration;
                elseif nodesConcentration_new(i) < 0
                obj.Links{j}.concentration(t+1) = 0;
                else
                obj.Links{j}.concentration(t+1) = nodesConcentration_new(i);
                end 
                
                if obj.Links{j}.isOutlet
                    % calculation for BreakThroughCurve at outlet of network
                    sumOfConcentration = sumOfConcentration + ...
                        obj.Nodes{node1Index}.concentration(t)*flowRate_node(node1Index);
                    sumOfFlowRate = sumOfFlowRate + flowRate_node(node1Index);
                end
             end              
            % calculate BreakThroughCurve at outlet of network
            flux_averagedConcentration(t,1) = sumOfConcentration / sumOfFlowRate / soluteConcentration;
            end
            
            obj.BreakThroughCurve_singlePhase(:,1) = timePlot;
            obj.BreakThroughCurve_singlePhase(:,2) = flux_averagedConcentration;
            plot(timePlot,flux_averagedConcentration,'*'); 
            title('Break Through Curve')
            xlabel('Time(s)')            
            ylabel('DimensionlessConcentration(-)')
%             ylim([0 1])
        end       
       
        %% Single phase Reactive transport_Raoof_2017
        function calculateReactiveTransport_SinglePhaseDiffusion(obj, inletPressure, outletPressure, soluteConcentration, poreVolumeInjected)
            % Flow & ReactiveTransport Based on Raoof paper 2017 
            
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
            timeStep = timeStep * 0.01;    
            
            flowRate_node = zeros(obj.numberOfNodes,1); 
            diffusion_node = zeros(obj.numberOfNodes,1); 
            N_total = obj.numberOfNodes+obj.numberOfLinks;
            Factor = zeros(N_total, N_total);
            Known = zeros(N_total, 1);
            
            % calculation of 3 Unknowns (concentration of nodes & 
            % concentration & adsorbedConcentration of links) in each timeStep 
            
             
            t = 0; 
            time = 0;
            timePlot = zeros(round(poreVolumeInjected/timeStep)+1 ,1);
            flux_averagedConcentration = zeros(round(poreVolumeInjected/timeStep)+1 ,1);
            
            while t<50 %time < poreVolumeInjected
                t = t+1;
                time = time + timeStep; 
                timePlot(t) = time;
                sumOfConcentration = 0;
                sumOfFlowRate = 0;
                
            % calculate concentration of nodes & links based on eq. 8 & 9
            for i = 1:N_total 
                
                if i <= obj.numberOfNodes % Element is Node
                    
                for j = 1:obj.Nodes{i}.connectionNumber 
                    
                    connectedLinkIndex = obj.Nodes{i}.connectedLinks(j);
                    connectedNodeIndex = obj.Nodes{i}.connectedNodes(j);
%                     adsorbedConcentration = obj.Links{connectedLinkIndex}.adsorbedConcentration(t);
%                     linksConcentration = obj.Links{connectedLinkIndex}.concentration(t);
                    
                    
                    diffusion_node(i) = diffusion_node(i) + diffusion_link(connectedLinkIndex);  
                    Factor(i, i+connectedLinkIndex) = -1 * timeStep / obj.Nodes{i}.volume *  diffusion_link(connectedLinkIndex);  
                    
                    if connectedNodeIndex ~= 0 && connectedNodeIndex ~= -1
                        
                        % determine link flowing into this node
                        if obj.Nodes{connectedNodeIndex}.waterPressure > obj.Nodes{i}.waterPressure 
                            
                            flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex);     
                            Factor(i, i+connectedLinkIndex) = Factor(i, i+connectedLinkIndex) +...
                                -1 * timeStep / obj.Nodes{i}.volume *  flowRate_link(connectedLinkIndex);   
                        end
                            
                    elseif connectedNodeIndex == -1
                                                
                            flowRate_node(i) = flowRate_node(i) + flowRate_link(connectedLinkIndex);     
                            Factor(i, i+connectedLinkIndex) = Factor(i, i+connectedLinkIndex) + ...
                                -1 * timeStep / obj.Nodes{i}.volume *  flowRate_link(connectedLinkIndex);  
                    end
                end
                Factor(i, i) = 1+ timeStep * (flowRate_node(i)+  diffusion_node(i)) / obj.Nodes{i}.volume  ;
                Known(i,1) = obj.Nodes{i}.concentration(t);
                
                else % Element is Link
                    j = i - obj.numberOfNodes;
                    node1Index = obj.Links{j}.pore1Index;
                    node2Index = obj.Links{j}.pore2Index;
                    
                    Factor(i, i) = 1 + timeStep * (flowRate_link(j)+  diffusion_link(j)) / obj.Links{j}.volume;
                    
                    if ~obj.Links{j}.isInlet && ~obj.Links{j}.isOutlet
                        
                        Known(i,1) = obj.Links{j}.concentration(t);
                        Factor(i, node1Index) = -1 * timeStep / obj.Links{j}.volume * diffusion_link(j);
                        Factor(i, node2Index) = -1 * timeStep / obj.Links{j}.volume * diffusion_link(j);
                        if obj.Nodes{node1Index}.waterPressure > obj.Nodes{node2Index}.waterPressure 
                            Factor(i, node1Index) = Factor(i, node1Index) - timeStep / obj.Links{j}.volume * flowRate_link(j);
                        else
                            Factor(i, node2Index) = Factor(i, node2Index) - timeStep / obj.Links{j}.volume * flowRate_link(j);
                        end
                        
                    elseif obj.Links{j}.isInlet
                        Known(i,1) = obj.Links{j}.concentration(t) + ...
                            soluteConcentration * timeStep / obj.Links{j}.volume * (flowRate_link(j)+ diffusion_link(j));
                        Factor(i, node2Index) = -1 * timeStep / obj.Links{j}.volume * diffusion_link(j);
                    else
                        Known(i,1) = obj.Links{j}.concentration(t);
                        Factor(i, node1Index) = -1 * timeStep / obj.Links{j}.volume * (flowRate_link(j) + diffusion_link(j));
                    end
                end
                
            end
             nodesConcentration_new = gmres (Factor, Known,[], 1e-10, 1000);
            
            % asign new concentration of nodes
            for i = 1:obj.numberOfNodes
                if nodesConcentration_new(i) > soluteConcentration
                obj.Nodes{i}.concentration(t+1) = soluteConcentration;
                else
                obj.Nodes{i}.concentration(t+1) = nodesConcentration_new(i);
                end
            end
            
            % calculate new concentration & adsorbedConcentration of links:
            % based on eq.3&4
            for i = 1:obj.numberOfLinks 
                j = obj.numberOfNodes+i;
                if nodesConcentration_new(j) > soluteConcentration
                obj.Links{i}.concentration(t+1) = soluteConcentration;
                else
                obj.Links{i}.concentration(t+1) = nodesConcentration_new(obj.numberOfNodes+i);
                end
                node1Index = obj.Links{i}.pore1Index;  
                
                if obj.Links{i}.isOutlet
                    % calculation for BreakThroughCurve at outlet of network
                    sumOfConcentration = sumOfConcentration + ...
                        obj.Nodes{node1Index}.concentration(t)*flowRate_node(node1Index);
                    sumOfFlowRate = sumOfFlowRate + flowRate_node(node1Index);
                end 
            end
            % calculate BreakThroughCurve at outlet of network
            flux_averagedConcentration(t) = sumOfConcentration / sumOfFlowRate / soluteConcentration;
            end 
            
            obj.BreakThroughCurve_singlePhase(:,1) = timePlot;
            obj.BreakThroughCurve_singlePhase(:,2) = flux_averagedConcentration;
            plot(timePlot,flux_averagedConcentration,'*');
            title('Break Through Curve')
            xlabel('Time(s)')            
            ylabel('DimensionlessConcentration(-)')
%             ylim([0 1])
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

