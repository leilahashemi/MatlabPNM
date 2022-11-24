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
        waterSaturation
        
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
        
        %% Network Properties calculation
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
                    obj.numOfCircularElements = obj.numOfCircularElements+1;
                elseif strcmp(obj.Nodes{ii}.geometry , 'Triangle')== 1
                    obj.numOfTriangularElements = obj.numOfTriangularElements+1;
                else
                    obj.numOfSquareElements = obj.numOfSquareElements+1;
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
        
        %% Pressure distribution calculation of single phase flow       
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
            nodesPressure = gmres(Factor, B,[], 1e-10, 1000);
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
               
        %% Flow rate calculation for each phase in the netwrok
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
        
        %% AbsolutePermeability
        function calculateAbsolutePermeability(obj, inletPressure, outletPressure) 
            calculateFlowRate(obj, inletPressure, outletPressure);                        
            unitConvertor = 1/0.987*10^15; % unit conversion from m2 to miliDarcy
            obj.absolutePermeability = unitConvertor * obj.velocity * obj.xDimension * obj.waterViscosity/ (inletPressure -outletPressure );
            format longE
            obj.absolutePermeability_m2 = obj.velocity * obj.xDimension * obj.waterViscosity/ (inletPressure -outletPressure );
        end  
        
    end
end

