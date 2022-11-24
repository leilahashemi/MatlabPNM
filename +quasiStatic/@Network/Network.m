classdef Network < handle & Fluids
    %Network Summary of this class goes here
    %   This class contain nodes and links of the network
    
    properties
        Nodes
        Links
        xDimension
        yDimension
        zDimension
        numberOfLinks
        numberOfNodes
        
        Porosity
        totalFlowRate
        absolutePermeability
        poreVolume
        Sw_drain
        Pc_drain_curve
        
    end
    
    methods
        %% Cunstructor function
        function obj = Network(fileName)
            %Network Construct an instance of this class
            %   Detailed explanation goes here
            
            % current path
            currentFoldet = pwd;
            networkFileFullPath = strcat(currentFoldet, '\Input\NetworkDataFile\' , fileName);

            
            % Opening the files
            link_1_fileID = fopen(strcat(networkFileFullPath, '_link1.dat'));
            obj.numberOfLinks = str2num(fgetl(link_1_fileID));
            link_2_fileID = fopen(strcat(networkFileFullPath, '_link2.dat'));
            
            node_2_fileID = fopen(strcat(networkFileFullPath, '_node2.dat'));            
            node_1_fileID = fopen(strcat(networkFileFullPath, '_node1.dat'));
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
        
    %%
    % Porosity calculation
    obj = calculatePorosity(obj)    
        % Pressure distribution calculation        
        pressureDistribution(obj, inletPressure, outletPressure)
        % Conductance Calculation
          calculateConductance(obj, Pc)      
        % Saturation and Conductance Calculation
        Sw_drain = calculateSaturations(obj, Pc)

        % This function calculates the flow rate for each phase in the netwrok
        calculateFlowRate(obj)
        % Calculate AbsolutePermeability 
        calculateAbsolutePermeability(obj)
        % Primary Drainage  
        PrimaryDrainage(obj)
        
        % vtk file generation
        vtkOutput(obj)
    end
end

