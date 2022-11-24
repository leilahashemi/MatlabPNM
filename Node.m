classdef Node < Element
    
    properties
        x_coordinate
        y_coordinate
        z_coordinate
        connectionNumber
        connectedNodes
        connectedLinks 
        imbThresholdPressure_PoreBodyFilling = nan;
    end
    
    methods
        function obj = Node(index,...
                            x_coordinate,...
                            y_coordinate,...
                            z_coordinate,...
                            connectionNumber,...
                            connectionData,...
                            volume,...
                            radius,...
                            shapeFactor,...
                            clayVolume)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.index = index;
            obj.x_coordinate = x_coordinate;
            obj.y_coordinate = y_coordinate;
            obj.z_coordinate = z_coordinate;
            obj.connectionNumber = connectionNumber;  
            obj.connectedNodes = connectionData(1:connectionNumber);
            if connectionData(connectionNumber + 1) == 1
                obj.isInlet = true;
            else
                obj.isInlet = false;
            end
            if connectionData(connectionNumber + 2) == 1
                obj.isOutlet = true;
            else
                obj.isOutlet = false;
            end
            obj.connectedLinks = connectionData(connectionNumber + 3: end);
            obj.volume = volume;
            obj.radius = radius;
            obj.shapeFactor = shapeFactor;
            obj.clayVolume = clayVolume;   
        end   
         %% PoreBodyFilling 
        %Patzek
        function calculateThresholdPressurePoreBodyFilling_P (obj,network) 
         % Based on Patzek: eqs 42-49
         obj.oilLayerExist(1,:) = nan;
         W = [0;0.72;0.45;1.2;1.5;5];
         attachedThroats = obj.connectedLinks;
         
         z = 0;% number of oil filled attached throats
         for i = 1:obj.connectionNumber
             if network.Links{attachedThroats(i)}.occupancy == 'B' 
                 z=z+1; 
             end
         end    
         
         if z == 0
             obj.imbThresholdPressure_PoreBodyFilling = nan;
         elseif z == 1
%              obj.calculateThresholdPressurePistonLike_Imbibition(network.Pc_drain_max);
            obj.imbThresholdPressure_PoreBodyFilling = obj.imbThresholdPressure_PistonLike;
         else            
            if z > 5
                w = W(6);
            else
                w = W(z);
            end             
             nominator = 0;
             denominator = 0;
             sumOfThroatRadius = 0;
             for ii = 1:z
                 randNumber = rand;
                 sumOfThroatRadius = sumOfThroatRadius + network.Links{oilFilledAttachedThroats(ii)}.radius;
                 nominator = nominator + randNumber * sumOfThroatRadius;
                 denominator = denominator + randNumber;
             end
             R_ave = (obj.radius + w * nominator / denominator)/cos(obj.advancingContactAngle);
             obj.imbThresholdPressure_PoreBodyFilling = 2*obj.sig_ow/R_ave; 
         end        
        end
        %Blunt1
        function calculateThresholdPressurePoreBodyFilling_Blunt1 (obj,network)
         % Based on Blunt1
         obj.oilLayerExist(1,:) = nan;
         W = [0;2.5;5;20;100];
         attachedThroats = obj.connectedLinks;
         
         z = 0;% number of oil filled attached throats
         for i = 1:obj.connectionNumber
             if network.Links{attachedThroats(i)}.occupancy == 'B' 
                 z=z+1; 
             end
         end    
         
         if z == 0
             obj.imbThresholdPressure_PoreBodyFilling = nan;
         elseif z == 1
%              obj.calculateThresholdPressurePistonLike_Imbibition(network.Pc_drain_max);
            obj.imbThresholdPressure_PoreBodyFilling = obj.imbThresholdPressure_PistonLike;
         else            
             if z > 5
                 w = W(5);
             else
                 w = W(z);
             end
             nominator = 0;
             
             for ii = 1:z
                  
                 randNumber = rand;
                 nominator = nominator + randNumber * w;
             end
             obj.imbThresholdPressure_PoreBodyFilling = 2*obj.sig_ow * cos(obj.advancingContactAngle)/(obj.radius + nominator); 
         end        
        end
        %Blunt2
        function calculateThresholdPressurePoreBodyFilling (obj,network)
         % Based on Blunt2
         obj.oilLayerExist(1,:) = nan;
         W =15000;
         attachedThroats = obj.connectedLinks;
         
         z = 0;% number of oil filled attached throats
         for i = 1:obj.connectionNumber
             if network.Links{attachedThroats(i)}.occupancy == 'B' 
                 z=z+1; 
             end
         end    
         
         if z == 0
             obj.imbThresholdPressure_PoreBodyFilling = nan;
         elseif z == 1
%              obj.calculateThresholdPressurePistonLike_Imbibition(network.Pc_drain_max);
            obj.imbThresholdPressure_PoreBodyFilling = obj.imbThresholdPressure_PistonLike;
         else             
            nominator = 0;
             for ii = 1:z
                 randNumber = rand;
                 nominator = nominator + randNumber * W;
             end
             obj.imbThresholdPressure_PoreBodyFilling = 2*obj.sig_ow * cos(obj.advancingContactAngle)/obj.radius - obj.sig_ow *nominator; 
        end
        
        end
        %Valvatne
        function calculateThresholdPressurePoreBodyFilling_V(obj,network)
         % Based on Valvatne 
         obj.oilLayerExist(1,:) = nan;
         W = 0.03/sqrt(network.absolutePermeability/1.01325E+15);
         attachedThroats = obj.connectedLinks; 
         
         z = 0;% number of oil filled attached throats
         for i = 1:obj.connectionNumber
             if network.Links{attachedThroats(i)}.occupancy == 'B' 
                 z=z+1; 
             end
         end    
         
         if z == 0
             obj.imbThresholdPressure_PoreBodyFilling = nan;
         elseif z == 1
%              obj.calculateThresholdPressurePistonLike_Imbibition(network.Pc_drain_max);
            obj.imbThresholdPressure_PoreBodyFilling = obj.imbThresholdPressure_PistonLike;
         else            
            nominator = 0;
             for ii = 1:z
                 randNumber = rand;
                 nominator = nominator + randNumber * W;
             end
             obj.imbThresholdPressure_PoreBodyFilling = 2*obj.sig_ow * cos(obj.advancingContactAngle)/obj.radius - obj.sig_ow *nominator; 
         end 
        end
        %Oren1
        function calculateThresholdPressurePoreBodyFilling_Oren1 (obj,network)
         % Based on Oren1
         obj.oilLayerExist(1,:) = nan;
         W = [0;0.5;1;2;5;10];
         attachedThroats = obj.connectedLinks;
         
         z = 0;% number of oil filled attached throats
         for i = 1:obj.connectionNumber
             if network.Links{attachedThroats(i)}.occupancy == 'B' 
                 z=z+1; 
             end
         end    
         
         if z == 0
             obj.imbThresholdPressure_PoreBodyFilling = nan;
         elseif z == 1
%              obj.calculateThresholdPressurePistonLike_Imbibition(network.Pc_drain_max);
            obj.imbThresholdPressure_PoreBodyFilling = obj.imbThresholdPressure_PistonLike;
         else                        
            nominator = 0;
             for ii = 1:z
                 if z > 5
                     w = W(6);
                 else
                     w = W(z);
                 end 
                 randNumber = rand;
                 nominator = nominator + randNumber * w * network.Links{oilFilledAttachedThroats(ii,1)}.radius;
             end
             obj.imbThresholdPressure_PoreBodyFilling = 2*obj.sig_ow * cos(obj.advancingContactAngle)/(obj.radius + nominator); 
         end        
        end
        %Oren2
        function calculateThresholdPressurePoreBodyFilling_O2 (obj,network) 
         % Based on Oren2
         obj.oilLayerExist(1,:) = nan;
         W = [0;0.5;1;2;5;10];
         attachedThroats = obj.connectedLinks;
         
         z = 0;% number of oil filled attached throats
         for i = 1:obj.connectionNumber
             if network.Links{attachedThroats(i)}.occupancy == 'B' 
                 z=z+1; 
             end
         end    
         
         if z == 0
             obj.imbThresholdPressure_PoreBodyFilling = nan;
         elseif z == 1
%              obj.calculateThresholdPressurePistonLike_Imbibition(network.Pc_drain_max);
            obj.imbThresholdPressure_PoreBodyFilling = obj.imbThresholdPressure_PistonLike;
         else                        
            nominator = 0;
             for ii = 1:z
                 if z > 5
                     w = W(6);
                 else
                     w = W(z);
                 end 
                 randNumber = rand;
                 nominator = nominator + randNumber * w * network.Links{oilFilledAttachedThroats(ii,1)}.radius;
             end
             obj.imbThresholdPressure_PoreBodyFilling = (1 + 2 * sqrt(pi * obj.shapeFactor))*obj.sig_ow * cos(obj.advancingContactAngle)/...
                 (obj.radius + nominator); 
         end        
        end
        %Piri
        function calculateThresholdPressurePoreBodyFilling_Piri (obj,network)
         % Based on Piri
         obj.oilLayerExist(1,:) = nan;
         W =0.03 * 10^(-6);
         attachedThroats = obj.connectedLinks;
         
         z = 0;% number of oil filled attached throats
         for i = 1:obj.connectionNumber
             if network.Links{attachedThroats(i)}.occupancy == 'B' 
                 z=z+1; 
             end
         end    
         
         if z == 0
             obj.imbThresholdPressure_PoreBodyFilling = nan;
         elseif z == 1
%              obj.calculateThresholdPressurePistonLike_Imbibition(network.Pc_drain_max);
            obj.imbThresholdPressure_PoreBodyFilling = obj.imbThresholdPressure_PistonLike;
         else                        
            nominator = 0;
             for ii = 1:z
                 randNumber = rand;
                 nominator = nominator + randNumber * W;
             end
             obj.imbThresholdPressure_PoreBodyFilling = 2*obj.sig_ow * cos(obj.advancingContactAngle)/obj.radius - obj.sig_ow *nominator; 
        end
        
        end
    end
end

