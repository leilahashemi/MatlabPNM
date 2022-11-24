classdef Element<Fluids
    
    properties
        % properties set from the link and node input files
        index
        radius
        shapeFactor
        volume
        clayVolume
        isInlet
        isOutlet
        
        % Calculated properties
        
        geometry % geometrical shape of the element
        surfaces
        points
        halfAngle1
        halfAngle2
        halfAngle3
        halfAngle4
        area   
        
        % Element conductance & area
        conductanceSinglePhase         
        waterConductance
        oilConductance        
        waterCrossSectionArea 
        oilCrossSectionArea  
        waterSaturation
        oilSaturation                          
        
        waterPressure
        oilPressure     
        
        occupancy = 'A';  % Element filled by Water
        b = zeros(1,4);
        L_AM = zeros(1,4);             
        oilLayerExist = nan(1,4);      
        waterCornerExist = nan(1,4);   
        isInvaded = false; % check for invasion in Imbibition
        
        % Capillary Threshold Pressure in Imbibition
        imbPressureTrapped = nan;
        imbThresholdPressure_SnapOff = nan;
        imbThresholdPressure_PistonLike = nan;
        imbThresholdPressure_LayerCollapse = nan(1,4);
        % Capillary Threshold Pressure in Drainage
        drainThresholdPressure_PistonLike 
                
        hingeAngles 
        newContactAngle
        recedingContactAngle = 0*pi/18;
        advancingContactAngle = 0*pi/18;
        
        control
        
        % reactive
        concentration = 0;
        adsorbedConcentration = 0;
    end    
    methods
        function calculateElementsProperties(obj) 
            % Geometry and conductance specification of the elements is
            % based of : Patzek, T. W., & Silin, D. B. (2001). Shape factor and hydraulic conductance in noncircular capillaries: I. One-phase creeping flow. Journal of Colloid and Interface Science. https://doi.org/10.1006/jcis.2000.7413
            % For ducts with square cross-sections, all four half-angles are equal to /4? and G = 1/16 . Circular ducts have no corners and G =1/ 4? . For simplicity, all ducts with shape factors between those of equilateral triangle and square can be mapped onto squares, and those with shape factors above 1/16 onto circles.
            % we'd better to insert star shapes later
            if  obj.shapeFactor <= sqrt(3) / 36
                obj.shapeFactor = max(obj.shapeFactor, 10^-7);
                obj.geometry = 'Triangle';
                betha2_min = atan((2 / sqrt(3)) * cos((acos(-12 * sqrt(3) * obj.shapeFactor)) / 3 + (4 * pi / 3)));
                betha2_max = atan((2 / sqrt(3)) * cos((acos(-12 * sqrt(3) * obj.shapeFactor)) / 3 ));
                % rand(0.25-0.75)
                rand_1 = 0.5*(rand+0.5);
                obj.halfAngle2 = betha2_min + rand_1 * (betha2_max - betha2_min);
                obj.halfAngle1 = -0.5 * obj.halfAngle2 + 0.5 * asin((tan(obj.halfAngle2) + 4 * obj.shapeFactor) * sin(obj.halfAngle2) / (tan(obj.halfAngle2) - 4 * obj.shapeFactor));
                obj.halfAngle3 = pi / 2 - obj.halfAngle1 - obj.halfAngle2;
                obj.halfAngle4 = nan;
                obj.area = obj.radius^2/4/obj.shapeFactor;                
                obj.conductanceSinglePhase = 3 * obj.area^2 * obj.shapeFactor /obj.waterViscosity / 5; 
                obj.surfaces = 5;
                obj.points = 6;
            elseif obj.shapeFactor > sqrt(3) / 36 && obj.shapeFactor <= 1 / 16
                obj.geometry = 'Square';
                obj.halfAngle1 = pi / 4;
                obj.halfAngle2 = pi / 4;
                obj.halfAngle3 = pi / 4;
                obj.halfAngle4 = pi / 4;
                obj.area = 4*obj.radius^2;                
                obj.conductanceSinglePhase = 0.5623 * obj.area^2 * obj.shapeFactor /obj.waterViscosity;
                obj.surfaces = 6;
                obj.points = 8;
            elseif obj.shapeFactor > 1 / 16
                obj.geometry = 'Circle';
                obj.halfAngle1 = nan;
                obj.halfAngle2 = nan;
                obj.halfAngle3 = nan;
                obj.halfAngle4 = nan;
                obj.area = pi*obj.radius^2;                
                obj.conductanceSinglePhase = 0.5 * obj.area^2 * obj.shapeFactor /obj.waterViscosity; 
                obj.surfaces = 8;
                obj.points = 6;
            end
                if obj.volume == 0
                    obj.volume = obj.area * obj.length;
                end
        end
         
        %% Drainage
        
        % Piston-Like
        % Based ob Al-Futaisi & Patzek 2001
        function calculateThresholdPressurePistonLike_drainage_Patzek (obj)
         % eqs 2-5
             if strcmp(obj.geometry , 'Circle')== 1
                 obj.drainThresholdPressure_PistonLike = 2*obj.sig_ow *cos(obj.recedingContactAngle)/obj.radius;
             else                 
                 nominator = 0;
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 for j = 1:4
                     if ~isnan(halfAngles(j)) 
                         E2 = cos(obj.recedingContactAngle + halfAngles(j)) *...
                             cos(obj.recedingContactAngle) / sin(halfAngles(j));
                         E0 = pi / 2 - obj.recedingContactAngle - halfAngles(j);
                         nominator = nominator +  (E2 - E0); 
                     end
                 end
                 obj.drainThresholdPressure_PistonLike = (obj.sig_ow / obj.radius)*...
                     cos(obj.recedingContactAngle)*(1+sqrt(1 -(4*obj.shapeFactor*...
                     nominator)/(cos(obj.recedingContactAngle)^2)));
             end
        end
        % Blunt 2017
        function calculateThresholdPressurePistonLike_drainage (obj)
             if strcmp(obj.geometry , 'Circle')== 1
                 obj.drainThresholdPressure_PistonLike = 2*obj.sig_ow *cos(obj.recedingContactAngle)/obj.radius;
             else % eq 3.10 - 3.12    
                 D = pi - 2 / 3 * obj.recedingContactAngle +3 * sin(obj.recedingContactAngle) * cos(obj.recedingContactAngle) - ...
                     cos(obj.recedingContactAngle)^2 / (4*obj.shapeFactor);
                 F_d = ( 1 + sqrt(1+4* obj.shapeFactor * D / cos(obj.recedingContactAngle)^2))/(1+2*sqrt(pi * obj.shapeFactor));
                 obj.drainThresholdPressure_PistonLike = (obj.sig_ow / obj.radius)*...
                     (1+2*sqrt(pi * obj.shapeFactor))*cos(obj.recedingContactAngle)* F_d; 
             end
        end
        
        % Conductance Calculation _ Drainage
        % Blunt 2017
        function calculateConductance_Drainage(obj, Pc)  
             
             if obj.occupancy == 'B' %element was occupied by oil
                 
                 if strcmp(obj.geometry , 'Circle')== 1                     
                 obj.waterCrossSectionArea = 0;    
                 obj.waterSaturation = 0;                 
                 obj.waterConductance = 0;
                 obj.oilCrossSectionArea = obj.area;
                 obj.oilSaturation = 1;
                 obj.oilConductance = obj.area^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                 
                 else
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 cornerArea = zeros(1,4);
                 cornerConductance = zeros(1,4);
                 %Raduis of Curvature
                 Rc = obj.sig_ow / abs(Pc);
                 for jj = 1:4
                     if ~isnan(halfAngles(jj)) && halfAngles(jj) < pi/2 - obj.recedingContactAngle
                         
                         obj.waterCornerExist(1,jj) = 1;
                         
                         % Apex & length of corner: eq 3.6-3.7
                         obj.b(jj) = obj.sig_ow / Pc * ...
                             cos(obj.recedingContactAngle + halfAngles(jj)) * sin(halfAngles(jj));                         
                         obj.L_AM(jj) = 2 * Rc * (pi / 2 - obj.recedingContactAngle - halfAngles(jj)); 
                         
                         % Area: eq 3.8
                         cornerArea(jj) = Rc ^2 * (cos(obj.recedingContactAngle)*...
                             cos(halfAngles(jj)+ obj.recedingContactAngle) / sin (halfAngles(jj))+ ...
                             obj.recedingContactAngle + halfAngles(jj) - pi/2);
                         
                         %Conductance
                         % Based on Piri_2005: eq B(10 - 15)
                         f = 1; % no-flow boundary condition suitable for oil-water interfaces
                         F1 = pi/2 - halfAngles(jj) - obj.recedingContactAngle;
                         F2 = cot(halfAngles(jj)) * cos(obj.recedingContactAngle) - sin(obj.recedingContactAngle); 
                         F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));
                         
                         if (obj.recedingContactAngle <= pi/2 - halfAngles(jj))  % Positive Curvature                       
                             cornerConductance(jj) = (cornerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                                 (F2 * cos(obj.recedingContactAngle) - F1) * F3 ^ 2) / ...
                                 (12 * obj.waterViscosity * ((sin(halfAngles(jj))) * ...
                                 (1 - F3) * (F2 + f * F1))^ 2);
                         elseif (obj.recedingContactAngle > pi/2 - halfAngles(jj)) % Negative Curvature
                             cornerConductance(jj) = (cornerArea(jj)^2 * tan(halfAngles(jj))* ...
                                 (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                                 (12 * obj.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
                         end
                     end
                     obj.waterCrossSectionArea = sum(cornerArea); 
                     obj.waterConductance = sum(cornerConductance);
                 end
                 obj.oilCrossSectionArea = obj.area - obj.waterCrossSectionArea; 
                 if strcmp(obj.geometry , 'Triangle')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                 end                 
                 obj.waterSaturation = obj.waterCrossSectionArea/obj.area;                 
                 obj.oilSaturation = obj.oilCrossSectionArea/obj.area;
                 end
             else
                 obj.waterCrossSectionArea = obj.area;
                 obj.waterConductance = obj.conductanceSinglePhase; 
                 obj.oilCrossSectionArea = 0;
                 obj.oilConductance = 0; 
                 obj.waterSaturation = 1;                 
                 obj.oilSaturation = 0;
             end
             
        end
        % Patzek_Piri
        function calculateConductance_Drainage_Patzek(obj, Pc)  
             
             if obj.occupancy == 'B' %element was occupied by oil
                 
                 if strcmp(obj.geometry , 'Circle')== 1                     
                 obj.waterCrossSectionArea = 0;    
                 obj.waterSaturation = 0;                 
                 obj.waterConductance = 0;
                 obj.oilCrossSectionArea = obj.area;
                 obj.oilSaturation = 1;
                 obj.oilConductance = obj.area^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                 
                 else
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 cornerArea = zeros(1,4);
                 cornerConductance = zeros(1,4);
                 %Raduis of Curvature
                 Rc = obj.sig_ow / abs(Pc);
                 for jj = 1:4
                     if ~isnan(halfAngles(jj))                       
                         obj.b(jj) = obj.sig_ow / Pc * ...
                             cos(obj.recedingContactAnglee + halfAngles(jj)) * sin(halfAngles(jj));  
                         % Area
                         % Based on Piri_2005: eq A4 & A5 
                         if (halfAngles(jj) + obj.recedingContactAnglee) == pi/2
                             cornerArea(jj) = (Rc * cos(obj.recedingContactAnglee + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             cornerArea(jj) = Rc ^2 * (cos(obj.recedingContactAnglee)*...
                                 (cot(halfAngles(jj)) * cos(obj.recedingContactAnglee) - sin(obj.recedingContactAnglee))+ ...
                                 obj.recedingContactAnglee + halfAngles(jj) - pi/2);
                         end
                         %Conductance
                         % Based on Piri_2005: eq B(10 - 15)
                         f = 1; % no-flow boundary condition suitable for oil-water interfaces
                         F1 = pi/2 - halfAngles(jj) - obj.recedingContactAnglee;
                         F2 = cot(halfAngles(jj)) * cos(obj.recedingContactAnglee) - sin(obj.recedingContactAnglee); 
                         F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));
                         
                         if (obj.recedingContactAnglee <= pi/2 - halfAngles(jj))  % Positive Curvature                       
                             cornerConductance(jj) = (cornerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                                 (F2 * cos(obj.recedingContactAnglee) - F1) * F3 ^ 2) / ...
                                 (12 * obj.waterViscosity * ((sin(halfAngles(jj))) * ...
                                 (1 - F3) * (F2 + f * F1))^ 2);
                         elseif (obj.recedingContactAnglee > pi/2 - halfAngles(jj)) % Negative Curvature
                             cornerConductance(jj) = (cornerArea(jj)^2 * tan(halfAngles(jj))* ...
                                 (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                                 (12 * obj.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
                         end
                     end
                     obj.waterCrossSectionArea = sum(cornerArea); 
                     obj.waterConductance = sum(cornerConductance);
                 end
                 obj.oilCrossSectionArea = obj.area - obj.waterCrossSectionArea; 
                 if strcmp(obj.geometry , 'Triangle')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                 end                 
                 obj.waterSaturation = obj.waterCrossSectionArea/obj.area;                 
                 obj.oilSaturation = obj.oilCrossSectionArea/obj.area;
                 end
             else
                 obj.waterCrossSectionArea = obj.area;
                 obj.waterConductance = obj.conductanceSinglePhase; 
                 obj.oilCrossSectionArea = 0;
                 obj.oilConductance = 0; 
                 obj.waterSaturation = 1;                 
                 obj.oilSaturation = 0;
             end
             
        end 
         % Valvatne
        function calculateConductance_Drainage_V(obj, Pc)  
            
             if obj.occupancy == 'B' %element was occupied by oil
                 
                 if strcmp(obj.geometry , 'Circle')== 1                     
                 obj.waterCrossSectionArea = 0;    
                 obj.waterSaturation = 0;                 
                 obj.waterConductance = 0;
                 obj.oilCrossSectionArea = obj.area;
                 obj.oilSaturation = 1;
                 obj.oilConductance = obj.area^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                 else
                     
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 cornerArea = zeros(1,4);
                 cornerConductance = zeros(1,4);
                  
                 for jj = 1:4
                     if ~isnan(halfAngles(jj)) 
                         
                         obj.b(jj) = obj.sig_ow / Pc * ...
                             cos(obj.recedingContactAnglee + halfAngles(jj)) * sin(halfAngles(jj)); 
                         % Area
                          % Based on Valvatne 3.45-3.49   
                         if obj.recedingContactAnglee + halfAngles(jj) < pi/2
                             cornerArea(jj) =  obj.b(jj) ^ 2 * sin(halfAngles(jj))*cos(halfAngles(jj));
                         else
                             cornerArea(jj) = (obj.b(jj)*sin(halfAngles(jj))/cos(halfAngles(jj)+obj.recedingContactAnglee)) ^2 ...
                                 * (cos(obj.recedingContactAnglee)*(cos(halfAngles(jj)+obj.recedingContactAnglee)/ sin(halfAngles(jj)))+ ...
                                 obj.recedingContactAnglee+ halfAngles(jj) - pi/2); 
                         end
                         %Conductance
                         Gstar = sin(halfAngles(jj))*cos(halfAngles(jj))/4/(1+sin(halfAngles(jj)))^2;
                         if obj.recedingContactAnglee + halfAngles(jj) > pi/2
                             Gc = cornerArea(jj)/...
                                 4/(obj.b(jj)*(1-sin(halfAngles(jj))/cos(halfAngles(jj)+obj.recedingContactAnglee)/...
                                 (obj.recedingContactAnglee + halfAngles(jj) - pi/2)))^2;
                         else
                             Gc = Gstar;
                         end                         
                         C = 0.364 + 0.28 * Gstar/Gc;
                         cornerConductance(jj) =  C * cornerArea(jj)^2*Gc/obj.waterViscosity;    
                     end
                     obj.waterCrossSectionArea = sum(cornerArea);
                     obj.waterConductance = sum(cornerConductance);
                 end
                 obj.oilCrossSectionArea = obj.area - obj.waterCrossSectionArea;
                 if strcmp(obj.geometry , 'Triangle')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                 end                 
                 obj.waterSaturation = obj.waterCrossSectionArea/obj.area;                 
                 obj.oilSaturation = obj.oilCrossSectionArea/obj.area;
                 end
             else
                 obj.waterCrossSectionArea = obj.area;
                 obj.waterConductance = obj.conductanceSinglePhase; 
                 obj.oilCrossSectionArea = 0;
                 obj.oilConductance = 0; 
                 obj.waterSaturation = 1;                 
                 obj.oilSaturation = 0;
             end
        end
         % Raeini_Code
        function calculateConductance_Drainage_VR(obj, Pc)  
            
             if obj.occupancy == 'B' %element was occupied by oil
                 
                 if strcmp(obj.geometry , 'Circle')== 1                     
                 obj.waterCrossSectionArea = 0;    
                 obj.waterSaturation = 0;                 
                 obj.waterConductance = 0;
                 obj.oilCrossSectionArea = obj.area;
                 obj.oilSaturation = 1;
                 obj.oilConductance = obj.area^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                 else
                     
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 cornerArea = zeros(1,4);
                 cornerConductance = zeros(1,4);
                  
                 for jj = 1:4
                     if ~isnan(halfAngles(jj))
                         
                         obj.b(jj) = obj.sig_ow / Pc * ...
                             cos(obj.recedingContactAnglee + halfAngles(jj)) * sin(halfAngles(jj));
                         if obj.b(jj) < 0
                                 obj.b(jj) = 0;
                         end
                         % Area                                               
                          % Based on Valvatne 3.45-3.49   
                         if obj.recedingContactAnglee + halfAngles(jj) < pi/2
                             cornerArea(jj) = sin(halfAngles(jj))*cos(halfAngles(jj));
                         else
                             cornerArea(jj) = (sin(halfAngles(jj))/cos(halfAngles(jj)+obj.recedingContactAnglee)) ^2 ...
                                 * (cos(obj.recedingContactAnglee)*(cos(halfAngles(jj)+obj.recedingContactAnglee)/ sin(halfAngles(jj)))+ ...
                                 obj.recedingContactAnglee+ halfAngles(jj) - pi/2); 
                         end
                         %Conductance
                         Gstar = sin(halfAngles(jj))*cos(halfAngles(jj))/4/(1+sin(halfAngles(jj)))^2;
                         if obj.recedingContactAnglee + halfAngles(jj) > pi/2
                             Gc = cornerArea(jj)/...
                                 4/((1-sin(halfAngles(jj))/cos(halfAngles(jj)+obj.recedingContactAnglee)/...
                                 (obj.recedingContactAnglee + halfAngles(jj) - pi/2)))^2;
                         else
                             Gc = Gstar;
                         end                         
                         C = 0.364 + 0.28 * Gstar/Gc;
                         cornerArea(jj) = (cornerArea(jj)* obj.b(jj)) ^ 2;
                         cornerConductance(jj) =  C * cornerArea(jj)^2*Gc  /obj.waterViscosity;                                                  
                     end
                     obj.waterCrossSectionArea = sum(cornerArea);
                     obj.waterConductance = sum(cornerConductance);
                 end
                 obj.oilCrossSectionArea = obj.area - obj.waterCrossSectionArea;
                 if strcmp(obj.geometry , 'Triangle')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                 end                 
                 obj.waterSaturation = obj.waterCrossSectionArea/obj.area;                 
                 obj.oilSaturation = obj.oilCrossSectionArea/obj.area;
                 end
             else
                 obj.waterCrossSectionArea = obj.area;
                 obj.waterConductance = obj.conductanceSinglePhase; 
                 obj.oilCrossSectionArea = 0;
                 obj.oilConductance = 0; 
                 obj.waterSaturation = 1;                 
                 obj.oilSaturation = 0;
             end
        end
        
        %TODO
        % Snap-off for Drainage!!
        
        %% Imbibition
        
        % Piston-Like & oilLayerExist   
        % Based on Patzek 2003
        function calculateThresholdPressurePistonLike_Imbibition_2003(obj, Pc_max_drainage)
         % calculateThresholdPressurePistonLike Summary of this method goes here
         % Detailed explanation goes here  
         
         if strcmp(obj.geometry , 'Circle')== 1
                 % Based on Al-Futaisi&Patzek_2001: eqs 2-5 & Piri_2005: eq C4 
                 obj.imbThresholdPressure_PistonLike = 2*obj.sig_ow *cos(obj.advancingContactAngle)/obj.radius;
                 obj.oilLayerExist(1,:) = nan;
         else
             % Based on  Al-Futaisi&Patzek_2001: eqs 2-4 & 6-10             
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4]; 
             hingingAngles = zeros(1,4);
             b_i = zeros(1,4);
             if strcmp(obj.geometry , 'Triangle')== 1
                 nc = 3;
             else
                 nc = 4;
             end
             nominator = 0; alpha = zeros(nc,1); 
             for i = 1:nc
                 if ~isnan(halfAngles(i))
                     nominator = nominator + cos(obj.recedingContactAngle + halfAngles(i));
                 end
             end
             maxAdvAngle = acos ((-4 * obj.shapeFactor * nominator)/...
                 ((obj.radius * Pc_max_drainage / obj.sig_ow) - cos(obj.recedingContactAngle)+...
                 4 * nc * obj.shapeFactor * sin(obj.recedingContactAngle)));
             rpd = obj.sig_ow / Pc_max_drainage;
             rp = rpd *2;
             if obj.advancingContactAngle <= maxAdvAngle %Spontaneous Imbibibtion
                 rp1 = 2 * rp;
                 rp2 = rp;
                 while abs(rp2 - rp1) > 10^-10
                     rp1 = rp2;  
                     nominator1 = 0; nominator2 =0;
                     for ii = 1:nc
                         if ~isnan(halfAngles(ii))
                             hingingAngles(ii) = acos((rpd / rp1)*cos(obj.recedingContactAngle + halfAngles(ii))) - halfAngles(ii);                            
                             if hingingAngles(ii) <= obj.advancingContactAngle 
                                 b_i(1,ii) = obj.b(ii);
                                 alpha(ii) = asin(b_i(1,ii)/rp1*sin(halfAngles(ii)));                                 
                             else
                                 b_i(1,ii) = rp1 * cos(obj.advancingContactAngle + halfAngles(ii))/sin(halfAngles(ii));
                                 alpha(ii) = pi/2-obj.advancingContactAngle- halfAngles(ii);
                             end
                             if b_i(1,ii) < 0
                                 b_i(1,ii) = 0;
                             end
                             hingingAngles(ii) = min(hingingAngles(ii) , obj.advancingContactAngle);
                             nominator1 = nominator1 + b_i(1,ii)*cos(hingingAngles(ii));
                             nominator2 = nominator2 + (pi/2 - hingingAngles(ii) - halfAngles(ii));
                         end
                     end
                     rp2 = ((obj.radius ^ 2 / 4 / obj.shapeFactor) -rp1*nominator1 + rp1^2 *nominator2) /...
                         (2*rp1 * sum(alpha) + cos(obj.advancingContactAngle)*...
                         ((obj.radius/2/obj.shapeFactor) - 2 * sum(b_i)));
                 end
                 obj.imbThresholdPressure_PistonLike = obj.sig_ow / rp2;
                 obj.oilLayerExist(1,:) = nan;
             elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi/2 + max(halfAngles) %Forced Imbibition
                 obj.imbThresholdPressure_PistonLike = 2 * obj.sig_ow * cos(obj.advancingContactAngle) / obj.radius;
                 obj.oilLayerExist(1,:) = nan;
             elseif obj.advancingContactAngle >= pi/2 + max(halfAngles) %Forced Imbibition 
                 obj.imbThresholdPressure_PistonLike = -calculateThresholdPressurePistonLike_drainage(obj, (pi - obj.advancingContactAngle));                  
             end
         end
        end
        % Based on Patzek 2001 with NR_Khazali
        function calculateThresholdPressurePistonLike_Imbibition(obj, Pc_max_drainage)
       
         if strcmp(obj.geometry , 'Circle')== 1
                 % Based on Al-Futaisi&Patzek_2001: eqs 2-5 & Piri_2005: eq C4 
                 obj.imbThresholdPressure_PistonLike = 2*obj.sig_ow *cos(obj.advancingContactAngle)/obj.radius;
                 obj.oilLayerExist(1,:) = nan;
         else
             % Based on  Al-Futaisi&Patzek_2001: eqs 2-4 & 6-10             
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];  
             if strcmp(obj.geometry , 'Triangle')== 1
                 nc = 3;
             else
                 nc = 4;
             end 
             
             nominator = 0;
             for i = 1:nc
                 if ~isnan(halfAngles(i))
                     nominator = nominator + cos(obj.recedingContactAngle + halfAngles(i));
                 end
             end
             a = (-4 * obj.shapeFactor * nominator)/...
                 ((obj.radius * Pc_max_drainage / obj.sig_ow) - cos(obj.recedingContactAngle)+...
                 12 * obj.shapeFactor * sin(obj.recedingContactAngle));
             if a >1  
                 a = 1; 
             elseif a< -1 
                 a = -1; 
             end
             maxAdvAngle = acos (a); 
             
             if obj.advancingContactAngle <= maxAdvAngle % Spontaneous imbibition
                 
                 maxIteration = 5000;it = 0;
                 rpd = obj.sig_ow / Pc_max_drainage;    
                 NR = 0; 
                 E0 = zeros(nc,1); E1 = zeros(nc,1); alpha = zeros(nc,1);hingingAngles = zeros(1,nc);b_i=zeros(1,nc);
                 rp1 =rpd*2;
                 rp2 = rpd;
                 while abs(rp2 - rp1) > 10^-10 && NR == 0 % fixed point iteration method
                     rp1 = rp2;                                         
                     for ii = 1:nc
                         t = (rpd / rp1)*cos(obj.recedingContactAngle + halfAngles(ii));
                         if t > 1 && t >= -1
                             t = 1;
                         elseif t < -1
                             t = -1;
                         end
                         hingingAngles(ii) = acos(t) - halfAngles(ii);
%                           
                         if ~isnan(halfAngles(ii)) && hingingAngles(ii) <= obj.advancingContactAngle
                             b_i(ii) = rpd * cos(obj.recedingContactAngle + halfAngles(ii))/ sin(halfAngles(ii));
                             part = b_i(ii)/rp1*sin(halfAngles(ii));
                             if part < -1 
                                 part = -1;
                             elseif part > 1 
                                 part = 1;
                             end
                             alpha(ii) = asin(part);
                         elseif ~isnan(halfAngles(ii)) && hingingAngles(ii) > obj.advancingContactAngle                             
                             b_i(ii) = rp1 * cos(obj.advancingContactAngle + halfAngles(ii))/ sin(halfAngles(ii));
                             alpha(ii) = pi/2 - obj.advancingContactAngle - halfAngles(ii);
                         end   
                         if b_i(ii) < 0
                             b_i(ii) = 0;
                         end
                         hingingAngles(ii) = min ( hingingAngles(ii) , obj.advancingContactAngle);
                         E0(ii) = pi/2 - hingingAngles(ii) - halfAngles(ii);
                         E1(ii) = b_i(ii) * cos(hingingAngles(ii)); 
                     end
                     rp2 = (obj.radius ^ 2 / 4 / obj.shapeFactor - rp1 * sum(E1) + rp1^2 * sum(E0))/...
                         (2*rp1 * sum(alpha) + (obj.radius/2/obj.shapeFactor - 2 * sum(b_i)) * cos(obj.advancingContactAngle));
                     
                     R_n = rp2; 
                  
                     if R_n == rp1 %|| abs(rp2-rp1)> 10^-10
                         it = it+1; 
                         NR = 1;
                     end                   
                 end   
                
                 if NR == 1
                 R_o =  rpd * 2;
                 R_n =  rpd;
                 while abs(R_o - R_n) > 10^-10  % NR method based on Khazali 2018
                     dU_dR = 0; dO_dR = 0; dN_dR = 0; dM_dR = 0; M = 0; N = 0; O = 0; U = 0; 
                     R_o = R_n;
                     for ii = 1:nc 
                         %-------------------------------------------------
                         t = rpd/R_o*cos(obj.recedingContactAngle + halfAngles(ii));
                         if t <= 1 && t >= -1
                             teta_H = acos(rpd/R_o*cos(obj.recedingContactAngle + halfAngles(ii)))-halfAngles(ii);
                         else
                             teta_H = obj.recedingContactAngle;
                         end
                         if teta_H <= obj.advancingContactAngle 
                             b_i = rpd * cos(obj.recedingContactAngle + halfAngles(ii)) / sin(halfAngles(ii));
                             t = b_i / R_o * sin( halfAngles(ii));
                             if t <= 1 && t >= -1
                                 alpha = asin (b_i / R_o * sin( halfAngles(ii)));
                             else
                                 alpha = 0;
                             end
                         else
                              b_i = R_o * cos(obj.advancingContactAngle + halfAngles(ii)) / sin(halfAngles(ii));
                              alpha = pi / 2 - obj.advancingContactAngle - halfAngles(ii);
                         end
                         teta_H = min (teta_H , obj.advancingContactAngle);
                         %-------------------------------------------------
                         if teta_H <= obj.advancingContactAngle
                             db_dR = 0;
                             t = (b_i/ R_o * sin(halfAngles(ii)))^2;
                             if t <= 1 && t >= -1
                             dalpha_dR = (R_o * db_dR - b_i) * sin( halfAngles(ii)) / ...
                                 (R_o ^ 2 * sqrt(1-(b_i/ R_o * sin(halfAngles(ii)))^2));
                             else
                                dalpha_dR = 0;
                             end
                         else
                             db_dR = cos(obj.advancingContactAngle + halfAngles(ii))/sin( halfAngles(ii));
                             dalpha_dR = 0;
                         end
                         t = rpd / R_o * cos(obj.recedingContactAngle+halfAngles(ii));
                         if t <= 1 && t >= -1
                             a = acos(rpd / R_o * cos(obj.recedingContactAngle+halfAngles(ii))) - halfAngles(ii);
                             if a >= obj.advancingContactAngle 
                                 dtetaH_dR = 0;
                             else
                                 t = (rpd/R_o * cos(obj.recedingContactAngle+halfAngles(ii)))^2;
                                 if t <= 1 && t >= -1
                                 dtetaH_dR = rpd * cos(obj.recedingContactAngle+halfAngles(ii))/...
                                     (R_o ^ 2 * sqrt(1-(rpd/R_o * cos(obj.recedingContactAngle+halfAngles(ii)))^2));
                                 else
                                    dtetaH_dR = rpd * cos(obj.recedingContactAngle+halfAngles(ii))/(R_o ^ 2 );
                                 end
                             end
                         else
                             dtetaH_dR = 0;
                         end
                         %-------------------------------------------------
                         dU_dR = dU_dR - 2 * cos(obj.advancingContactAngle) * db_dR;
                         dO_dR = dO_dR + 2 * alpha + 2 * R_o * dalpha_dR;
                         dN_dR = dN_dR + 2 * R_o * (pi/2 - teta_H - halfAngles(ii)) - R_o ^2 * dtetaH_dR;
                         dM_dR = dM_dR + R_o * (dtetaH_dR * b_i * sin(teta_H) - db_dR * cos(teta_H)) - b_i * cos(teta_H); 
                         M = M - R_o * b_i * cos(teta_H) ;
                         N = N + R_o ^ 2 * (pi/2 - teta_H - halfAngles(ii));
                         O = O + 2 * R_o * alpha;
                         U = U + (obj.radius /2 / obj.shapeFactor - 2 * b_i ) * cos(obj.advancingContactAngle);
                     end                     
                     A = obj.radius ^ 2 / 4 / obj.shapeFactor;
                     F_R = (A + M + N) / (O + U)- R_o;
                     dF_R = ((dM_dR+dN_dR)*(O + U)-(dO_dR + dU_dR)*(M + N + A)) / (O + U)^2 - 1;
                     R_n = R_o - F_R / dF_R; 
                 end 
                 if it > maxIteration  
                     err = 2 * abs(R_o - R_n) / (abs(R_o) + abs(R_n)+0.001); 
                     fprintf('err %f\n',err); 
                 end
                 end
                 obj.imbThresholdPressure_PistonLike = obj.sig_ow / R_n; 
             elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi/2 + max(halfAngles) % Forced imbibition
                 obj.imbThresholdPressure_PistonLike = 2 * obj.sig_ow * cos(obj.advancingContactAngle) / obj.radius; 
             elseif obj.advancingContactAngle >= pi/2 + max(halfAngles) % Forced imbibition
                 obj.imbThresholdPressure_PistonLike = -calculateThresholdPressurePistonLike_drainage(obj, (pi - obj.advancingContactAngle)); 
             end 
         end
        end 
        % Based on Raeini 2019 with NR
        function calculateThresholdPressurePistonLike_Imbibition_R(obj, Pc_max_drainage)
       
         if strcmp(obj.geometry , 'Circle')== 1
                 % Based on Al-Futaisi&Patzek_2001: eqs 2-5 & Piri_2005: eq C4 
                 obj.imbThresholdPressure_PistonLike = 2*obj.sig_ow *cos(obj.advancingContactAngle)/obj.radius;
                 obj.oilLayerExist(1,:) = nan;
         else
             % Based on  Al-Futaisi&Patzek_2001: eqs 2-4 & 6-10             
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];  
             if strcmp(obj.geometry , 'Triangle')== 1
                 nc = 3;
             else
                 nc = 4;
             end 
             
             nominator = 0;
             for i = 1:nc
                 if ~isnan(halfAngles(i))
                     nominator = nominator + cos(obj.recedingContactAngle + halfAngles(i));
                 end
             end
             a = (-4 * obj.shapeFactor * nominator)/...
                 ((obj.radius * Pc_max_drainage / obj.sig_ow) - cos(obj.recedingContactAngle)+...
                 12 * obj.shapeFactor * sin(obj.recedingContactAngle));
             if a >1  
                 a = 1; 
             elseif a< -1 
                 a = -1; 
             end
             maxAdvAngle = acos (a);  
             if obj.advancingContactAngle <= maxAdvAngle % Spontaneous imbibition  
                     newPc = 1.1 * obj.sig_ow * 2 * cos(obj.advancingContactAngle) / obj.radius; 
                     MAX_NEWT_ITR = 50000;
                     for itr = 1:MAX_NEWT_ITR
                         sumOne = 0; sumTwo = 0; sumThree = 0; sumFour = 0;
                         oldPc = newPc;
                         for i = 1:nc  
                             hingTeta = obj.advancingContactAngle;
                             b_i = obj.sig_ow / oldPc * cos(hingTeta +  halfAngles(i))/sin(halfAngles(i)); 
                             if b_i< 0
                                 b_i = 0;
                             end
                             partus = b_i * sin(halfAngles(i)) * oldPc / obj.sig_ow;
                             if partus > 1 
                                 partus = 1;
                             elseif partus < -1
                                 partus = -1;
                             end
                             sumOne = sumOne + b_i * cos(hingTeta);
                             sumTwo = sumTwo + pi / 2 - hingTeta - halfAngles(i);
                             sumThree = sumThree + asin(partus);
                             sumFour = sumFour + b_i;
                         end
                         a = 2 * sumThree - sumTwo;
                         bb = cos(obj.advancingContactAngle) * obj.radius / 2 / obj.shapeFactor - 2 * sumFour + sumOne;
                         c = -obj.area;
                         root_s = bb ^ 2 - 4 * a  * c;
                         if root_s > 0
                             newPc2 = obj.sig_ow * 2 * a / ( -bb + sqrt(root_s));
                         else
                             newPc2 = obj.sig_ow * 2 * a / (-bb);
                         end
                         newPc = newPc2;
                         err = 2 * abs(newPc - oldPc) / (abs(oldPc) + abs(newPc)+0.001);
                         if err < 10^-11
                             break;
                         end
                     end
                     obj.imbThresholdPressure_PistonLike = newPc;
                     if err <0.1 && err > 0.0001
                         fprintf('err %f\n',obj.index);
                     end                     
             elseif obj.advancingContactAngle < pi/2 + max(halfAngles) % Forced imbibition
                 obj.imbThresholdPressure_PistonLike = 2 * obj.sig_ow * cos(obj.advancingContactAngle) / obj.radius;
             elseif obj.advancingContactAngle >= pi/2 + max(halfAngles) % Forced imbibition
                 obj.imbThresholdPressure_PistonLike = -calculateThresholdPressurePistonLike_drainage(obj, (pi - obj.advancingContactAngle));
             end   
         end
        end  
        
        % oilLayer existance check        
        function oilLayerExistance(obj)
       
         if strcmp(obj.geometry , 'Circle')== 1 
                 obj.oilLayerExist(1,:) = nan;
         else   
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];             
             for i = 1:4
                 if ~isnan(halfAngles(i)) && obj.advancingContactAngle >= pi/2 + halfAngles(i)
                     obj.oilLayerExist(1,i) = 1;
                 end
             end
         end
        end  
        
        % TODO 
        % CornerAPEX
        
        % Snap-off 
        % Patzek
        function calculateThresholdPressureSnapOff_Patzek(obj,Pc_max_drainage)
            
         if strcmp(obj.geometry , 'Circle')== 1
             obj.imbThresholdPressure_SnapOff = nan;
         else
             % Based on  Al-Futaisi&Patzek_2003: eqs 12-14          
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];             
             maxAdvAngle = pi/2 - min(halfAngles);
             if strcmp(obj.geometry , 'Triangle')== 1
                 nc = 3; 
             if obj.advancingContactAngle < maxAdvAngle %Spontaneous Imbibition
                 
                 rso = zeros(1,nc);
                 for edge = 1:nc   
                     
                     rso1 = obj.sig_ow / obj.imbThresholdPressure_PistonLike; 
                     rso2 = rso1*2;
                     
                     while abs(rso1 - rso2) > 10^-8
                         rso1 = rso2;
                         hingeAngle_ii = acos((obj.sig_ow / rso1)*...
                             cos(obj.recedingContactAngle + halfAngles(1))/Pc_max_drainage) - halfAngles(1); 
                         if hingeAngle_ii <= obj.advancingContactAngle
                             E1_i = cos(obj.recedingContactAngle + halfAngles(1))/sin(halfAngles(1));
                         else
                             E1_i = cos(obj.advancingContactAngle + halfAngles(1))/sin(halfAngles(1));
                         end                             
                         
                         hingeAngle_jj = acos((obj.sig_ow / rso1)*...
                             cos(obj.recedingContactAngle + halfAngles(2)) / Pc_max_drainage) - halfAngles(2); 
                         if hingeAngle_jj <= obj.advancingContactAngle
                             E1_j = cos(obj.recedingContactAngle + halfAngles(2))/sin(halfAngles(2));
                         else
                             E1_j = cos(obj.advancingContactAngle + halfAngles(2))/sin(halfAngles(2));
                         end 
                                                  
                         hingeAngle_kk = acos((obj.sig_ow / rso1)*...
                             cos(obj.recedingContactAngle + halfAngles(3)) / Pc_max_drainage) - halfAngles(3); 
                         if hingeAngle_kk <= obj.advancingContactAngle
                             E1_k = cos(obj.recedingContactAngle + halfAngles(3))/sin(halfAngles(3));
                         else
                             E1_k = cos(obj.advancingContactAngle + halfAngles(3))/sin(halfAngles(3));
                         end  
                         
                         rso12 = obj.radius *(cot(halfAngles(1))+ cot(halfAngles(2))) / (E1_i + E1_j);
                         rso23 = obj.radius *(cot(halfAngles(2))+ cot(halfAngles(3))) / (E1_j + E1_k);
                         rso31 = obj.radius *(cot(halfAngles(3))+ cot(halfAngles(1))) / (E1_k + E1_i);
                         r2 = min(rso12,rso23);
                         rso = min(r2,rso31);
                     end
                 end
                 obj.imbThresholdPressure_SnapOff = obj.sig_ow / rso;
                 
                 % Forced imbibition part
             elseif obj.advancingContactAngle == maxAdvAngle
                 obj.imbThresholdPressure_SnapOff = 0;
             elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
                 obj.imbThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + min(halfAngles))/...
                     cos(obj.recedingContactAngle + min(halfAngles)); 
             elseif obj.advancingContactAngle >= pi - min(halfAngles)
                 obj.imbThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
             end
             else 
                 % elemnt is square : Piri
                 if obj.advancingContactAngle <= pi/4
                     %eq C34
                     obj.imbThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cot(pi/4)*cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     obj.imbThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(obj.recedingContactAngle + pi/4); 
                 elseif obj.advancingContactAngle > 3*pi/4
                     obj.imbThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end  
             end
         end         
        end
        % Valvatne        
        function calculateThresholdPressureSnapOff_V(obj,Pc_max_drainage)
         if strcmp(obj.geometry , 'Circle')== 1
             obj.imbThresholdPressure_SnapOff = nan;
         else
             % Based on  Al-Futaisi&Patzek_2003: eqs 12-14          
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];             
             maxAdvAngle = pi/2 - min(halfAngles);
             if strcmp(obj.geometry , 'Triangle')== 1
                 if obj.advancingContactAngle < maxAdvAngle %Spontaneous Imbibition
                     Pc_a = obj.sig_ow / obj.radius *(cos(obj.advancingContactAngle)-...
                         2*sin(obj.advancingContactAngle)/(cot(obj.halfAngle1)+cot(obj.halfAngle2)));
                     Pc_b = obj.sig_ow / obj.radius *...
                         (cos(obj.advancingContactAngle)*cot(obj.halfAngle1)-sin(obj.advancingContactAngle)+...
                         cos(obj.recedingContactAngle)*cot(obj.halfAngle3)-sin(obj.recedingContactAngle))/...
                         ((cot(obj.halfAngle1)+cot(obj.halfAngle2)));                 
                     obj.imbThresholdPressure_SnapOff = max(Pc_a,Pc_b);                 
                     % Forced imbibition part
                 elseif obj.advancingContactAngle == maxAdvAngle
                     obj.imbThresholdPressure_SnapOff = 0;
                 elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
                     obj.imbThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + min(halfAngles))/...
                         cos(obj.recedingContactAngle + min(halfAngles));
                 elseif obj.advancingContactAngle >= pi - min(halfAngles)
                     obj.imbThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end
             else 
                 % elemnt is square
                 if obj.advancingContactAngle <= pi/4
                     %eq C34
                     obj.imbThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cot(pi/4)*cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     obj.imbThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(obj.recedingContactAngle + pi/4);
                     if obj.imbThresholdPressure_SnapOff > 0
                         obj.imbThresholdPressure_SnapOff = -1 *obj.imbThresholdPressure_SnapOff;
                     end
                 elseif obj.advancingContactAngle > 3*pi/4
                     obj.imbThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end  
             end
         end         
        end
        % Raeini TODO       
        function calculateThresholdPressureSnapOff_Raeini(obj,Pc_max_drainage)
         if strcmp(obj.geometry , 'Circle')== 1
             obj.imbThresholdPressure_SnapOff = nan;
         else
             % Based on  Al-Futaisi&Patzek_2003: eqs 12-14          
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];             
             maxAdvAngle = pi/2 - min(halfAngles);
             if strcmp(obj.geometry , 'Triangle')== 1
                 if obj.advancingContactAngle < maxAdvAngle %Spontaneous Imbibition
%                      b_2 = obj.sig_ow / Pc_max_drainage * cos(obj.recedingContactAngle+obj.halfAngle1)/sin(obj.halfAngle1);
%                      snapOffPrsOne = obj.sig_ow *(cos(obj.advancingContactAngle+obj.halfAngle1)...
%                          / ((obj.radius/tan(obj.halfAngle1)+obj.radius/tan(obj.halfAngle3)-b_3)*sin(obj.halfAngle1)));
%                      oldPc = Pc_max_drainage; errorVal = 1;
%                      L0pL3 = (obj.radius/tan(obj.halfAngle1)+obj.radius/tan(obj.halfAngle2));
%                      for i= 1:1000
%                          apexDist3 = 0; teta3 = obj.advancingContactAngle; 
%                          b_1 = ;
                     Pc_b = obj.sig_ow / obj.radius *...
                         (cos(obj.advancingContactAngle)*cot(obj.halfAngle1)-sin(obj.advancingContactAngle)+...
                         cos(obj.recedingContactAngle)*cot(obj.halfAngle3)-sin(obj.recedingContactAngle))/...
                         ((cot(obj.halfAngle1)+cot(obj.halfAngle2)));                 
                     obj.imbThresholdPressure_SnapOff = max(Pc_a,Pc_b);                 
                     % Forced imbibition part
                 elseif obj.advancingContactAngle == maxAdvAngle
                     obj.imbThresholdPressure_SnapOff = 0;
                 elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
                     obj.imbThresholdPressure_SnapOff = Pc_max_drainage*cos(min(halfAngles))/...
                         cos(min(halfAngles)); 
                 elseif obj.advancingContactAngle >= pi - min(halfAngles)
                     obj.imbThresholdPressure_SnapOff = -Pc_max_drainage/cos( min(halfAngles));
                 end
             else 
                 % elemnt is square
                 if obj.advancingContactAngle <= pi/4
                     %eq C34
                     obj.imbThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cot(pi/4)*cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     obj.imbThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(pi/4);
                     if obj.imbThresholdPressure_SnapOff > 0
                         obj.imbThresholdPressure_SnapOff = -1 *obj.imbThresholdPressure_SnapOff;
                     end
                 elseif obj.advancingContactAngle > 3*pi/4
                     obj.imbThresholdPressure_SnapOff = -Pc_max_drainage/cos(min(halfAngles));
                 end  
             end
         end         
        end
        % Zolfaghari 
        function calculateThresholdPressureSnapOff(obj,Pc_max_drainage)
            
         if strcmp(obj.geometry , 'Circle')== 1
             obj.imbThresholdPressure_SnapOff = nan;
         else
             r_dr = obj.sig_ow / Pc_max_drainage;
             % Based on  Zolfaghari_2014: eqs 4.31    
             if strcmp(obj.geometry , 'Triangle')== 1 
                 A = zeros(1,7);
                 if obj.advancingContactAngle <= pi/2 - obj.halfAngle1 %Spontaneous Imbibition
                     A(1) = cos(obj.advancingContactAngle) - 2 * sin(obj.advancingContactAngle)/ ...
                         (cot(obj.halfAngle1)+cot(obj.halfAngle2));
                     A(2) = cos(obj.advancingContactAngle) - 2 * sin(obj.advancingContactAngle)/ ...
                         (cot(obj.halfAngle1)+cot(obj.halfAngle3));
                     A(3) = cos(obj.advancingContactAngle) - 2 * sin(obj.advancingContactAngle)/ ...
                         (cot(obj.halfAngle2)+cot(obj.halfAngle3));
                     A(4) = (cos(obj.advancingContactAngle)*cot(obj.halfAngle1)-sin(obj.advancingContactAngle)) / ...
                         (cot(obj.halfAngle1)+cot(obj.halfAngle3)- r_dr * obj.radius * ...
                         cos(obj.recedingContactAngle + obj.halfAngle3)/sin(obj.halfAngle3));
                     A(5) = (cos(obj.advancingContactAngle)*cot(obj.halfAngle2)-sin(obj.advancingContactAngle)) / ...
                         (cot(obj.halfAngle2)+cot(obj.halfAngle3)- r_dr * obj.radius * ...
                         cos(obj.recedingContactAngle + obj.halfAngle3)/sin(obj.halfAngle3));
                     if obj.advancingContactAngle > pi/2 - obj.halfAngle2
                         A(6) = Pc_max_drainage * cos(obj.advancingContactAngle + obj.halfAngle2) / ...
                             cos(obj.recedingContactAngle + obj.halfAngle2);
                     elseif obj.advancingContactAngle > pi/2 - obj.halfAngle3
                         A(7) = Pc_max_drainage * cos(obj.advancingContactAngle + obj.halfAngle3) / ...
                             cos(obj.recedingContactAngle + obj.halfAngle3);
                     end
                     max_A = max(A);
                     obj.imbThresholdPressure_SnapOff = obj.sig_ow /obj.radius * max_A; 
                 else  % Forced imbibition part
                     obj.imbThresholdPressure_SnapOff = Pc_max_drainage * cos(obj.advancingContactAngle + obj.halfAngle1) / ...
                             cos(obj.recedingContactAngle + obj.halfAngle1);
                 end
             
             else % elemnt is square 
                 if obj.advancingContactAngle <= pi/4 %Spontaneous Imbibition 
                     obj.imbThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));  
                     % Forced imbibition part
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     obj.imbThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(obj.recedingContactAngle + pi/4); 
                 elseif obj.advancingContactAngle > 3*pi/4
                     obj.imbThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + pi/4);
                 end                
             end
         end         
        end 
        
        % LayerCollapse Piri 
        % Piri
        function imbThresholdPressure_LayerCollapse = calculateThresholdPressureLayerCollapse_Piri(obj, Pc_max_drainage)
         
         if strcmp(obj.geometry , 'Circle')== 1
             imbThresholdPressure_LayerCollapse = nan(1,4);
         else
              halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];  
              if strcmp(obj.geometry , 'Triangle')== 1
                  nc = 3;
              elseif strcmp(obj.geometry , 'Square')== 1
                  nc = 4;
              end
              imbThresholdPressure_LayerCollapse = nan(1, 4); hingingAngles = zeros(1,4);   b_i = zeros(nc , 1);
              for i = 1:nc 
                  hingingAngles(i) = obj.advancingContactAngle;
                  h = hingingAngles(i)/2;
                  b_i(i) = obj.sig_ow / obj.imbThresholdPressure_PistonLike * cos(hingingAngles(i) + halfAngles(i))/sin(halfAngles(i));
                  while abs(hingingAngles(i) - h) > 10 ^ -5  
                      Pc_n = obj.sig_ow *...
                          (3*(sin(halfAngles(i)))^2+ 4*sin(halfAngles(i))*cos(hingingAngles(i))+(cos(hingingAngles(i)))^2)/...
                          (b_i(i)*(cos(halfAngles(i))*sin(halfAngles(i))*(2*sin(halfAngles(i))+cos(hingingAngles(i)))+...
                          (sin(halfAngles(i)))^2 * ...
                          sqrt(4*(cos(halfAngles(i)))^2-3-(cos(hingingAngles(i)))^2-4*sin(halfAngles(i))*cos(hingingAngles(i)))));                       
                      h = hingingAngles(i);
                      hingingAngles(i) = acos((Pc_n/Pc_max_drainage)*cos(obj.recedingContactAngle + halfAngles(i))) - halfAngles(i);  
                      hingingAngles(i) = min (hingingAngles(i), obj.advancingContactAngle);
                      b_i(i) = obj.sig_ow / Pc_n * cos(hingingAngles(i) + halfAngles(i))/sin(halfAngles(i));
                  end
                  imbThresholdPressure_LayerCollapse(i) = Pc_n;
              end
              
         end
        end        
        % Zolfaghari
        function imbThresholdPressure_LayerCollapse = calculateThresholdPressureLayerCollapse(obj, Pc_max_drainage)
         
         if strcmp(obj.geometry , 'Circle')== 1
             imbThresholdPressure_LayerCollapse = nan(1,4);
         else
              halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];  
              if strcmp(obj.geometry , 'Triangle')== 1
                  nc = 3;
              elseif strcmp(obj.geometry , 'Square')== 1
                  nc = 4;
              end
              imbThresholdPressure_LayerCollapse = nan(1, 4); hingingAngles = zeros(1,4);   b_i = zeros(nc , 1);
              for i = 1:nc   
                  b_i(i) = obj.sig_ow / Pc_max_drainage * cos(hingingAngles(i) + halfAngles(i))/sin(halfAngles(i));  
                  imbThresholdPressure_LayerCollapse(i) = obj.sig_ow / b_i(i)*...
                      (cot(halfAngles(i))+ (2*sin(halfAngles(i))*cos(obj.advancingContactAngle))- ...
                      sqrt((sin(obj.advancingContactAngle))^2-4*sin(halfAngles(i))^2-4*sin(halfAngles(i))*cos(obj.advancingContactAngle))); 
              end
              
         end
        end
        
        % Conductance Calculation _ Imbibition 
        % Based on Valvatne eqs.
        function  calculateConductance_Imbibition_V(obj, network, Pc)
             Pc = abs(Pc);
             
           if ~any(obj.oilLayerExist)
               
            if obj.occupancy == 'B' 
                
             if strcmp(obj.geometry , 'Circle')== 1
                  obj.waterCrossSectionArea = 0;
                  obj.waterConductance = 0;
                  obj.oilCrossSectionArea = obj.area;
                  obj.oilConductance = obj.oilCrossSectionArea^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                  obj.waterSaturation = 0;
                  obj.oilSaturation = 1;
             else
                 if strcmp(obj.geometry , 'Triangle')== 1
                     nc = 3;
                 else
                     nc = 4;
                 end
                 
                 R = network.sig_ow / Pc; %Raduis of Curvature                 
                 R_min = network.sig_ow / network.Pc_drain_max;
                 
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 cornerArea = zeros(1,nc);
                 cornerConductance = zeros(1,nc);
                 for jj = 1:nc
                     part = R_min * cos(obj.recedingContactAngle+halfAngles(jj))/R; 
                     obj.hingeAngles(jj) = acos(part)-halfAngles(jj);
                     obj.hingeAngles(jj) = min(obj.hingeAngles(jj), obj.advancingContactAngle); 
                     obj.b(jj) = R * cos(obj.hingeAngles(jj) + halfAngles(jj))/ sin(halfAngles(jj));
                     if obj.b(jj) < 0
                         obj.b(jj) = 0;
                     else
                     
                     % Area
                     % Based on Valvatne 3.45-3.49   
                     if obj.hingeAngles(jj) + halfAngles(jj) < pi/2
                         cornerArea(jj) =  obj.b(jj) ^ 2 * sin(halfAngles(jj))*cos(halfAngles(jj));
                     else
                         cornerArea(jj) = (obj.b(jj)*sin(halfAngles(jj))/cos(halfAngles(jj)+obj.hingeAngles(jj))) ^2 ...
                             * (cos(obj.hingeAngles(jj))*(cos(halfAngles(jj)+obj.hingeAngles(jj))/ sin(halfAngles(jj)))+ ...
                             obj.hingeAngles(jj) + halfAngles(jj) - pi/2); 
                     end
                     %Conductance
                     Gstar = sin(halfAngles(jj))*cos(halfAngles(jj))/4/(1+sin(halfAngles(jj)))^2;
                     if obj.hingeAngles(jj)  + halfAngles(jj) > pi/2
                         Gc = cornerArea(jj)/...
                             4/(obj.b(jj)*(1-sin(halfAngles(jj))/cos(halfAngles(jj)+obj.hingeAngles(jj))/...
                             (obj.hingeAngles(jj)  + halfAngles(jj) - pi/2)))^2;
                     else
                         Gc = Gstar;
                     end
                     C = 0.364 + 0.28 * Gstar/Gc;
                     cornerConductance(jj) =  C * cornerArea(jj)^2*Gc/obj.waterViscosity;                      
                     end 
                 end
                 
                 obj.waterConductance = sum(cornerConductance);
                 obj.waterCrossSectionArea = sum(cornerArea);
                 obj.oilCrossSectionArea = obj.area - obj.waterCrossSectionArea;
                 
                 if strcmp(obj.geometry , 'Triangle')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     obj.oilConductance = obj.oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                 end                 
                 obj.waterSaturation = obj.waterCrossSectionArea / obj.area;
                 obj.oilSaturation = obj.oilCrossSectionArea / obj.area;
             end             
            else
                 obj.waterCrossSectionArea = obj.area;
                 obj.waterConductance = obj.conductanceSinglePhase; 
                 obj.oilCrossSectionArea = 0;
                 obj.oilConductance = 0; 
                 obj.waterSaturation = 1;
                 obj.oilSaturation = 0;
            end
            
           else
                 if strcmp(obj.geometry , 'Triangle')== 1
                     nc = 3;
                 else
                     nc = 4;
                 end
                 
                 R = network.sig_ow / Pc; %Raduis of Curvature
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 inerArea = zeros(1,nc);
                 outerArea = zeros(1,nc);
                 layerArea = zeros(1,nc);
                 inerConductance = zeros(1,nc);
                 layerConductance = zeros(1,nc);
                 for jj = 1:nc
                     if ~isnan(obj.oilLayerExist(1,jj)) %if the layer exist in the corner
                         
                         R_min = network.sig_ow / network.Pc_drain_max;
                         part = R_min * cos(obj.recedingContactAngle+halfAngles(jj))/R; 
                         thetaHingRec = acos(part)-halfAngles(jj);  
                         thetaHingRec = min(thetaHingRec, obj.advancingContactAngle);
                         
                         % Area of corner water
                         if (halfAngles(jj) + thetaHingRec) == pi/2
                             inerArea(jj) = (R * cos(thetaHingRec + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             inerArea(jj) = R ^2 * (cos(thetaHingRec)*...
                                 (cot(halfAngles(jj)) * cos(thetaHingRec) - sin(thetaHingRec))+ ...
                                 thetaHingRec+ halfAngles(jj) - pi/2);
                         end
                         % Conductance of corner water
                         if thetaHingRec <= pi/2 - halfAngles(jj) % Positive Curvature       
                             f = 1; % no-flow boundary condition suitable for oil-water interfaces
                             F1 = pi/2 - halfAngles(jj) - thetaHingRec;
                             F2 = cot(halfAngles(jj)) * cos(thetaHingRec) - sin(thetaHingRec); 
                             F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                             
                             inerConductance(jj) = (inerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                                 (F2 * cos(thetaHingRec) - F1) * F3 ^ 2) / ...
                                 (12 * network.waterViscosity * ((sin(halfAngles(jj))) * ...
                                 (1 - F3) * (F2 + f * F1))^ 2);
                         elseif (thetaHingRec > pi/2 - halfAngles(jj)) % Negative Curvature
                             inerConductance(jj) = (inerArea(jj)^2 * tan(halfAngles(jj))* ...
                                 (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                                 (12 * network.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
                         end
                         thetaHingAdv = pi - obj.advancingContactAngle;
                         if (halfAngles(jj) + thetaHingAdv) == pi/2
                             outerArea(jj) = (R * cos(thetaHingAdv + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             outerArea(jj) = R ^2 * (cos(thetaHingAdv)*...
                                 (cot(halfAngles(jj)) * cos(thetaHingAdv) - sin(thetaHingAdv))+ ...
                                 thetaHingAdv+ halfAngles(jj) - pi/2);
                         end
                         % Area of oil layer
                         layerArea(jj) = outerArea(jj)-inerArea(jj);
                         % Conductance of oil layer
                         F3_layer = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                         layerConductance(jj) = (layerArea(jj)^3 * (1 - sin(halfAngles(jj)))^2 * ...
                             tan(halfAngles(jj)) * F3_layer ^ 2) / ...
                             (12 * network.oilViscosity *outerArea(jj)* (sin(halfAngles(jj)))^2 * ...
                             (1 - F3_layer) * (1 + F3_layer - (1- F3) * sqrt(inerArea(jj)/outerArea(jj))));
                     end
                 end
                 % Center water area and conductance
                 centerWaterArea = obj.area - sum(outerArea);
                 if strcmp(obj.geometry , 'Triangle')== 1
                     centerWaterConductance = 3 *obj.radius^2*centerWaterArea/20/network.waterViscosity;
                 else                     
                     centerWaterConductance = 0.5623 * obj.shapeFactor * centerWaterArea^2/network.waterViscosity;
                 end
                 
                 obj.waterCrossSectionArea = obj.area - sum(layerArea);
                 obj.waterConductance = centerWaterConductance + sum(inerConductance);
                 obj.oilCrossSectionArea = sum(layerArea);
                 obj.oilConductance = sum(layerConductance);
                 obj.waterSaturation = obj.waterCrossSectionArea / obj.area;
                 obj.oilSaturation = obj.oilCrossSectionArea / obj.area;
           end           
            if obj.waterSaturation > 1
                % Control 
                fprintf('sat %f %4.0d  %f \n',Pc, obj.index,obj.waterSaturation);  
                obj.control(1:3) = [Pc, obj.imbThresholdPressure_PistonLike, obj.imbThresholdPressure_SnapOff];
                
                     obj.waterCrossSectionArea = obj.area;
                     obj.waterConductance = obj.conductanceSinglePhase; 
                     obj.oilCrossSectionArea = 0;
                     obj.oilConductance = 0; 
                     obj.waterSaturation = 1;
                     obj.oilSaturation = 0; 
            end
         end   
        % Based on Valvatne & Raeini
        function  calculateConductance_Imbibition(obj, network, Pc)
             Pc = abs(Pc);
             
           if ~any(obj.oilLayerExist)
               
            if obj.occupancy == 'B' 
                
             if strcmp(obj.geometry , 'Circle')== 1
                  obj.waterCrossSectionArea = 0;
                  obj.waterConductance = 0;
                  obj.oilCrossSectionArea = obj.area;
                  obj.oilConductance = obj.oilCrossSectionArea^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                  obj.waterSaturation = 0;
                  obj.oilSaturation = 1;
             else
                 if strcmp(obj.geometry , 'Triangle')== 1
                     nc = 3;
                 else
                     nc = 4;
                 end
%                  conAng = obj.advancingContactAngle;
                 R = network.sig_ow / Pc; %Raduis of Curvature                 
                 R_min = network.sig_ow / network.Pc_drain_max;
                 
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 cornerArea = zeros(1,nc);
                 cornerConductance = zeros(1,nc);
                 for jj = 1:nc
                     
                     part = R_min * cos(obj.recedingContactAngle+halfAngles(jj))/R; 
                     obj.hingeAngles(jj) = acos( part)-halfAngles(jj);
                     obj.hingeAngles(jj) = min(obj.hingeAngles(jj), obj.advancingContactAngle);  
                     obj.b(jj) = R * cos(obj.hingeAngles(jj) + halfAngles(jj))/ sin(halfAngles(jj)); 
                     
                     if obj.b(jj) < 0
                         obj.b(jj) = 0;
                     else 
                                          
                     % Area
                     % Based on Valvatne 3.45-3.49   
                     if abs(obj.hingeAngles(jj) + halfAngles(jj) - pi/2) < 0.01
                         cornerArea(jj) = sin(halfAngles(jj))*cos(halfAngles(jj));
                     else
                         cornerArea(jj) = (sin(halfAngles(jj))/cos(halfAngles(jj)+obj.hingeAngles(jj))) ^2 ...
                             * (cos(obj.hingeAngles(jj))*cos(halfAngles(jj)+obj.hingeAngles(jj))/ sin(halfAngles(jj))+ ...
                             obj.hingeAngles(jj) + halfAngles(jj) - pi/2); 
                     end
                     %Conductance
                     Gstar = (sin(halfAngles(jj))*cos(halfAngles(jj)))/(4*(1+sin(halfAngles(jj)))^2);
                     Gc = Gstar;
                     if abs(obj.hingeAngles(jj) + halfAngles(jj) - pi/2) > 0.01
                         Gc = cornerArea(jj)/...
                             (4*(1-(sin(halfAngles(jj))/cos(halfAngles(jj)+obj.hingeAngles(jj)))*...
                             (obj.hingeAngles(jj) + halfAngles(jj) - pi/2))^2);
                     end
                     C = 0.364 + 0.28 * Gstar/Gc;
                     cornerArea(jj) = cornerArea(jj) * obj.b(jj)^2;
                     cornerConductance(jj) =  C * cornerArea(jj)^2 * Gc /obj.waterViscosity;                      
                     end
                 end
                 
                 obj.waterConductance = sum(cornerConductance);
                 obj.waterCrossSectionArea = sum(cornerArea);              
                 obj.waterSaturation = min(obj.waterCrossSectionArea / obj.area, 1);
                 obj.oilSaturation = max(1-obj.waterSaturation,0);
                 
                 if strcmp(obj.geometry , 'Triangle')== 1
                     obj.oilConductance = obj.area^2 * 3  * obj.shapeFactor / obj.oilViscosity/5 *obj.oilSaturation;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     obj.oilConductance = obj.area^2 *0.5623 * obj.shapeFactor /obj.oilViscosity*obj.oilSaturation;
                 end    
             end             
            else
                 obj.waterCrossSectionArea = obj.area;
                 obj.waterConductance = obj.conductanceSinglePhase; 
                 obj.oilCrossSectionArea = 0;
                 obj.oilConductance = 0; 
                 obj.waterSaturation = 1;
                 obj.oilSaturation = 0;
            end
            
           else
                 if strcmp(obj.geometry , 'Triangle')== 1
                     nc = 3;
                 else
                     nc = 4;
                 end
                 
                 R = network.sig_ow / Pc; %Raduis of Curvature
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 inerArea = zeros(1,nc);
                 outerArea = zeros(1,nc);
                 layerArea = zeros(1,nc);
                 inerConductance = zeros(1,nc);
                 layerConductance = zeros(1,nc);
                 for jj = 1:nc
                     if ~isnan(obj.oilLayerExist(1,jj)) %if the layer exist in the corner
                         
                         R_min = network.sig_ow / network.Pc_drain_max;
                         part = R_min * cos(obj.recedingContactAngle+halfAngles(jj))/R;
                         if part > 1
                             part = 1;
                         elseif part < -1
                             part = -1;
                         end
                         thetaHingRec = acos( part)-halfAngles(jj); 
                         thetaHingRec = min (thetaHingRec, obj.advancingContactAngle);                          
                         
                         % Area of corner water
                         if (halfAngles(jj) + thetaHingRec) == pi/2
                             inerArea(jj) = (R * cos(thetaHingRec + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             inerArea(jj) = R ^2 * (cos(thetaHingRec)*...
                                 (cot(halfAngles(jj)) * cos(thetaHingRec) - sin(thetaHingRec))+ ...
                                 thetaHingRec+ halfAngles(jj) - pi/2);
                         end
                         % Conductance of corner water
                         if thetaHingRec <= pi/2 - halfAngles(jj) % Positive Curvature       
                             f = 1; % no-flow boundary condition suitable for oil-water interfaces
                             F1 = pi/2 - halfAngles(jj) - thetaHingRec;
                             F2 = cot(halfAngles(jj)) * cos(thetaHingRec) - sin(thetaHingRec); 
                             F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                             
                             inerConductance(jj) = (inerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                                 (F2 * cos(thetaHingRec) - F1) * F3 ^ 2) / ...
                                 (12 * network.waterViscosity * ((sin(halfAngles(jj))) * ...
                                 (1 - F3) * (F2 + f * F1))^ 2);
                         elseif (thetaHingRec > pi/2 - halfAngles(jj)) % Negative Curvature
                             inerConductance(jj) = (inerArea(jj)^2 * tan(halfAngles(jj))* ...
                                 (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                                 (12 * network.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
                         end
                         thetaHingAdv = pi - obj.advancingContactAngle;
                         if (halfAngles(jj) + thetaHingAdv) == pi/2
                             outerArea(jj) = (R * cos(thetaHingAdv + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             outerArea(jj) = R ^2 * (cos(thetaHingAdv)*...
                                 (cot(halfAngles(jj)) * cos(thetaHingAdv) - sin(thetaHingAdv))+ ...
                                 thetaHingAdv+ halfAngles(jj) - pi/2);
                         end
                         % Area of oil layer
                         layerArea(jj) = outerArea(jj)-inerArea(jj);
                         % Conductance of oil layer
                         F3_layer = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                         layerConductance(jj) = (layerArea(jj)^3 * (1 - sin(halfAngles(jj)))^2 * ...
                             tan(halfAngles(jj)) * F3_layer ^ 2) / ...
                             (12 * network.oilViscosity *outerArea(jj)* (sin(halfAngles(jj)))^2 * ...
                             (1 - F3_layer) * (1 + F3_layer - (1- F3) * sqrt(inerArea(jj)/outerArea(jj))));
                     end
                 end
                 % Center water area and conductance
                 centerWaterArea = obj.area - sum(outerArea);
                 if strcmp(obj.geometry , 'Triangle')== 1
                     centerWaterConductance = 3 *obj.radius^2*centerWaterArea/20/network.waterViscosity;
                 else                     
                     centerWaterConductance = 0.5623 * obj.shapeFactor * centerWaterArea^2/network.waterViscosity;
                 end
                 
                 obj.waterCrossSectionArea = obj.area - sum(layerArea);
                 obj.waterConductance = centerWaterConductance + sum(inerConductance);
                 obj.oilCrossSectionArea = sum(layerArea);
                 obj.oilConductance = sum(layerConductance);
                 obj.waterSaturation = obj.waterCrossSectionArea / obj.area;
                 obj.oilSaturation = obj.oilCrossSectionArea / obj.area;
           end            
            if obj.waterSaturation > 1 
                % Control 
                fprintf('sat %f %4.0d  %f \n',Pc, obj.index,obj.waterSaturation);  
                obj.control(1:3) = [Pc, obj.imbThresholdPressure_PistonLike, obj.imbThresholdPressure_SnapOff];
                
                     obj.waterCrossSectionArea = obj.area;
                     obj.waterConductance = obj.conductanceSinglePhase; 
                     obj.oilCrossSectionArea = 0;
                     obj.oilConductance = 0; 
                     obj.waterSaturation = 1;
                     obj.oilSaturation = 0;
            end
        end 
        % Based on Patzek & Piri eqs. with hinging angle
        function calculateConductance_Imbibition_KH(obj, network, Pc)
             Pc = abs(Pc); 
             if ~any(obj.oilLayerExist) 
                 
                 if obj.occupancy == 'B' 
                     
                     if strcmp(obj.geometry , 'Circle')== 1 
                         
                         obj.waterCrossSectionArea = 0;
                         obj.waterConductance = 0; 
                         obj.oilCrossSectionArea = obj.area;
                         obj.oilConductance = obj.conductanceSinglePhase; 
                         obj.waterSaturation = 0;
                         obj.oilSaturation = 1;
                     else
                         
                         halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                         cornerArea = zeros(1,4);    cornerConductance = zeros(1,4);
                         % Based on  Al-Futaisi&Patzek_2003: eqs 12-14             
                         for i = 1:4     
                             if ~isnan(halfAngles(i)) && ~isnan(obj.waterCornerExist(i))  
                                 %Raduis of Curvature
                                 rso1 = obj.sig_ow / Pc;  
                                 rpd = obj.sig_ow / network.Pc_drain_max;  
%                                  r_adv = abs(rpd * cos(obj.recedingContactAngle+halfAngles(i))/cos(obj.advancingContactAngle+halfAngles(i)));
%                                  if rso1 > r_adv
%                                      rso1 = r_adv;
%                                  end
                                 part = rpd * cos(obj.recedingContactAngle+halfAngles(i))/rso1; 
                                 obj.hingeAngles(i) = acos(part) - halfAngles(i); 
                                 obj.hingeAngles(i) = min (obj.hingeAngles(i), obj.advancingContactAngle); 
                                 obj.b(i) = rso1*cos(obj.hingeAngles(i) + halfAngles(i))/sin(halfAngles(i)); 
                                 if obj.b(i) > obj.radius / tan(halfAngles(i))
                                     obj.b(i) = obj.radius / tan(halfAngles(i));
                                 end
                                 % Area
                                 % Based on Piri_2005: eq A4 & A5 
                                 if (obj.hingeAngles(i) + obj.recedingContactAngle) == pi/2
                                     cornerArea(i) = (rso1 * cos(obj.hingeAngles(i) + halfAngles(i))/ sin (halfAngles(i)))^2 *...
                                         sin(halfAngles(i)) * cos(halfAngles(i));
                                 else
                                     cornerArea(i) = rso1 ^2 * (cos(obj.hingeAngles(i))*...                  
                                         (cot(halfAngles(i)) * cos(obj.hingeAngles(i)) - sin(obj.hingeAngles(i)))+ ...
                                         obj.hingeAngles(i) + halfAngles(i) - pi/2);
                                 end
                                 %Conductance
                                 % Based on Piri_2005: eq B(10 - 15)
                                 f = 1; % no-flow boundary condition suitable for oil-water interfaces
                                 F1 = pi/2 - halfAngles(i) - obj.hingeAngles(i);
                                 F2 = cot(halfAngles(i)) * cos(obj.hingeAngles(i)) - sin(obj.hingeAngles(i)); 
                                 F3 = (pi/2 - halfAngles(i)) * tan(halfAngles(i));
                                 
                                 if (obj.hingeAngles(i) <= pi/2 - halfAngles(i))  % Positive Curvature                       
                                     cornerConductance(i) = (cornerArea(i)^2 * (1 - sin(halfAngles(i)))^2 * ...
                                         (F2 * cos(obj.hingeAngles(i)) - F1) * F3 ^ 2) / ...
                                         (12 * obj.waterViscosity * ((sin(halfAngles(i))) * ...
                                         (1 - F3) * (F2 + f * F1))^ 2);
                                 elseif (obj.hingeAngles(i) > pi/2 - halfAngles(i)) % Negative Curvature
                                     cornerConductance(i) = (cornerArea(i)^2 * tan(halfAngles(i))* ...
                                         (1 - sin(halfAngles(i)))^2 * F3 ^ 2) / ...
                                         (12 * obj.waterViscosity *(sin(halfAngles(i)))^2*(1 - F3) * (1 + f * F3)^ 2);
                                 end     
                             end
                             
                             obj.waterCrossSectionArea = sum(cornerArea);
                             obj.waterConductance = sum(cornerConductance);
                         end
                         obj.oilCrossSectionArea = obj.area - obj.waterCrossSectionArea;
                         if strcmp(obj.geometry , 'Triangle')== 1
                             obj.oilConductance = obj.oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                         elseif strcmp(obj.geometry , 'Square')== 1
                             obj.oilConductance = obj.oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                         end
                         obj.waterSaturation = obj.waterCrossSectionArea/obj.area;                         
                         obj.oilSaturation = obj.oilCrossSectionArea/obj.area;
                     end                     
                 else 
                     obj.waterCrossSectionArea = obj.area;
                     obj.waterConductance = obj.conductanceSinglePhase; 
                     obj.oilCrossSectionArea = 0;
                     obj.oilConductance = 0; 
                     obj.waterSaturation = 1;
                     obj.oilSaturation = 0;
                 end
                 
             else
                 if strcmp(obj.geometry , 'Triangle')== 1
                     nc = 3;
                 else
                     nc = 4;
                 end
                 
                 R = network.sig_ow / Pc; %Raduis of Curvature
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 inerArea = zeros(1,nc);
                 outerArea = zeros(1,nc);
                 layerArea = zeros(1,nc);
                 inerConductance = zeros(1,nc);
                 layerConductance = zeros(1,nc);
                 for jj = 1:nc
                     if ~isnan(obj.oilLayerExist(1,jj)) %if the layer exist in the corner
                         
                         R_min = network.sig_ow / network.Pc_drain_max;                         
                         part = R_min * cos(obj.recedingContactAngle+halfAngles(jj))/R; 
                         thetaHingRec = acos( part)-halfAngles(jj);                          
                         thetaHingRec = min(thetaHingRec, obj.advancingContactAngle); 
                         
                         % Area of corner water
                         if (halfAngles(jj) + thetaHingRec) == pi/2
                             inerArea(jj) = (R * cos(thetaHingRec + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             inerArea(jj) = R ^2 * (cos(thetaHingRec)*...
                                 (cot(halfAngles(jj)) * cos(thetaHingRec) - sin(thetaHingRec))+ ...
                                 thetaHingRec+ halfAngles(jj) - pi/2);
                         end                         
                               
                             f = 1; % no-flow boundary condition suitable for oil-water interfaces
                             F1 = pi/2 - halfAngles(jj) - thetaHingRec;
                             F2 = cot(halfAngles(jj)) * cos(thetaHingRec) - sin(thetaHingRec); 
                             F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                             
                         % Conductance of corner water
                         if thetaHingRec <= pi/2 - halfAngles(jj) % Positive Curvature 
                             
                             inerConductance(jj) = (inerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                                 (F2 * cos(thetaHingRec) - F1) * F3 ^ 2) / ...
                                 (12 * network.waterViscosity * ((sin(halfAngles(jj))) * ...
                                 (1 - F3) * (F2 + f * F1))^ 2);
                         elseif (thetaHingRec > pi/2 - halfAngles(jj)) % Negative Curvature
                             inerConductance(jj) = (inerArea(jj)^2 * tan(halfAngles(jj))* ...
                                 (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                                 (12 * network.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
                         end
                         
                         thetaHingAdv = pi - obj.advancingContactAngle;
                         if (halfAngles(jj) + thetaHingAdv) == pi/2
                             outerArea(jj) = (R * cos(thetaHingAdv + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             outerArea(jj) = R ^2 * (cos(thetaHingAdv)*...
                                 (cot(halfAngles(jj)) * cos(thetaHingAdv) - sin(thetaHingAdv))+ ...
                                 thetaHingAdv+ halfAngles(jj) - pi/2);
                         end
                         % Area of oil layer
                         layerArea(jj) = outerArea(jj)-inerArea(jj);
                         % Conductance of oil layer
                         F3_layer = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                         layerConductance(jj) = (layerArea(jj)^3 * (1 - sin(halfAngles(jj)))^2 * ...
                             tan(halfAngles(jj)) * F3_layer ^ 2) / ...
                             (12 * network.oilViscosity *outerArea(jj)* (sin(halfAngles(jj)))^2 * ...
                             (1 - F3_layer) * (1 + F3_layer - (1- F3) * sqrt(inerArea(jj)/outerArea(jj))));
                     end
                 end
                 % Center water area and conductance
                 centerWaterArea = obj.area - sum(outerArea);
                 if strcmp(obj.geometry , 'Triangle')== 1
                     centerWaterConductance = 3 *obj.radius^2*centerWaterArea/20/network.waterViscosity;
                 else                     
                     centerWaterConductance = 0.5623 * obj.shapeFactor * centerWaterArea^2/network.waterViscosity;
                 end
                 
                 obj.waterCrossSectionArea = obj.area - sum(layerArea);
                 obj.waterConductance = centerWaterConductance + sum(inerConductance);
                 obj.oilCrossSectionArea = sum(layerArea);
                 obj.oilConductance = sum(layerConductance);
                 obj.waterSaturation = obj.waterCrossSectionArea/obj.area;                 
                 obj.oilSaturation = obj.oilCrossSectionArea/obj.area;
             end
            
            if obj.waterSaturation > 1
                % Control 
                fprintf('sat %f %4.0d  %f \n',Pc, obj.index,obj.waterSaturation);  
                obj.control(1:3) = [Pc, obj.imbThresholdPressure_PistonLike, obj.imbThresholdPressure_SnapOff];
                
                     obj.waterCrossSectionArea = obj.area;
                     obj.waterConductance = obj.conductanceSinglePhase; 
                     obj.oilCrossSectionArea = 0;
                     obj.oilConductance = 0; 
                     obj.waterSaturation = 1;
                     obj.oilSaturation = 0;
            end
         end         
        % Based on Zolfaghari
        function calculateConductance_Imbibition_Z(obj, network, Pc)
             Pc = abs(Pc); 
             
             if ~any(obj.oilLayerExist) 
                 
                 if obj.occupancy == 'B'  
                     if strcmp(obj.geometry , 'Circle')== 1 
                         
                         obj.waterCrossSectionArea = 0;
                         obj.waterConductance = 0; 
                         obj.oilCrossSectionArea = obj.area;
                         obj.oilConductance = obj.conductanceSinglePhase; 
                         obj.waterSaturation = 0;
                         obj.oilSaturation = 1;
                     else 
                         halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                         cornerArea = zeros(1,4);    cornerConductance = zeros(1,4);
                         % Based on  Al-Futaisi&Patzek_2003: eqs 12-14             
                         for i = 1:4     
                             if ~isnan(halfAngles(i)) && ~isnan(obj.waterCornerExist(i))   
                                 %Raduis of Curvature
                                 rso1 = obj.sig_ow / Pc; 
                                 rpd = obj.sig_ow / network.Pc_drain_max; 
                                 part = rpd * cos(obj.recedingContactAngle+halfAngles(i))/rso1; 
                                 obj.hingeAngles(i) = acos(part) - halfAngles(i);                                 
                                 obj.hingeAngles(i) = min (obj.hingeAngles(i), obj.advancingContactAngle); 
                                 obj.b(i) = rso1 * cos(obj.hingeAngles(i) + halfAngles(i))/sin(halfAngles(i));  
                                 % Area
                                     cornerArea(i) = rso1^2 * (obj.hingeAngles(i) + halfAngles(i) - pi/2 ...
                                         +  cos(obj.hingeAngles(i)) * cos(obj.hingeAngles(i) + halfAngles(i))/ sin (halfAngles(i)));
                                 %Conductance
                                 % Based on Piri_2005: eq B(10 - 15)
                                 f = 1; % no-flow boundary condition suitable for oil-water interfaces
                                 F1 = pi/2 - halfAngles(i) - obj.hingeAngles(i);
                                 F2 = cot(halfAngles(i)) * cos(obj.hingeAngles(i)) - sin(obj.hingeAngles(i)); 
                                 F3 = (pi/2 - halfAngles(i)) * tan(halfAngles(i));
                                 
                                 if (obj.hingeAngles(i) <= pi/2 - halfAngles(i))  % Positive Curvature                       
                                     cornerConductance(i) = (cornerArea(i)^2 * (1 - sin(halfAngles(i)))^2 * ...
                                         (F2 * cos(obj.hingeAngles(i)) - F1) * F3 ^ 2) / ...
                                         (12 * obj.waterViscosity * ((sin(halfAngles(i))) * ...
                                         (1 - F3) * (F2 + f * F1))^ 2);
                                 elseif (obj.hingeAngles(i) > pi/2 - halfAngles(i)) % Negative Curvature
                                     cornerConductance(i) = (cornerArea(i)^2 * tan(halfAngles(i))* ...
                                         (1 - sin(halfAngles(i)))^2 * F3 ^ 2) / ...
                                         (12 * obj.waterViscosity *(sin(halfAngles(i)))^2*(1 - F3) * (1 + f * F3)^ 2);
                                 end 
                             end
                             
                             obj.waterCrossSectionArea = sum(cornerArea);
                             obj.waterConductance = sum(cornerConductance);
                         end
                         obj.oilCrossSectionArea = obj.area - obj.waterCrossSectionArea;
                         if strcmp(obj.geometry , 'Triangle')== 1
                             obj.oilConductance = obj.oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                         elseif strcmp(obj.geometry , 'Square')== 1
                             obj.oilConductance = obj.oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                         end
                         obj.waterSaturation = obj.waterCrossSectionArea/obj.area;
                         obj.oilSaturation = obj.oilCrossSectionArea/obj.area;
                     end                     
                 else 
                     obj.waterCrossSectionArea = obj.area;
                     obj.waterConductance = obj.conductanceSinglePhase; 
                     obj.oilCrossSectionArea = 0;
                     obj.oilConductance = 0; 
                     obj.waterSaturation = 1;
                     obj.oilSaturation = 0;
                 end
                 
             else
                 if strcmp(obj.geometry , 'Triangle')== 1
                     nc = 3;
                 else
                     nc = 4;
                 end
                 
                 R = network.sig_ow / Pc; %Raduis of Curvature
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 inerArea = zeros(1,nc);
                 outerArea = zeros(1,nc);
                 layerArea = zeros(1,nc);
                 inerConductance = zeros(1,nc);
                 layerConductance = zeros(1,nc);
                 for jj = 1:nc
                     if ~isnan(obj.oilLayerExist(1,jj)) %if the layer exist in the corner
                         
                         R_min = network.sig_ow / network.Pc_drain_max;                    
                         part = R_min * cos(obj.recedingContactAngle+halfAngles(jj))/R;
                         if part > 1
                             part = 1;
                         elseif part < -1
                             part = -1;
                         end
                         thetaHingRec = acos( part)-halfAngles(jj);         
                         thetaHingRec = min(thetaHingRec, obj.advancingContactAngle); 
                         
                         % Area of corner water
                         if (halfAngles(jj) + thetaHingRec) == pi/2
                             inerArea(jj) = (R * cos(thetaHingRec + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             inerArea(jj) = R ^2 * (cos(thetaHingRec)*...
                                 (cot(halfAngles(jj)) * cos(thetaHingRec) - sin(thetaHingRec))+ ...
                                 thetaHingRec+ halfAngles(jj) - pi/2);
                         end                         
                               
                             f = 1; % no-flow boundary condition suitable for oil-water interfaces
                             F1 = pi/2 - halfAngles(jj) - thetaHingRec;
                             F2 = cot(halfAngles(jj)) * cos(thetaHingRec) - sin(thetaHingRec); 
                             F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                             
                         % Conductance of corner water
                         if thetaHingRec <= pi/2 - halfAngles(jj) % Positive Curvature 
                             
                             inerConductance(jj) = (inerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                                 (F2 * cos(thetaHingRec) - F1) * F3 ^ 2) / ...
                                 (12 * network.waterViscosity * ((sin(halfAngles(jj))) * ...
                                 (1 - F3) * (F2 + f * F1))^ 2);
                         elseif (thetaHingRec > pi/2 - halfAngles(jj)) % Negative Curvature
                             inerConductance(jj) = (inerArea(jj)^2 * tan(halfAngles(jj))* ...
                                 (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                                 (12 * network.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
                         end
                         
                         thetaHingAdv = pi - obj.advancingContactAngle;
                         if (halfAngles(jj) + thetaHingAdv) == pi/2
                             outerArea(jj) = (R * cos(thetaHingAdv + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             outerArea(jj) = R ^2 * (cos(thetaHingAdv)*...
                                 (cot(halfAngles(jj)) * cos(thetaHingAdv) - sin(thetaHingAdv))+ ...
                                 thetaHingAdv+ halfAngles(jj) - pi/2);
                         end
                         % Area of oil layer
                         layerArea(jj) = outerArea(jj)-inerArea(jj);
                         % Conductance of oil layer
                         F3_layer = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                         layerConductance(jj) = (layerArea(jj)^3 * (1 - sin(halfAngles(jj)))^2 * ...
                             tan(halfAngles(jj)) * F3_layer ^ 2) / ...
                             (12 * network.oilViscosity *outerArea(jj)* (sin(halfAngles(jj)))^2 * ...
                             (1 - F3_layer) * (1 + F3_layer - (1- F3) * sqrt(inerArea(jj)/outerArea(jj))));
                     end
                 end
                 % Center water area and conductance
                 centerWaterArea = obj.area - sum(outerArea);
                 if strcmp(obj.geometry , 'Triangle')== 1
                     centerWaterConductance = 3 *obj.radius^2*centerWaterArea/20/network.waterViscosity;
                 else                     
                     centerWaterConductance = 0.5623 * obj.shapeFactor * centerWaterArea^2/network.waterViscosity;
                 end
                 
                 obj.waterCrossSectionArea = obj.area - sum(layerArea);
                 obj.waterConductance = centerWaterConductance + sum(inerConductance);
                 obj.oilCrossSectionArea = sum(layerArea);
                 obj.oilConductance = sum(layerConductance);
                 obj.waterSaturation = obj.waterCrossSectionArea/obj.area;                 
                 obj.oilSaturation = obj.oilCrossSectionArea/obj.area;
             end
            if obj.waterSaturation > 1
                % Control 
                fprintf('sat %f %4.0d  %f \n',Pc, obj.index,obj.waterSaturation);  
                obj.control(1:3) = [Pc, obj.imbThresholdPressure_PistonLike, obj.imbThresholdPressure_SnapOff];
                
                     obj.waterCrossSectionArea = obj.area;
                     obj.waterConductance = obj.conductanceSinglePhase; 
                     obj.oilCrossSectionArea = 0;
                     obj.oilConductance = 0; 
                     obj.waterSaturation = 1;
                     obj.oilSaturation = 0;
            end
        end 
    end
     
end

