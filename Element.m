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
        
        % Element conductance 
        conductance        
        cylindricalConductance
        nodeLinkSystemConductance_W
        nodeLinkSystemConductance_O
        
        geometry % geometrical shape of the element
        halfAngle1
        halfAngle2
        halfAngle3
        halfAngle4
        area
        
        %It's better to define the phase conductances and the area of each
        %fluid as a structure for this each element        
        waterCrossSectionArea
        oilCrossSectionArea
        gasCrossSectionArea 
        
        waterConductance
        oilConductance 
        gasConductance
        
        concentration = 0; 
        
        hingeAngles 
        recedingContactAngle = 0;
        advancingContactAngle = 6*pi/18;  
        wettabilityAlteration = 0*pi/18;
        
        
        waterPressure
        oilPressure
        gasPressure
        thresholdPressure % Capillary Threshold Pressure in Drainage
        occupancy = 'A';  % Element filled by Water
        isInvaded = false; % check for invasion in Imbibition
        
        ThresholdPressure_PistonLike = nan;
        ThresholdPressure_SnapOff = nan;
        ThresholdPressure_LayerCollapse = nan(1,4);
        oilLayerExist = nan(1,4);
        b = zeros(1,4);
        
        % reactive properties        
        adsorbedConcentration = 0;
        adsorbedConcentration_SolidFluid = 0;
        adsorbedConcentration_FluidFluid = 0;
    end    
    
    methods        
        %% Drainage
        %% Piston-Like
        function ThresholdPressure = calculateThresholdPressurePistonLike (obj)
         % calculateThresholdPressurePistonLike Summary of this method goes here
         % Detailed explanation goes here  
         % Based ob Al-Futaisi & Patzek 2001: eqs 2-5
             if strcmp(obj.geometry , 'Circle')== 1
                 ThresholdPressure = 2*obj.sig_ow *cos(obj.recedingContactAngle)/obj.radius;
             else                 
                 nominator = 0;
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 for j = 1:4
                     if ~isnan(halfAngles(j)) &&  halfAngles(j) < pi/2 - obj.recedingContactAngle
                         E2 = cos(obj.recedingContactAngle + halfAngles(j)) *...
                             cos(obj.recedingContactAngle) / sin(halfAngles(j));
                         E0 = pi / 2 - obj.recedingContactAngle - halfAngles(j);
                         nominator = nominator +  (E2 - E0); 
                     end
                 end
                 ThresholdPressure = (obj.sig_ow / obj.radius)*...
                     cos(obj.recedingContactAngle)*(1+sqrt(1 -(4*obj.shapeFactor*...
                     nominator)/(cos(obj.recedingContactAngle)^2)));
             end
        end
        %% Conductance Calculation _ Drainage
         function [waterCrossSectionArea, waterConductance, oilCrossSectionArea, oilConductance] = ...
                 calculateConductance_Drainage(obj, Pc, Pc_max_drainage)  
             
             if obj.occupancy == 'B' %element was occupied by oil
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 cornerArea = zeros(1,4);
                 cornerConductance = zeros(1,4);
                 %Raduis of Curvature
                 Rc = obj.sig_ow / abs(Pc);
                 for jj = 1:4
                     if ~isnan(halfAngles(jj))                         
                         obj.b(jj) = obj.sig_ow / Pc_max_drainage * ...
                             cos(obj.recedingContactAngle + halfAngles(jj)) * sin(halfAngles(jj));  
                             if obj.b(jj) < 0
                                 obj.b(jj) = 0;
                             end
                         % Area
                         % Based on Piri_2005: eq A4 & A5 
                         if (halfAngles(jj) + obj.recedingContactAngle) == pi/2
                             cornerArea(jj) = (Rc * cos(obj.recedingContactAngle + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                                 sin(halfAngles(jj)) * cos(halfAngles(jj));
                         else
                             cornerArea(jj) = Rc ^2 * (cos(obj.recedingContactAngle)*...
                                 (cot(halfAngles(jj)) * cos(obj.recedingContactAngle) - sin(obj.recedingContactAngle))+ ...
                                 obj.recedingContactAngle + halfAngles(jj) - pi/2);
                         end
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
                     waterCrossSectionArea = sum(cornerArea);
                     if waterCrossSectionArea>obj.area
                         waterCrossSectionArea = obj.area;
                     end
                     waterConductance = sum(cornerConductance);
                 end
                 oilCrossSectionArea = obj.area - waterCrossSectionArea;
                 if strcmp(obj.geometry , 'Circle')== 1
                     oilConductance = oilCrossSectionArea^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                 elseif strcmp(obj.geometry , 'Triangle')== 1
                     oilConductance = oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     oilConductance = oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                 end
             else
                 waterCrossSectionArea = obj.area;
                 waterConductance = obj.conductance; 
                 oilCrossSectionArea = 0;
                 oilConductance = 0; 
             end
         end 
         function [waterCrossSectionArea, waterConductance, oilCrossSectionArea, oilConductance] = ...
                 calculateConductance_Drainage_Valvatne(obj, Pc_max_drainage)  
            
             if obj.occupancy == 'B' %element was occupied by oil
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 cornerArea = zeros(1,4);
                 cornerConductance = zeros(1,4);
                  
                 for jj = 1:4
                     if ~isnan(halfAngles(jj))
                         
                         obj.b(jj) = obj.sig_ow / Pc_max_drainage * ...
                             cos(obj.recedingContactAngle + halfAngles(jj)) * sin(halfAngles(jj));
                         if obj.b(jj) < 0
                                 obj.b(jj) = 0;
                         end
                         % Area
                         % Based on Valvatne 3.45-3.49                          
                         cornerArea(jj) = (obj.b(jj)*sin(halfAngles(jj))/cos(halfAngles(jj)+obj.recedingContactAngle)) ^2 ...
                                 * (cos(obj.recedingContactAngle)*(cos(halfAngles(jj)+obj.recedingContactAngle)/ sin(halfAngles(jj)))+ ...
                                 obj.recedingContactAngle + halfAngles(jj) - pi/2); 
                         %Conductance
                         Gc = cornerArea(jj)/...
                             4*obj.b(jj)^2/(1-sin(halfAngles(jj))/cos(halfAngles(jj)+obj.recedingContactAngle)/...
                             (obj.recedingContactAngle + halfAngles(jj) - pi/2))^2;
                         Gstar = sin(halfAngles(jj))*cos(halfAngles(jj))/4/(1+sin(halfAngles(jj)))^2;
                         C = 0.364 + 0.28 * Gstar/Gc;
                          cornerConductance(jj) =  C * cornerArea(jj)^2*Gc/obj.waterViscosity;
                     end
                     waterCrossSectionArea = sum(cornerArea);
                     if waterCrossSectionArea>obj.area
                         waterCrossSectionArea = obj.area;
                     end
                     waterConductance = sum(cornerConductance);
                 end
                 oilCrossSectionArea = obj.area - waterCrossSectionArea;
                 if strcmp(obj.geometry , 'Circle')== 1
                     oilConductance = oilCrossSectionArea^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                 elseif strcmp(obj.geometry , 'Triangle')== 1
                     oilConductance = oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     oilConductance = oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                 end
             else
                 waterCrossSectionArea = obj.area;
                 waterConductance = obj.conductance; 
                 oilCrossSectionArea = 0;
                 oilConductance = 0; 
             end
         end
        %% Imbibition
        %% Piston-Like & oilLayerExist
        % Based on Patzek 2003
        function [ThresholdPressure_PistonLike, oilLayerExist] = calculateThresholdPressurePistonLike_Imbibition_2003(obj, Pc_max_drainage)
         % calculateThresholdPressurePistonLike Summary of this method goes here
         % Detailed explanation goes here  
         
         if strcmp(obj.geometry , 'Circle')== 1
                 % Based on Al-Futaisi&Patzek_2001: eqs 2-5 & Piri_2005: eq C4 
                 ThresholdPressure_PistonLike = 2*obj.sig_ow *cos(obj.advancingContactAngle)/obj.radius;
                 oilLayerExist(1,:) = nan;
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
                             hingingAngles(ii) = min(obj.hingeAngles(ii) , obj.advancingContactAngle);
                             nominator1 = nominator1 + b_i(1,ii)*cos(hingingAngles(ii));
                             nominator2 = nominator2 + (pi/2 - hingingAngles(ii) - halfAngles(ii));
                         end
                     end
                     rp2 = ((obj.radius ^ 2 / 4 / obj.shapeFactor) -rp1*nominator1 + rp1^2 *nominator2) /...
                         (2*rp1 * sum(alpha) + cos(obj.advancingContactAngle)*...
                         ((obj.radius/2/obj.shapeFactor) - 2 * sum(b_i)));
                 end
                 ThresholdPressure_PistonLike = obj.sig_ow / rp2;
                 oilLayerExist(1,:) = nan;
             elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi/2 + max(halfAngles) %Forced Imbibition
                 ThresholdPressure_PistonLike = 2 * obj.sig_ow * cos(obj.advancingContactAngle) / obj.radius;
                 oilLayerExist(1,:) = nan;
             elseif obj.advancingContactAngle >= pi/2 + max(halfAngles)
                 ThresholdPressure_PistonLike = -calculateThresholdPressurePistonLike(obj, (pi - obj.advancingContactAngle)); %Forced Imbibition
                 oilLayerExist(1,:) = 1;
             end
         end
        end
        % Based on Patzek 2001
        function [ThresholdPressure_PistonLike, oilLayerExist] = calculateThresholdPressurePistonLike_Imbibition(obj, Pc_max_drainage)
       
         if strcmp(obj.geometry , 'Circle')== 1
                 % Based on Al-Futaisi&Patzek_2001: eqs 2-5 & Piri_2005: eq C4 
                 ThresholdPressure_PistonLike = 2*obj.sig_ow *cos(obj.advancingContactAngle)/obj.radius;
                 oilLayerExist(1,:) = nan;
         else
             % Based on  Al-Futaisi&Patzek_2001: eqs 2-4 & 6-10             
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];  
             b_i = zeros(1,4);
             if strcmp(obj.geometry , 'Triangle')== 1
                 nc = 3;
             else
                 nc = 4;
             end
             E0 = zeros(nc,1); E1 = zeros(nc,1); alpha = zeros(nc,1);hingingAngles = zeros(1,nc);
             nominator = 0;
             for i = 1:nc
                 if ~isnan(halfAngles(i))
                     nominator = nominator + cos(obj.recedingContactAngle + halfAngles(i));
                 end
             end
             maxAdvAngle = acos ((-4 * obj.shapeFactor * nominator)/...
                 ((obj.radius * Pc_max_drainage / obj.sig_ow) - cos(obj.recedingContactAngle)+...
                 4 * nc * obj.shapeFactor * sin(obj.recedingContactAngle)));             
             rpd = obj.sig_ow / Pc_max_drainage;             
             r_initialGuess = rpd;
             
             if obj.advancingContactAngle <= maxAdvAngle % Spontaneous imbibition
                 rp1 = 2 * r_initialGuess;
                 rp2 = r_initialGuess;
                 while abs(rp2 - rp1) > 10^-10 % fixed point iteration method
                     rp1 = rp2;                                         
                     for ii = 1:nc
                         hingingAngles(ii) = acos((rpd / rp1)*cos(obj.recedingContactAngle + halfAngles(ii))) - halfAngles(ii);
                         if ~isnan(halfAngles(ii)) && hingingAngles(ii) <= obj.advancingContactAngle
                             b_i(ii) = rpd * cos(obj.recedingContactAngle + halfAngles(ii))/ sin(halfAngles(ii));
                             alpha(ii) = asin(obj.b(ii)/rp1*sin(halfAngles(ii)));
                         elseif ~isnan(halfAngles(ii)) && hingingAngles(ii) > obj.advancingContactAngle                             
                             b_i(ii) = rp1 * cos(obj.advancingContactAngle + halfAngles(ii))/ sin(halfAngles(ii));
                             alpha(ii) = pi/2 - obj.advancingContactAngle - halfAngles(ii);
                         end 
                         if b_i(ii) < 0
                             b_i(ii) = 0;
                         end
                         hingingAngles(ii) = min ( hingingAngles(ii) , obj.advancingContactAngle);
                         E0(ii) = pi/2 - hingingAngles(ii) - halfAngles(ii);
                         E1(ii) = obj.b(ii) * cos(hingingAngles(ii)); 
                     end
                     rp2 = (obj.radius ^ 2 / 4 / obj.shapeFactor - rp1 * sum(E1) + rp1^2 * sum(E0))/...
                         (2*rp1 * sum(alpha) + (obj.radius/2/obj.shapeFactor - 2 * sum(obj.b)) * cos(obj.advancingContactAngle));
                 end 
                 ThresholdPressure_PistonLike = obj.sig_ow / rp2;
                 oilLayerExist(1,:) = nan;
             elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi/2 + max(halfAngles) % Forced imbibition
                 ThresholdPressure_PistonLike = 2 * obj.sig_ow * cos(obj.advancingContactAngle) / obj.radius;
                 oilLayerExist(1,:) = nan;
             elseif obj.advancingContactAngle >= pi/2 + max(halfAngles) % Forced imbibition
                 ThresholdPressure_PistonLike = -calculateThresholdPressurePistonLike(obj, (pi - obj.advancingContactAngle));
                 oilLayerExist(1,:) = 1;
             end
         end
        end 
        %% Snap-off        
        function ThresholdPressure_SnapOff = calculateThresholdPressureSnapOff_P(obj,Pc_max_drainage)
            
         if strcmp(obj.geometry , 'Circle')== 1
             ThresholdPressure_SnapOff = nan;
         else
             % Based on  Al-Futaisi&Patzek_2003: eqs 12-14          
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];             
             maxAdvAngle = pi/2 - min(halfAngles);
             if strcmp(obj.geometry , 'Triangle')== 1
                 nc = 3; 
             if obj.advancingContactAngle < maxAdvAngle %Spontaneous Imbibition
                 
                 rso = zeros(1,nc);
                 for edge = 1:nc   
                     
                     rso1 = obj.sig_ow / obj.ThresholdPressure_PistonLike; 
                     rso2 = rso1*2;
                     
                     while abs(rso1 - rso2) > 10^-12
                         rso1 = rso2;
                         hingeAngle_ii = acos((obj.sig_ow / rso1)*...
                             cos(obj.recedingContactAngle + halfAngles(1))/Pc_max_drainage) - halfAngles(1); 
                         if hingeAngle_ii <= obj.advancingContactAngle
                             E1_i = cos(obj.recedingContactAngle + halfAngles(1))/sin(halfAngles(1));
                         else
                             E1_i = cos(obj.advancingContactAngle + halfAngles(1))/sin(halfAngles(1));
                         end                             
                         
                         hingeAngle_jj = acos((obj.sig_ow / rso1)*...
                             cos(obj.recedingContactAngle + halfAngles(jj)) / Pc_max_drainage) - halfAngles(jj); 
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
                         rso = min(rso12,rso23,rso31);
                     end
                 end
                 ThresholdPressure_SnapOff = obj.sig_ow / rso;
                 
                 % Forced imbibition part
             elseif obj.advancingContactAngle == maxAdvAngle
                 ThresholdPressure_SnapOff = 0;
             elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
                 ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + min(halfAngles))/...
                     cos(obj.recedingContactAngle + min(halfAngles));
             elseif obj.advancingContactAngle >= pi - min(halfAngles)
                 ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
             end
             else 
                 % elemnt is square : Piri
                 if obj.advancingContactAngle <= pi/4
                     %eq C34
                     ThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cot(pi/4)*cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(obj.recedingContactAngle + pi/4);
                 elseif obj.advancingContactAngle > 3*pi/4
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end  
             end
         end         
        end
        % Piri for square elements
        function ThresholdPressure_SnapOff = calculateThresholdPressureSnapOff_Piri(obj,Pc_max_drainage)
            
         if strcmp(obj.geometry , 'Circle')== 1
             ThresholdPressure_SnapOff = nan;
         else
             % Based on  Al-Futaisi&Patzek_2003: eqs 12-14          
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];             
             maxAdvAngle = pi/2 - min(halfAngles);
             if strcmp(obj.geometry , 'Triangle')== 1
                 nc = 3;
             
                 if obj.advancingContactAngle < maxAdvAngle %Spontaneous Imbibition

                     rso = zeros(1,nc);
                     for edge = 1:nc         

                         ii = edge;                         
                         if ii == nc
                             jj = 1;                             
                         else
                             jj = edge + 1;
                         end

                         rso1 = obj.sig_ow / Pc_max_drainage; 
                         rso2 = rso1*2;

                         while abs(rso1 - rso2) > 10^-10
                             rso1 = rso2;
                             hingeAngle_ii = acos((obj.sig_ow / rso1)*...
                                 cos(obj.recedingContactAngle + halfAngles(ii))/Pc_max_drainage) - halfAngles(ii);                         
                             hingeAngle_jj = acos((obj.sig_ow / rso1)*...
                                 cos(obj.recedingContactAngle + halfAngles(jj)) / Pc_max_drainage) - halfAngles(jj);

                             if hingeAngle_ii <= obj.advancingContactAngle
                                 E1_i = cos(hingeAngle_ii + halfAngles(jj))/sin(halfAngles(jj));
                                 if hingeAngle_jj > obj.advancingContactAngle
                                     E1_j = cos(obj.advancingContactAngle + halfAngles(jj)) / sin(halfAngles(jj));

                                     rso2 = obj.radius *(cot(halfAngles(ii))+ cot(halfAngles(jj)))/(E1_i+E1_j);
                                end
                             else
                                 E1_i = cos(obj.advancingContactAngle + halfAngles(ii))/sin(halfAngles(ii));
                                 if hingeAngle_jj > obj.advancingContactAngle
                                     E1_j = cos(obj.advancingContactAngle + halfAngles(jj)) / sin(halfAngles(jj));

                                     rso2 = obj.radius *(cot(halfAngles(ii))+ cot(halfAngles(jj)))/(E1_i+E1_j);
                                 else
                                     E1_j = cos(hingeAngle_jj + halfAngles(jj))/sin(halfAngles(jj));       

                                     rso2 = obj.radius *(cot(halfAngles(ii))+ cot(halfAngles(jj)))/(E1_i+E1_j);
                                end
                             end                    
                             rso(edge) = rso2;
                         end
                     end
                     ThresholdPressure_SnapOff = obj.sig_ow / min(rso);

                     % Forced imbibition part
                 elseif obj.advancingContactAngle == maxAdvAngle
                     ThresholdPressure_SnapOff = 0;
                 elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + min(halfAngles))/...
                         cos(obj.recedingContactAngle + min(halfAngles));
                 elseif obj.advancingContactAngle >= pi - min(halfAngles)
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end
             
             else % elemnt is square
                 if obj.advancingContactAngle <= pi/4
                     %eq C34
                     ThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cot(pi/4)*cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(obj.recedingContactAngle + pi/4);
                 elseif obj.advancingContactAngle > 3*pi/4
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end                 
             end
         end         
        end
        % Boujelben
        function ThresholdPressure_SnapOff = calculateThresholdPressureSnapOff_Boujelben (obj,Pc_max_drainage)
            
         if strcmp(obj.geometry , 'Circle')== 1
             ThresholdPressure_SnapOff = nan;
         else
             % Based on  Al-Futaisi&Patzek_2003: eqs 12-14          
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];             
             maxAdvAngle = pi/2 - min(halfAngles);
             if strcmp(obj.geometry , 'Triangle')== 1
                 if obj.advancingContactAngle < maxAdvAngle %Spontaneous Imbibition
                     ThresholdPressure_SnapOff = obj.sig_ow /obj.radius;
                     %Forced imbibition part
                 elseif obj.advancingContactAngle == maxAdvAngle
                     ThresholdPressure_SnapOff = 0;
                 elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + min(halfAngles))/...
                         cos(obj.recedingContactAngle + min(halfAngles));
                 elseif obj.advancingContactAngle >= pi - min(halfAngles)
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end
             else
                 % elemnt is square
                 if obj.advancingContactAngle <= pi/4
                     %eq C34
                     ThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cot(pi/4)*cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(obj.recedingContactAngle + pi/4);
                 elseif obj.advancingContactAngle > 3*pi/4
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end  
             end
         end         
        end
        % Valvatne        
        function ThresholdPressure_SnapOff = calculateThresholdPressureSnapOff(obj,Pc_max_drainage)
         if strcmp(obj.geometry , 'Circle')== 1
             ThresholdPressure_SnapOff = nan;
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
                         cos(obj.halfAngle3)*cot(obj.halfAngle3)-sin(obj.halfAngle3))/((cot(obj.halfAngle1)+cot(obj.halfAngle2)));                 
                     ThresholdPressure_SnapOff = max(Pc_a,Pc_b);                 
                     % Forced imbibition part
                 elseif obj.advancingContactAngle == maxAdvAngle
                     ThresholdPressure_SnapOff = 0;
                 elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + min(halfAngles))/...
                         cos(obj.recedingContactAngle + min(halfAngles));
                 elseif obj.advancingContactAngle >= pi - min(halfAngles)
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end
             else 
                 % elemnt is square
                 if obj.advancingContactAngle <= pi/4
                     %eq C34
                     ThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cot(pi/4)*cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(obj.recedingContactAngle + pi/4);
                 elseif obj.advancingContactAngle > 3*pi/4
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end  
             end
         end         
        end
        % Ryazanov       
        function ThresholdPressure_SnapOff = calculateThresholdPressureSnapOff_R(obj,Pc_max_drainage)
         if strcmp(obj.geometry , 'Circle')== 1
             ThresholdPressure_SnapOff = nan;
         else
             % Based on  Al-Futaisi&Patzek_2003: eqs 12-14          
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];             
             maxAdvAngle = pi/2 - min(halfAngles);
             if strcmp(obj.geometry , 'Triangle')== 1 
%                  L_dr = obj.sig_ow / Pc_max_drainage * cos(obj.recedingContactAngle + obj.halfAngle1)/sin(obj.halfAngle1);
                 if obj.advancingContactAngle < maxAdvAngle %Spontaneous Imbibition
                     ThresholdPressure_SnapOff = obj.radius * sin(pi/3)/cos(obj.advancingContactAngle + obj.halfAngle1); 
                     % Forced imbibition part
                 elseif obj.advancingContactAngle == maxAdvAngle
                     ThresholdPressure_SnapOff = 0;
                 elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + min(halfAngles))/...
                         cos(obj.recedingContactAngle + min(halfAngles));
                 elseif obj.advancingContactAngle >= pi - min(halfAngles)
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end
%                      % Forced imbibition part
%                  elseif obj.advancingContactAngle == maxAdvAngle
%                      ThresholdPressure_SnapOff = 0;
%                  elseif obj.advancingContactAngle > maxAdvAngle && obj.advancingContactAngle < pi - min(halfAngles)
%                      ThresholdPressure_SnapOff = L_dr * sin(obj.halfAngle1)/ cos(obj.advancingContactAngle + min(halfAngles));
%                  elseif obj.advancingContactAngle >= pi - min(halfAngles)
%                      ThresholdPressure_SnapOff = -L_dr* sin(obj.halfAngle1);
%                  end
             else 
                 % elemnt is square
                 if obj.advancingContactAngle <= pi/4
                     %eq C34
                     ThresholdPressure_SnapOff = obj.sig_ow / obj.radius * ...
                         (cot(pi/4)*cos(obj.advancingContactAngle)-sin(obj.advancingContactAngle));
                 elseif obj.advancingContactAngle > pi/4 && obj.advancingContactAngle <= 3*pi/4
                     ThresholdPressure_SnapOff = Pc_max_drainage*cos(obj.advancingContactAngle + pi/4)/...
                         cos(obj.recedingContactAngle + pi/4);
                 elseif obj.advancingContactAngle > 3*pi/4
                     ThresholdPressure_SnapOff = -Pc_max_drainage/cos(obj.recedingContactAngle + min(halfAngles));
                 end  
             end
         end         
        end
        %% LayerCollapse Piri
        % Ver_1 based on Piri eqs.
        function ThresholdPressure_LayerCollapse = calculateThresholdPressureLayerCollapse (obj, Pc_max_drainage)
          Pc = Pc_max_drainage;  
         if strcmp(obj.geometry , 'Circle')== 1
             ThresholdPressure_LayerCollapse = nan(1,4);
         else
              halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4]; 
              if strcmp(obj.geometry , 'Triangle')== 1
                  nc = 3;
              elseif strcmp(obj.geometry , 'Square')== 1
                  nc = 4;
              end
              ThresholdPressure_LayerCollapse = nan(1,4); hingingAngles = zeros(1,4); b_i = zeros(1,4);
              for i = 1:nc
                  if obj.advancingContactAngle > pi/2 + halfAngles(i) % convex oil layer is impossible
                      
                  hingingAngles(i) = acos(Pc/Pc_max_drainage*cos(obj.recedingContactAngle + halfAngles(i))) - halfAngles(i);
                  b_i(i) = obj.sig_ow * cos(obj.recedingContactAngle + halfAngles(i))/sin(halfAngles(i))/Pc_max_drainage;
                  ThresholdPressure_LayerCollapse(i) = obj.sig_ow *...
                      (3*(sin(halfAngles(i)))^2+ 4*sin(halfAngles(i))*cos(obj.advancingContactAngle)+(cos(obj.advancingContactAngle))^2)/...
                      (b_i(i)*(cos(halfAngles(i))*sin(halfAngles(i))*(2*sin(halfAngles(i))+cos(obj.advancingContactAngle))+...
                      (sin(halfAngles(i)))^2 * ...
                      sqrt(4*(cos(halfAngles(i)))^2-3-(cos(obj.advancingContactAngle))^2-4*sin(halfAngles(i))*cos(obj.advancingContactAngle))));
                  end
              end
         end
        end
        % New based on Piri
        function ThresholdPressure_LayerCollapse = calculateThresholdPressureLayerCollapse_2(obj, Pc_max_drainage)
         
         if strcmp(obj.geometry , 'Circle')== 1
             ThresholdPressure_LayerCollapse = nan(1,4);
         else
              halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4]; 
%               maxRadius = obj.sig_ow / Pc_max_drainage;              
              if strcmp(obj.geometry , 'Triangle')== 1
                  nc = 3;
              elseif strcmp(obj.geometry , 'Square')== 1
                  nc = 4;
              end
              ThresholdPressure_LayerCollapse = nan(1, 4); hingingAngles = zeros(1,4);   b_i = zeros(nc , 1);
              for i = 1:nc
%                   Pc = Pc_max_drainage;
                  hingingAngles(i) = obj.advancingContactAngle;
%                   h = obj.hingeAngles/2;
%                   while abs(obj.hingeAngles - h) > 0.1
%                       if  obj.hingeAngles <= obj.advancingContactAngle
                          b_i(i) = obj.sig_ow / Pc_max_drainage * cos(obj.recedingContactAngle + halfAngles(i))/sin(halfAngles(i));
%                       else
%                           b_i(i) = obj.sig_ow / Pc * cos(obj.advancingContactAngle + halfAngles(i))/sin(halfAngles(i));
%                       end
                      Pc = obj.sig_ow *...
                          (3*(sin(halfAngles(i)))^2+ 4*sin(halfAngles(i))*cos(hingingAngles(i))+(cos(hingingAngles(i)))^2)/...
                          (b_i(i)*(cos(halfAngles(i))*sin(halfAngles(i))*(2*sin(halfAngles(i))+cos(hingingAngles(i)))+...
                          (sin(halfAngles(i)))^2 * ...
                          sqrt(4*(cos(halfAngles(i)))^2-3-(cos(hingingAngles(i)))^2-4*sin(halfAngles(i))*cos(hingingAngles(i)))));                       
%                       h = hingingAngles(i);
%                       hingingAngles(i) = acos((Pc/Pc_max_drainage)*cos(obj.advancingContactAngle + halfAngles(i))) - halfAngles(i);                      
%                   end
                  ThresholdPressure_LayerCollapse(i) = Pc;
              end
              
         end
        end
         %% Conductance Calculation _ Imbibition 
         % Based on Valvatne eqs.
         function [waterCrossSectionArea, waterConductance, oilCrossSectionArea, oilConductance] =...
                calculateConductance_Imbibition(obj, network, Pc)
             Pc = abs(Pc);
             
           if ~any(obj.oilLayerExist)
               
            if obj.occupancy == 'B' 
                
             if strcmp(obj.geometry , 'Circle')== 1
                  waterCrossSectionArea = 0;
                  waterConductance = 0;
                  oilCrossSectionArea = obj.area;
                  oilConductance = oilCrossSectionArea^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
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
                         
                         obj.hingeAngles(jj) = acos( R_min * ...
                             cos(obj.recedingContactAngle+halfAngles(jj))/R)-halfAngles(jj);
                         
                         if obj.hingeAngles(jj) > obj.advancingContactAngle
                             obj.hingeAngles(jj) = obj.advancingContactAngle;
                             bi = R * cos(obj.advancingContactAngle + halfAngles(jj))/...
                                 sin(halfAngles(jj));
                         else
                             bi = R_min * cos(obj.recedingContactAngle + halfAngles(jj))/...
                                 sin(halfAngles(jj));
                         end
                         if bi<0
                             bi = 0;
                         end
                         % Area of corner water
                         % Based on Valvatne
                         cornerArea(jj) = (bi * sin(halfAngles(jj))/cos(obj.hingeAngles(jj) + halfAngles(jj)))^2*...
                             (cos(obj.hingeAngles(jj)) * cos(obj.hingeAngles(jj)+ halfAngles(jj))/sin(halfAngles(jj)) +...
                             + obj.hingeAngles(jj) + halfAngles(jj)- pi/2);
                               
                             f = 1; % no-flow boundary condition suitable for oil-water interfaces
                             F1 = pi/2 - halfAngles(jj) - obj.hingeAngles(jj);
                             F2 = cot(halfAngles(jj)) * cos(obj.hingeAngles(jj)) - sin(obj.hingeAngles(jj)); 
                             F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));  
                             
                         % Conductance of corner water
                         if obj.hingeAngles(jj) <= pi/2 - halfAngles(jj) % Positive Curvature 
                             
                             cornerConductance(jj) = (cornerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                                 (F2 * cos(obj.hingeAngles(jj)) - F1) * F3 ^ 2) / ...
                                 (12 * network.waterViscosity * ((sin(halfAngles(jj))) * ...
                                 (1 - F3) * (F2 + f * F1))^ 2);
                         elseif (obj.hingeAngles(jj) > pi/2 - halfAngles(jj)) % Negative Curvature
                             cornerConductance(jj) = (cornerArea(jj)^2 * tan(halfAngles(jj))* ...
                                 (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                                 (12 * network.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
                         end
                 end
                 
                 waterConductance = sum(cornerConductance);
                 waterCrossSectionArea = sum(cornerArea);
                 if waterCrossSectionArea > obj.area && isnan(obj.ThresholdPressure_SnapOff)
                     waterCrossSectionArea = obj.area - R_min ^ 2 * pi;
                 elseif waterCrossSectionArea > obj.area
                     waterCrossSectionArea = obj.area;
                 end               
                 oilCrossSectionArea = obj.area - waterCrossSectionArea;
                 
                 if strcmp(obj.geometry , 'Triangle')== 1
                     oilConductance = oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                 elseif strcmp(obj.geometry , 'Square')== 1
                     oilConductance = oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                 end
             end
             else
                 waterCrossSectionArea = obj.area;
                 waterConductance = obj.conductance; 
                 oilCrossSectionArea = 0;
                 oilConductance = 0; 
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
                     if ~isnan(obj.ThresholdPressure_LayerCollapse(1,jj)) %if the layer exist in the corner
                         
                         R_min = network.sig_ow / network.Pc_drain_max;
                         thetaHingRec = acos( R_min * ...
                             cos(obj.recedingContactAngle+halfAngles(jj))/R)-halfAngles(jj);
                         
                         if thetaHingRec > obj.advancingContactAngle
                             thetaHingRec = obj.advancingContactAngle;
                         end                         
                         
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
                 if strcmp(obj.geometry , 'Circle')== 1
                     centerWaterConductance = 0.5 * obj.shapeFactor * centerWaterArea^2/network.waterViscosity;
                 elseif strcmp(obj.geometry , 'Triangle')== 1
                     centerWaterConductance = 3 *obj.radius^2*centerWaterArea/20/network.waterViscosity;
                 else                     
                     centerWaterConductance = 0.5623 * obj.shapeFactor * centerWaterArea^2/network.waterViscosity;
                 end
                 
                 waterCrossSectionArea = obj.area - sum(layerArea);
                 waterConductance = centerWaterConductance + sum(inerConductance);
                 oilCrossSectionArea = sum(layerArea);
                 oilConductance = sum(layerConductance);
            end
        end 
         % Based on Patzek eqs. with hinging angle
         function [waterCrossSectionArea, waterConductance, oilCrossSectionArea, oilConductance] =...
                calculateConductance_Imbibition_Patzek(obj, network, Pc)
             Pc = abs(Pc); 
             
             if ~any(obj.oilLayerExist) 
                 
                 if obj.occupancy == 'B' 
                     
                     if strcmp(obj.geometry , 'Circle')== 1 
                         
                         waterCrossSectionArea = 0;
                         waterConductance = 0; 
                         oilCrossSectionArea = obj.area;
                         oilConductance = obj.conductance; 
                     else
                         
                         halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                         cornerArea = zeros(1,4);    cornerConductance = zeros(1,4);
                         % Based on  Al-Futaisi&Patzek_2003: eqs 12-14             
                         for i = 1:4     
                             if ~isnan(halfAngles(i))    
                                 %Raduis of Curvature
                                 rso1 = obj.sig_ow / Pc;                   
                                 obj.hingeAngles(i) = acos((obj.sig_ow / rso1)*...
                                     cos(obj.recedingContactAngle + halfAngles(i))/Pc_max_drainage) - halfAngles(i); 
                                 if obj.hingeAngles(i) <= obj.advancingContactAngle
                                     obj.b(i) = cos(obj.recedingContactAngle + halfAngles(i))/sin(halfAngles(i));
                                 else
                                     obj.b(i) = cos(obj.advancingContactAngle + halfAngles(i))/sin(halfAngles(i));
                                 end
                                 if obj.b(i) < 0
                                     obj.b(i) = 0;
                                 end
                                 obj.hingeAngles(i) = min (obj.hingeAngles(i), obj.advancingContactAngle);
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
                             
                             waterCrossSectionArea = sum(cornerArea);
                             if waterCrossSectionArea>obj.area
                                 waterCrossSectionArea = obj.area;
                             end
                             waterConductance = sum(cornerConductance);
                         end
                         oilCrossSectionArea = obj.area - waterCrossSectionArea;
                         if strcmp(obj.geometry , 'Circle')== 1
                             oilConductance = oilCrossSectionArea^2 * 0.5 * obj.shapeFactor /obj.oilViscosity;
                         elseif strcmp(obj.geometry , 'Triangle')== 1
                             oilConductance = oilCrossSectionArea^2 * 3  * obj.shapeFactor / obj.oilViscosity/5;
                         elseif strcmp(obj.geometry , 'Square')== 1
                             oilConductance = oilCrossSectionArea^2 *0.5623 * obj.shapeFactor /obj.oilViscosity;
                         end
                     end                     
                 else 
                     waterCrossSectionArea = obj.area;
                     waterConductance = obj.conductance; 
                     oilCrossSectionArea = 0;
                     oilConductance = 0; 
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
                     if ~isnan(obj.ThresholdPressure_LayerCollapse(1,jj)) %if the layer exist in the corner
                         
                         R_min = network.sig_ow / network.Pc_drain_max;
                         thetaHingRec = acos( R_min * ...
                             cos(obj.recedingContactAngle+halfAngles(jj))/R)-halfAngles(jj);
                         
                         if thetaHingRec > obj.advancingContactAngle
                             thetaHingRec = obj.advancingContactAngle;
                         end                         
                         
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
                 if strcmp(obj.geometry , 'Circle')== 1
                     centerWaterConductance = 0.5 * obj.shapeFactor * centerWaterArea^2/network.waterViscosity;
                 elseif strcmp(obj.geometry , 'Triangle')== 1
                     centerWaterConductance = 3 *obj.radius^2*centerWaterArea/20/network.waterViscosity;
                 else                     
                     centerWaterConductance = 0.5623 * obj.shapeFactor * centerWaterArea^2/network.waterViscosity;
                 end
                 
                 waterCrossSectionArea = obj.area - sum(layerArea);
                 waterConductance = centerWaterConductance + sum(inerConductance);
                 oilCrossSectionArea = sum(layerArea);
                 oilConductance = sum(layerConductance);
            end
        end 
    end
end

