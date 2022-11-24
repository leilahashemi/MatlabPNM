function PLOTDRAIN_IMB(network)
%% Fluid info
fprintf('===============================Ntework Information=====================================\n');
fprintf('Water Viscosity:                        %f \n', network.waterViscosity); 
fprintf('Water Viscosity:                        %f \n', network.oilViscosity); 
fprintf('Interfacial Tension:                    %f \n', network.sig_ow); 
fprintf('======================================================================================\n\n'); 
% 
%             figure('name','Primary Drainage & Secondary Imbibition Cappilary Pressure Curves',...
%                 'units','normalized','outerposition',[0 0 1 1])         
%             subplot(2,2,[1 3]);
%             grid on
%             plot(network.DrainageData(:,1),network.DrainageData(:,2),'-b')  
%             hold on             
%             run('CARB_cycle1_drain.m')
%             plot(Res_draincycle_1(:,1),Res_draincycle_1(:,2),'-r') 
%             hold on         
%             run('CARB_draincycle_1.m')
%             plot(Res_draincycle_1(:,1),Res_draincycle_1(:,2),'-g')            
%             legend('MatlabPNM Drainage','Location','North','Imperial Raeini Drainage','Location','North','Imperial Valvatne Drainage','Location','North') 
%             title('Drainage Cappilary Pressure Curves')
%             xlabel('Sw')
%             xlim([0 1])
%             ylabel('Pc (Pa)')   
%                    
%             subplot(2,2,[2 4]);
%             grid on
%             plot(network.ImbibitionData(:,1),network.ImbibitionData(:,2),'-b')  
%             hold on             
%             run('CARB_cycle2_imb.m')
%             plot(Res_imb_2(:,1),Res_imb_2(:,2),'-r') 
%             hold on         
%             run('CARB_imbcycle_2.m')
%             plot(Res_imb_2(:,1),Res_imb_2(:,2),'-g')            
%             legend('MatlabPNM Imbibition','Location','North','Imperial Raeini Imbibition','Location','North','Imperial Valvatne Imbibition','Location','North') 
%             title('Imbibition Cappilary Pressure Curves')
%             xlabel('Sw')
%             xlim([0 1])
%             ylabel('Pc (Pa)')
% 
%             figure('name','Primary Drainage & Secondary Imbibition Relative Permeability Curves',...
%                 'units','normalized','outerposition',[0 0 1 1]) 
%             subplot(2,2,[1 2]);   
%             plot(network.DrainageData(:,1),network.DrainageData(:,3),'-b',network.DrainageData(:,1),network.DrainageData(:,4),'-b')     
%             hold on 
%             run('CARB_draincycle_1.m')
%             plot(Res_draincycle_1(:,1),Res_draincycle_1(:,3),'-g',Res_draincycle_1(:,1),Res_draincycle_1(:,4),'-g')
%             hold on 
%             run('CARB_cycle1_drain.m')
%             plot(Res_draincycle_1(:,1),Res_draincycle_1(:,3),'-r',Res_draincycle_1(:,1),Res_draincycle_1(:,4),'-r')         
%             xlabel('Sw')
%             xlim([0 1])
%             ylabel('Reative Permeability')
%             ylim([0 1])
%             legend('Water Relative Permeability MatlabPNM','Oil Relative PermeabilityMatlabPNM','Location','North',...
%                 'Water Relative Permeability Imperial Valvatne','Oil Relative Permeability Imperial Valvatne','Location','North',...
%                 'Water Relative Permeability Imperial Raeini','Oil Relative Permeability Imperial Raeini','Location','North') 
%             title('Drainage Relative Permeability Curves')
%             
%             subplot(2,2,[3 4]);   
%             plot(network.ImbibitionData(:,1),network.ImbibitionData(:,3),'-b',network.ImbibitionData(:,1),network.ImbibitionData(:,4),'-b')     
%             hold on 
%             run('CARB_cycle2_imb.m')
%             plot(Res_imb_2(:,1),Res_imb_2(:,3),'-r',Res_imb_2(:,1),Res_imb_2(:,4),'-r')
%             hold on 
%             run('CARB_imbcycle_2.m')
%             plot(Res_imb_2(:,1),Res_imb_2(:,3),'-g',Res_imb_2(:,1),Res_imb_2(:,4),'-g')         
%             xlabel('Sw')
%             xlim([0 1])
%             ylabel('Reative Permeability')
%             ylim([0 1])
%             legend('Water Relative Permeability MatlabPNM','Oil Relative PermeabilityMatlabPNM','Location','North',...
%                 'Water Relative Permeability Imperial Raeini','Oil Relative Permeability Imperial Raeini','Location','North',...
%                 'Water Relative Permeability Imperial Valvatne','Oil Relative Permeability Imperial Valvatne','Location','North') 
%             title('Imbibition Relative Permeability Curves')
            %==============================================================================================================
            figure('name','MatlabPNM Cappilary Pressure & Relative Permeability Curves',...
                'units','normalized','outerposition',[0 0 1 1])         
            subplot(2,2,[1 3]); 
            plot(network.DrainageData(:,1),network.DrainageData(:,2),'-r')  
            hold on               
            plot(network.ImbibitionData(:,1),network.ImbibitionData(:,2),'-b') 
            legend('MatlabPNM Drainage','Location','North','MatlabPNM Imbibition','Location','North')
            xlabel('Sw')
            xlim([0 1])
            ylabel('Pc (Pa)')   
            
            subplot(2,2,2);   
            plot(network.DrainageData(:,1),network.DrainageData(:,3),'-b',network.DrainageData(:,1),network.DrainageData(:,4),'-r')  
            hold on           
            xlabel('Sw')
            xlim([0 1])
            ylabel('Reative Permeability')
            ylim([0 1])
            legend('Water Relative Permeability Drainage','Oil Relative Permeability Drainage','Location','North') 
            title('Drainage Relative Permeability Curves')
            
            subplot(2,2,4); 
            plot(network.ImbibitionData(:,1),network.ImbibitionData(:,3),'-b',network.ImbibitionData(:,1),network.ImbibitionData(:,4),'-r')  
            hold on             
            xlabel('Sw')
            xlim([0 1])
            ylabel('Reative Permeability')
            ylim([0 1])
            legend('Water Relative Permeability Imbibition','Oil Relative Permeability Imbibition','Location','North') 
            title('Imbibition Relative Permeability Curves')
            
            %====================================================================================================================
            %             %% Plot 3D 
%             figure('name','Nodes Secondary Imbibition')
%              N = nan(obj.numberOfNodes ,4);
%              M = nan(obj.numberOfNodes ,3);
%              for i = 1:obj.numberOfNodes                 
%                  
%                   if obj.Nodes{i}.isInvaded  
%                   N (i,1:3)=[obj.Nodes{i}.x_coordinate,obj.Nodes{i}.y_coordinate,obj.Nodes{i}.z_coordinate];
%                   N (i,4)=obj.Nodes{i}.radius*20000000;
%                   if obj.Nodes{i}.isInlet
%                       M(i,:)=[0,1,0];
%                   else
%                       M(i,:) = [0,1,1];
%                       filledThroats = 0;
%                         for j = 1:obj.Nodes{i}.connectionNumber                            
%                             if (obj.Links{obj.Nodes{i}.connectedLinks(j)}.occupancy == 'A')                                 
%                                 filledThroats = filledThroats + 1;
%                             end
%                         end 
%                       if obj.Nodes{i}.ThresholdPressure_PoreBodyFilling > Pc_imb && filledThroats ~= 0
%                          M(i,:) = [0,0,0];
%                       end
%                   end 
%                   end
%              end 
%             figure('name',' Links Secondary Imbibition')
%              B = nan(obj.numberOfLinks ,4);
%              D = nan(obj.numberOfLinks ,3);
%              for i = 1:obj.numberOfLinks                 
%                   pore1Index = obj.Links{i}.pore1Index;
%                   pore2Index = obj.Links{i}.pore2Index;
%                   if obj.Links{i}.isInvaded  
%                   if obj.Links{i}.isInlet
%                       B (i,1:3)=[obj.Nodes{pore2Index}.x_coordinate,obj.Nodes{pore2Index}.y_coordinate,obj.Nodes{pore2Index}.z_coordinate];
%                       B (i,4)=obj.Links{i}.radius*20000000;                      
%                       D(i,:)=[0,1,0];
%                   else
%                       D(i,:) = [0,1,1]; 
%                       if obj.Links{i}.ThresholdPressure_SnapOff > Pc_imb %|| obj.Links{i}.ThresholdPressure_PistonLike > Pc_imb
%                          D(i,:) = [0,0,0];
%                       end
%                       if obj.Links{i}.isOutlet
%                           B (i,1:3)=[obj.Nodes{pore1Index}.x_coordinate+obj.Links{i}.linkLength/2,obj.Nodes{pore1Index}.y_coordinate,obj.Nodes{pore1Index}.z_coordinate];
%                           B (i,4)=obj.Links{i}.radius*20000000;
%                       else
%                           B (i,1:3)=[(obj.Nodes{pore1Index}.x_coordinate-obj.Nodes{pore2Index}.x_coordinate)/2,...
%                               (obj.Nodes{pore1Index}.y_coordinate+obj.Nodes{pore2Index}.y_coordinate)/2,...
%                               (obj.Nodes{pore1Index}.z_coordinate+obj.Nodes{pore2Index}.z_coordinate)/2];
%                           B (i,4)=obj.Links{i}.radius*20000000;
%                       end
%                           
%                   end 
%                   end
%              end
%             scatter3(B(:,1),B(:,2),B(:,3),B(:,4),D,'filled');
end