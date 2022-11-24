function PLOTDRAIN(network)
%% Fluid info
fprintf('===============================Ntework Information=====================================\n');
fprintf('Water Viscosity:                        %f \n', network.waterViscosity); 
fprintf('Water Viscosity:                        %f \n', network.oilViscosity); 
fprintf('Interfacial Tension:                    %f \n', network.sig_ow); 
fprintf('======================================================================================\n\n'); 

%             figure('name','Primary Dranage Cappilary Pressure & Relative Permeability Curves',...
%                 'units','normalized','outerposition',[0 0 1 1])         
%             subplot(2,2,[1 3]);
%             grid on
%             plot(network.DrainageData(:,1),network.DrainageData(:,2),'-r')  
%             hold on             
%             run('CARB_cycle1_drain.m')
%             plot(Res_draincycle_1(:,1),Res_draincycle_1(:,2),'-g') 
%             hold on         
%             run('CARB_draincycle_1.m')
%             plot(Res_draincycle_1(:,1),Res_draincycle_1(:,2),'-k')            
%             legend('MatlabPNM','Location','North','Imperial Raeini','Location','North','Imperial Valvatne','Location','North') 
%             title('Drainage Cappilary Pressure Curves')
%             xlabel('Sw')
%             xlim([0 1])
%             ylabel('Pc (Pa)')   
%             
%             subplot(2,2,2);            
%             xlabel('Sw')
%             xlim([0 1])
%             ylabel('Reative Permeability')
%             ylim([0 1])
%             plot(network.DrainageData(:,1),network.DrainageData(:,3),'-r',network.DrainageData(:,1),network.DrainageData(:,4),'-r')     
%             hold on 
%             plot(Res_draincycle_1(:,1),Res_draincycle_1(:,3),'-g',Res_draincycle_1(:,1),Res_draincycle_1(:,4),'-g')
%             hold on 
%             run('CARB_cycle1_drain.m')
%             plot(Res_draincycle_1(:,1),Res_draincycle_1(:,3),'-k',Res_draincycle_1(:,1),Res_draincycle_1(:,4),'-k')
%             legend('Water Relative Permeability MatlabPNM','Oil Relative PermeabilityMatlabPNM','Location','North',...
%                 'Water Relative Permeability Imperial Valvatne','Oil Relative Permeability Imperial Valvatne','Location','North',...
%                 'Water Relative Permeability Imperial Raeini','Oil Relative Permeability Imperial Raeini','Location','North') 
%             title('Drainage Relative Permeability Curves')
  %==============================================================================================================
            figure('name','MatlabPNM Cappilary Pressure & Relative Permeability Curves',...
                'units','normalized','outerposition',[0 0 1 1])         
            subplot(2,2,[1 3]); 
            plot(network.DrainageData(:,1),network.DrainageData(:,2),'-r')   
            legend('MatlabPNM Drainage','Location','North')
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
             
            
end