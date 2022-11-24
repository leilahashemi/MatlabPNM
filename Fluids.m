classdef Fluids < handle
    %Fluids: contains information related to oil, gas and water 
    
    properties
        waterViscosity
        oilViscosity
        gasViscosity
        sig_ow
    end
    
    
    methods 
         function obj = Fluids()
             obj.waterViscosity = 0.001;
             obj.oilViscosity = 0.001; 
             obj.sig_ow = 10e-3; % N/m

         end
    end
end

