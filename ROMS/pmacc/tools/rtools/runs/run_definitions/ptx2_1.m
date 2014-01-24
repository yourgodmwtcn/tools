classdef ptx2_1 < run_parent
    % 6/5/2012 Parker MacCready
    
    methods
         
        function setScoord(rn)
            rn.theta_s = 4;
            rn.theta_b = 0.8;
            rn.tcline =  0;
            rn.N   =     27;
            rn.Vtransform = 1;
            rn.Vstretching = 1;
        end 
                        
    end % methods
    
end % classdef