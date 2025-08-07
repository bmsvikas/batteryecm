function [dydt] = ECMODE_Meatech(t,y,CurrentFit,k,ODE)

    % Current at t
    Currnow  = CurrentFit(t);
    
    pt = ODE.nRCpairs;
    k=k';
    
    % Rc pais voltages, using tau and R
    
      Tau = k(1:pt);
      Rc=k(pt+1:2*pt);
      
      B = (Currnow .* Rc )./  k(1:pt);
   
      dydt= -y./Tau + B  ;
    


end
