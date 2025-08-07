function [OCV_SOCs,OCV_Temps,OCV_values,Cellcapa]= Get_OCV(cycle)

%% data from excell sheet A-123 aging test, OCV_values matrix with len(OCV_SOC,OCV_Temp
%OCV_temps and OCV_SOCs are line vector

OCV_SOCs     = flip([100 95 90 80 70 60 50 40 30 20 15 10 5 2]);

switch cycle

    case '0' %cycle 0 (0 km)

        OCV_values   =flip([ 3.4469 3.3316 3.3292 3.3274 3.3039 3.2892 3.286 3.2838 3.2636 3.2338 3.2132 3.1967 3.1749 3.077]);
        OCV_Temps     =23;
        Cellcapa     = 2.455;

    case '8' %cycle 8 (300)
        OCV_values    = flip([3.4515 3.3325 3.3301 3.3282 3.3024 3.2894 3.286 3.2838 3.263 3.2343 3.2149 3.1975 3.1837 3.104]);
        OCV_Temps     =23;
        Cellcapa     = 2.4309;

end

end