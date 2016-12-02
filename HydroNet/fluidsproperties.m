function [ rho,mu,lamda,cp]= fluidsproperties(fluid, temperature, weight_fraction )
%fluidsproperties returns the properties of the fluids at a temperature
%   Calculate the properties aproximation given a temperature and a
%   concentration	of :
%               Air
% 				Water
% 				Steam
% 				Oil
% 				EthyleneGlycol
% 				DiethyleneGlycol
% 				TriethyleneGlycol
% 				PropyleneGlycol
% 				DipropyleneGlycol
% 				EthyleneGycol_Coolant
% 				DiethyleneGlycol_Coolant
% 				TriethyleneGlycol_Coolant
% 				PropyleneGlycol_Coolant
% 				DipropyleneGlycol_Coolant
%   Be careful when typing, if the name is not a matxh to any of the above
%   it won't work.


switch fluid
    
    case 'Oil'

    rho_Oil=-16572.74497 - 38.71871842 * temperature + 0.054641551 * temperature ^ 2 ...
        - 0.0000345327 * temperature ^ 3 + 4398.316892 * log(temperature);
	mu_Oil= 0.0000748 * exp(1005.2 / (temperature - 157.45));
	lamda_Oil=-13.78897843 - 0.028315003 * temperature + 0.0000380794 * temperature ^ 2 ...
    - 0.0000000225321 * temperature ^ 3 + 3.437962885 * log(temperature);
    cp_Oil=-34040.12013 - 71.24978467 * temperature + 0.107361063 * temperature ^ 2 ...
        - 0.0000663966 * temperature ^ 3 + 8670.542844 * log(temperature);
    rho=rho_Oil;
    mu=mu_Oil;
    lamda=lamda_Oil;
    cp=cp_Oil;
    
    case 'Water'
        if temperature>373.15
            temperature=373.15;
        end
         if temperature<273.15
            temperature=273.15;
        end

	rho_Water=1002.17 - 0.116189*(temperature -273.15)- (0.358024 * 10^-2)* (temperature-273.15)^2 + ...
        (0.373667 *10^-5)*(temperature-273.15)^3;
	LN_mu_Water= -3.758023 + 590.9808/((temperature-273.15) + 137.2645);
    mu_Water=exp(LN_mu_Water);
	lamda_Water=0.570990 + (0.167156 * 10^-2)*(temperature-273.15) -(0.609054* 10^-5)* (temperature-273.15)^2;
    T_m = 325.65;
    T1 = temperature-T_m;
    A1 = 4.1814288700;
    A2 = 0.0003202506;
    A3 = 7.309494E-06;
    A4 = 7.847426E-09;
    A5 = 1.327386E-09;
    A6 = -3.752825E-11;
    A7 = 3.691599E-13;
	cp_Water = A1 + A2*T1 + A3*T1^2 + A4*T1^3 + A5*T1^4 + A6*T1^5 + A7*T1^6;
    
    rho=rho_Water;
    mu=mu_Water;
    lamda=lamda_Water;
    cp=cp_Water*1000;
    
    case 'Steam'
        
        T_m = 625;
        T1 = temperature-T_m;
        A1 = 0.3505343750;
        A2 = -0.0005597328;
        A3 = 1.077121E-06;
        A4 = -1.941725E-09;
        rho_Steam=A1 + A2*T1 + A3*T1^2 + A4*T1^3;
        mu_Steam=exp(658.25*(1/temperature-1/283.16));
        A1 = 0.0442012500;
        A2 = 8.409538E-05;
        A3 = 1.303030E-08;
        A4 = -1.087801E-11;
        lamda_Steam=A1 + A2*T1 + A3*T1^2 + A4*T1^3;
        cp_Steam=7.701+4.595E-4*temperature+2.521E-6*temperature^2-0.859E-9*temperature^3;
               
        rho=rho_Steam;
        mu=mu_Steam/1000;
        lamda=lamda_Steam;
        cp=cp_Steam*4.1868*1000/18.015;

    case 'Air'
        
        R_Air=286.9; %J/kg*k R
    
	rho_Air=-6.77476761509903 + -4.53248600005648E-03 * temperature + 2.8723517898272E-06 * temperature ^ 2 ...
        - 7.9439617758404E-10 * temperature ^ 3 + 1.32166671126734 * log(temperature) ...
        + 458.230937786361 / temperature;
	mu_Air=0.0000014615 * (temperature ^ 1.5) / (temperature + 110.4);
	lamda_Air= -0.042199691369186 + 1.49433610434619E-05 * temperature + 5.00268427493186E-08 * temperature ^ 2 ...
        - 2.72058945048655E-11 * temperature ^ 3 + 1.04662104631254E-02 * log(temperature) + 0.151780981767906 / temperature;
	cp_Air=( 2522.88 - 10.4199 * sqrt(temperature) - 67227.1 / sqrt(temperature) + 917124.4 / temperature - 4174853.6 * temperature ^ (-1.5))+ R_Air;
    
    rho=rho_Air;
    mu=mu_Air;
    lamda=lamda_Air;
    cp=cp_Air;
    
    case 'EthyleneGlycol'
    
    rho_EtyhleneGlycol=1127.68-0.65816*temperature-6.1765*10^-4*temperature^2;
	LN_mu_EthyleneGlycol=-3.61359+986.519/(temperature+127.861);
    mu_EthyleneGlycol=exp(LN_mu_EthyleneGlycol);
	lamda_EthyleneGlycol=0.24658+2.5372*10^-4*temperature-1.3186*10^-6*temperature^2;
	 T_m = 50;
        X_m = 49.52;
        X = 100;
        T1 = temperature-T_m;
        X1 = X-X_m;
        A1 = 3.4496924300;
        A2 = -0.0175106131;
        A3 = -6.500687E-05;
        A4 = 3.509995E-08;
        A5 = 1.027679E-08;
        A6 = 1.410706E-10;
        A7 = 0.0052575507;
        A8 = 4.848820E-05;
        A9 = -1.744302E-06;
        A10 = 1.245819E-09;      
        A11 = 3.107839E-10;
        A12 = -2.730481E-05;
        A13 = 1.752306E-07;
        A14 = 1.323768E-08;
        A15 = -1.467148E-10;
        A16 = -5.182718E-09;
        A17 =  3.069004E-10;
        A18 = -6.120226E-12;
        cp_EthyleneGlycol= A1 + A2*X1 + A3*X1^2 + A4*X1^3 + A5*X1^4 + A6*X1^5 + A7*T1 + A8*T1*X1 ...        % Adjusted fromCalorimetric investigation of excess molar heat capacities for water + ethylene glycol fromT = 273.15 to T = 373.15 K
            + A9*T1*X1^2 + A10*T1*X1^3 + A11*T1*X1^4 + A12*T1^2 + A13*T1^2*X1 + A14*T1^2*X1^2 ...
            + A15*T1^2*X1^3 + A16*T1^3 + A17*T1^3*X1 + A18*T1^3*X1^2; 
    
    rho=rho_EtyhleneGlycol;
    mu=mu_EthyleneGlycol;
    lamda=lamda_EthyleneGlycol;
    cp=cp_EthyleneGlycol*1000;
    
    
    case 'DiethyleneGlycol'

	rho_DiethyleneGlycol=1132.35-0.67950*temperature-4.7565*10^-4*temperature^2;
	LM_mu_DiethyleneGlycol=-3.25001+901.095/(temperature+110.695);
    mu_DiethyleneGlycol=exp(LM_mu_DiethyleneGlycol);
	lamda_DiethyleneGlycol=0.19365+1.9938*10^-4*temperature-1.0584*10^-6*temperature^2;
	cp_DiethyleneGlycol=NaN;       
    
    rho=rho_DiethyleneGlycol;
    mu=mu_DiethyleneGlycol;
    lamda=lamda_DiethyleneGlycol;
    cp=cp_DiethyleneGlycol;
    
    case 'TriethyleneGlycol'

	rho_TriethyleneGlycol=1139.48-0.71040*temperature-4.3663*10^-4*temperature^2;
	LM_mu_TriethyleneGlycol=-3.11771+914.766/(temperature+110.068);
    mu_TriethyleneGlycol=exp(LM_mu_TriethyleneGlycol);
	lamda_TriethyleneGlycol=0.18890+1.1485*10^-4*temperature-8.4807*10^-7*temperature^2;       
	cp_TriethyleneGlycol=NaN;
    
    rho=rho_TriethyleneGlycol;
    mu=mu_TriethyleneGlycol;
    lamda=lamda_TriethyleneGlycol;
    cp=cp_TriethyleneGlycol;
    
    case 'PropyleneGlycol'

	rho_PropyleneGlycol=1003.7-0.20062*temperature-2.512E-3*temperature^2-147.12 ...
            -1.1024*temperature +2.6902E-3*temperature^2 ...
            -99.617 +0.63102*temperature ...
            -1.1267E-3*temperature^2;
	LN_mu_PropyleneGlycol=-3.9701+1000.8/(temperature+104.10);
    mu_PropyleneGlycol=exp(LN_mu_PropyleneGlycol);
	lamda_PropyleneGlycol=0.19116+1.1999*10^-4*temperature-9.2459*10^-7*temperature^2;
	T_m = 64.24;
    X_m = 45.59;
    X = weight_fraction*100;
    T1 = temperature-T_m;
    X1 = X-X_m;
    A1 = 3.8091256400;
    A2 = -0.0157459941;
    A3 = -0.0001604888;
    A4 = 1.678128E-06;
    A5 = 1.903529E-08;
    A6 = -3.416605E-10;
    A7 = 0.0029358287;
    A8 = 0.0001147783;
    A9 = -2.663112E-07;
    A10 = -2.184980E-08;
    A11 = 1.384410E-10;
    A12 = -2.012514E-06;
    A13 = 7.783426E-08;
    A14 = 4.042004E-09;
    A15 = -9.158986E-11;
    A16 = 7.375548E-08;
    A17 = -4.809393E-09;
    A18 = 5.103901E-11;
    cp_PropyleneGlycol= A1 + A2*X1 + A3*X1^2 + A4*X1^3 + A5*X1^4 + A6*X1^5 + A7*T1 + A8*T1*X1 ...        % Adjusted from http://dowac.custhelp.com/app/answers/detail/a_id/7470/~/propylene-glycols---specific-heat-values
        + A9*T1*X1^2 + A10*T1*X1^3 + A11*T1*X1^4 + A12*T1^2 + A13*T1^2*X1 + A14*T1^2*X1^2 ...
        + A15*T1^2*X1^3 + A16*T1^3 + A17*T1^3*X1 + A18*T1^3*X1^2;     
    
    rho=rho_PropyleneGlycol;
    mu=mu_PropyleneGlycol;
    lamda=lamda_PropyleneGlycol;
    cp=cp_PropyleneGlycol*4.1868*1000;
    
     case 'DipropyleneGlycol'

	rho_DipropyleneGlycol=1003.8-0.20207*temperature-2.5073E-3*temperature^2 ...
            +168.58-1.1338*temperature ...
            +2.3810*10^-3*temperature^2-134.05 ...
            +0.62288*temperature-5.8641*10^-4*temperature^2;
	LN_mu_DipropyleneGlycol=-3.6944+878.20/(temperature+84.119);
    mu_DipropyleneGlycol=exp(LN_mu_DipropyleneGlycol);
	lamda_DipropyleneGlycol=0.15428+7.0784*10^-5*temperature-6.1800*10^-7*temperature^2;
    T_m = 64.24;
    X_m = 45.59;
    X = 100;
    T1 = temperature-T_m;
    X1 = X-X_m;
    A1 = 4.0419793600;
    A2 = -0.0100217627;
    A3 = -0.0002168138;
    A4 = 1.461539E-06;  
    A5 =  3.410252E-08;
    A6 =  -6.763198E-10;
    A7 =  0.0009645226;
    A8 =  3.689494E-05;
    A9 =  1.240151E-06;
    A10 = -3.761298E-10;
    A11 = -3.041988E-10;
    A12 = -6.402469E-07;
    A13 =  3.419809E-08;
    A14 =  1.505691E-09;
    A15 = -3.158085E-10;
    A16 =  3.413751E-08;
    A17 =  1.732337E-09;
    A18 = -4.340647E-11;
	cp_DipropyleneGlycol= A1 + A2*X1 + A3*X1^2 + A4*X1^3 + A5*X1^4 + A6*X1^5 + A7*T1 + A8*T1*X1 + A9*T1*X1^2 ...        % Adjusted from http://dowac.custhelp.com/app/answers/detail/a_id/7470/~/propylene-glycols---specific-heat-values
        + A10*T1*X1^3 + A11*T1*X1^4 +  A12*T1^2 + A13*T1^2*X1 + A14*T1^2*X1^2 + A15*T1^2*X1^3 +  A16*T1^3 ...
        + A17*T1^3*X1 + A18*T1^3*X1^2;
    
    rho=rho_DipropyleneGlycol;
    mu=mu_DipropyleneGlycol;
    lamda=lamda_DipropyleneGlycol;
    cp=cp_DipropyleneGlycol*4.1868*1000;
    

    case 'EthyleneGlycol_Coolant'
         
        rho_Water=1002.17 - 0.116189*(temperature -273.15)- (0.358024 * 10^-2)* (temperature-273.15)^2 + ...
        (0.373667 *10^-5)*(temperature-273.15)^3;
        LN_mu_Water= -3.758023 + 590.9808/((temperature-273.15) + 137.2645);
        lamda_Water=0.570990 + (0.167156 * 10^-2)*(temperature-273.15) -(0.609054* 10^-5)* (temperature-273.15)^2;
    
        rho_EthyleneGlycol=1127.68-0.65816*temperature-6.1765*10^-4*temperature^2;
        LN_mu_EthyleneGlycol=-3.61359+986.519/(temperature+127.861);
        lamda_EthyleneGlycol=0.24658+2.5372*10^-4*temperature+-1.3186*10^-6*temperature^2;
        
        rho_EthyleneGlycolMix=weight_fraction*rho_EthyleneGlycol+(1-weight_fraction)*rho_Water...
            +(rho_EthyleneGlycol-rho_Water)*weight_fraction*(1-weight_fraction)*(0.30590+0.13781*weight_fraction-1.8961*10^-3*temperature);
        LN_mu_EthyleneGlycolMix=weight_fraction*LN_mu_EthyleneGlycol+(1-weight_fraction)*LN_mu_Water...
            +(LN_mu_EthyleneGlycol-LN_mu_Water)*weight_fraction*(1-weight_fraction)*(-0.165301-0.287325*weight_fraction+1.10978*10^-3*temperature);
        mu_EthyleneGlycolMix=exp(LN_mu_EthyleneGlycolMix);
        lamda_EthyleneGlycolMix=weight_fraction*lamda_EthyleneGlycol+(1-weight_fraction)*lamda_Water...
            +(lamda_EthyleneGlycol-lamda_Water)*weight_fraction*(1-weight_fraction)*(0.14219+0.38715*weight_fraction-6.6551*10^-4*temperature);
        T_m = 50;
        X_m = 49.52;
        X = weight_fraction*100;
        T1 = temperature-T_m;
        X1 = X-X_m;
        A1 = 3.4496924300;
        A2 = -0.0175106131;
        A3 = -6.500687E-05;
        A4 = 3.509995E-08;
        A5 = 1.027679E-08;
        A6 = 1.410706E-10;
        A7 = 0.0052575507;
        A8 = 4.848820E-05;
        A9 = -1.744302E-06;
        A10 = 1.245819E-09;      
        A11 = 3.107839E-10;
        A12 = -2.730481E-05;
        A13 = 1.752306E-07;
        A14 = 1.323768E-08;
        A15 = -1.467148E-10;
        A16 = -5.182718E-09;
        A17 =  3.069004E-10;
        A18 = -6.120226E-12;
        cp_EthyleneGlycolMix= A1 + A2*X1 + A3*X1^2 + A4*X1^3 + A5*X1^4 + A6*X1^5 + A7*T1 + A8*T1*X1 ...        % Adjusted fromCalorimetric investigation of excess molar heat capacities for water + ethylene glycol fromT = 273.15 to T = 373.15 K
            + A9*T1*X1^2 + A10*T1*X1^3 + A11*T1*X1^4 + A12*T1^2 + A13*T1^2*X1 + A14*T1^2*X1^2 ...
            + A15*T1^2*X1^3 + A16*T1^3 + A17*T1^3*X1 + A18*T1^3*X1^2;

        rho=rho_EthyleneGlycolMix;
        mu=mu_EthyleneGlycolMix;
        lamda=lamda_EthyleneGlycolMix;
        cp=cp_EthyleneGlycolMix*1000;
    
     case 'DiethyleneGlycol_Coolant'
         
        rho_Water=1002.17 - 0.116189*(temperature -273.15)- (0.358024 * 10^-2)* (temperature-273.15)^2 + ...
        (0.373667 *10^-5)*(temperature-273.15)^3;
        LN_mu_Water= -3.758023 + 590.9808/((temperature-273.15) + 137.2645);
        lamda_Water=0.570990 + (0.167156 * 10^-2)*(temperature-273.15) -(0.609054* 10^-5)* (temperature-273.15)^2;
    
        rho_DiethyleneGlycol=1132.35-0.67950*temperature-4.7565*10^-4*temperature^2;
        LM_mu_DiethyleneGlycol=-3.25001+901.095/(temperature+110.695);
        lamda_DiethyleneGlycol=0.19365+1.9938*10^-4*temperature-1.0584*10^-6*temperature^2;

        rho_DiethyleneGlycolMix=weight_fraction*rho_DiethyleneGlycol+(1-weight_fraction)*rho_Water...
            +(rho_DiethyleneGlycol-rho_Water)*weight_fraction*(1-weight_fraction)*(0.90820-0.26348*weight_fraction-3.3787*10^-3*temperature);
        LN_mu_DiethyleneGlycolMix=weight_fraction*LM_mu_DiethyleneGlycol+(1-weight_fraction)*LN_mu_Water...
            +(LM_mu_DiethyleneGlycol-LN_mu_Water)*weight_fraction*(1-weight_fraction)*(-0.364260+0.334513*weight_fraction+9.24070*10^-4*temperature);
        mu_DiethyleneGlycolMix=exp(LN_mu_DiethyleneGlycolMix);
        lamda_DiethyleneGlycolMix=weight_fraction*lamda_DiethyleneGlycol+(1-weight_fraction)*lamda_Water...
            +(lamda_DiethyleneGlycol-lamda_Water)*weight_fraction*(1-weight_fraction)*(0.63029+-0.19822*weight_fraction--1.2787*10^-3*temperature);
        cp_DiethyleneGlycolMix=NaN;

        rho=rho_DiethyleneGlycolMix;
        mu=mu_DiethyleneGlycolMix;
        lamda=lamda_DiethyleneGlycolMix;
        cp=cp_DiethyleneGlycolMix;
    
     case 'TriethyleneGlycol_Coolant'
         
        rho_Water=1002.17 - 0.116189*(temperature -273.15)- (0.358024 * 10^-2)* (temperature-273.15)^2 + ...
        (0.373667 *10^-5)*(temperature-273.15)^3;
        LN_mu_Water= -3.758023 + 590.9808/((temperature-273.15) + 137.2645);
        lamda_Water=0.570990 + (0.167156 * 10^-2)*(temperature-273.15) -(0.609054* 10^-5)* (temperature-273.15)^2;
    
        rho_TriethyleneGlycol=1139.48-0.71040*temperature-4.3663*10^-4*temperature^2;
        LM_mu_TriethyleneGlycol=-3.11771+914.766/(temperature+110.068);
        lamda_TriethyleneGlycol=0.18890+1.1485*10^-4*temperature-8.4807*10^-7*temperature^2;
        
        rho_TriethyleneGlycolMix=weight_fraction*rho_TriethyleneGlycol+(1-weight_fraction)*rho_Water...
            +(rho_TriethyleneGlycol-rho_Water)*weight_fraction*(1-weight_fraction)*(1.1712-0.52694*weight_fraction-3.8797*10^-3*temperature);
        LN_mu_TriethyleneGlycolMix=weight_fraction*LM_mu_TriethyleneGlycol+(1-weight_fraction)*LN_mu_Water...
            +(LM_mu_TriethyleneGlycol-LN_mu_Water)*weight_fraction*(1-weight_fraction)*(-0.727092+1.21086*weight_fraction-1.36642*10^-3*temperature);
        mu_TriethyleneGlycolMix=exp(LN_mu_TriethyleneGlycolMix);
        lamda_TriethyleneGlycolMix=weight_fraction*lamda_TriethyleneGlycol+(1-weight_fraction)*lamda_Water...
            +(lamda_TriethyleneGlycol-lamda_Water)*weight_fraction*(1-weight_fraction)*(0.11107+0.62975*weight_fraction-1.5995*10^-3*temperature);
        cp_TriethyleneGlycolMix=NaN;

        rho=rho_TriethyleneGlycolMix;
        mu=mu_TriethyleneGlycolMix;
        lamda=lamda_TriethyleneGlycolMix;
        cp=cp_TriethyleneGlycolMix;

     case 'PropyleneGlycol_Coolant'
         
        LN_mu_Water= -3.758023 + 590.9808/((temperature-273.15) + 137.2645);
        lamda_Water=0.570990 + (0.167156 * 10^-2)*(temperature-273.15) -(0.609054* 10^-5)* (temperature-273.15)^2;
    
        LN_mu_PropyleneGlycol=-3.9701+1000.8/(temperature+104.10);
        lamda_PropyleneGlycol=0.19116+1.1999*10^-4*temperature-9.2459*10^-7*temperature^2;
        
        rho_PropyleneGlycolMix=1003.7-0.20062*temperature-2.512E-3*temperature^2-147.12*weight_fraction ...
            -1.1024*temperature*weight_fraction +2.6902E-3*weight_fraction*temperature^2 ...
            -99.617*weight_fraction^2 +0.63102*weight_fraction^2*temperature ...
            -1.1267E-3*weight_fraction^2*temperature^2;
        LN_mu_PropyleneGlycolMix=weight_fraction*LN_mu_PropyleneGlycol+(1-weight_fraction)*LN_mu_Water...
            +(LN_mu_PropyleneGlycol-LN_mu_Water)*weight_fraction*(1-weight_fraction)*(1.5232-5.0007*weight_fraction-9.8106*10^-4*temperature+3.2452*weight_fraction^2);
        mu_PropyleneGlycolMix=exp(LN_mu_PropyleneGlycolMix);
        lamda_PropyleneGlycolMix=weight_fraction*lamda_PropyleneGlycol+(1-weight_fraction)*lamda_Water...
            +(lamda_PropyleneGlycol-lamda_Water)*weight_fraction*(1-weight_fraction)*(0.3622+09.0345*10^-2*weight_fraction--2.0935*10^-4*temperature);
        T_m = 64.24;
        X_m = 45.59;
        X = weight_fraction*100;
        T1 = temperature-T_m;
        X1 = X-X_m;
        A1 = 3.8091256400;
        A2 = -0.0157459941;
        A3 = -0.0001604888;
        A4 = 1.678128E-06;
        A5 = 1.903529E-08;
        A6 = -3.416605E-10;
        A7 = 0.0029358287;
        A8 = 0.0001147783;
        A9 = -2.663112E-07;
        A10 = -2.184980E-08;
        A11 = 1.384410E-10;
        A12 = -2.012514E-06;
        A13 = 7.783426E-08;
        A14 = 4.042004E-09;
        A15 = -9.158986E-11;
        A16 = 7.375548E-08;
        A17 = -4.809393E-09;
        A18 = 5.103901E-11;
        cp_PropyleneGlycolMix= A1 + A2*X1 + A3*X1^2 + A4*X1^3 + A5*X1^4 + A6*X1^5 + A7*T1 + A8*T1*X1 ...        % Adjusted from http://dowac.custhelp.com/app/answers/detail/a_id/7470/~/propylene-glycols---specific-heat-values
            + A9*T1*X1^2 + A10*T1*X1^3 + A11*T1*X1^4 + A12*T1^2 + A13*T1^2*X1 + A14*T1^2*X1^2 ...
            + A15*T1^2*X1^3 + A16*T1^3 + A17*T1^3*X1 + A18*T1^3*X1^2;

        rho=rho_PropyleneGlycolMix;
        mu=mu_PropyleneGlycolMix;
        lamda=lamda_PropyleneGlycolMix;
        cp=cp_PropyleneGlycolMix*4.1868*1000;

     case 'DipropyleneGlycol_Coolant'
         
        LN_mu_Water= -3.758023 + 590.9808/((temperature-273.15) + 137.2645);
        lamda_Water=0.570990 + (0.167156 * 10^-2)*(temperature-273.15) -(0.609054* 10^-5)* (temperature-273.15)^2;
    
        LN_mu_DipropyleneGlycol=-3.6944+878.20/(temperature+84.119);
        lamda_DipropyleneGlycol=0.15428+7.0784*10^-5*temperature-6.1800*10^-7*temperature^2;
        
        rho_DipropyleneGlycolMix=1003.8-0.20207*temperature-2.5073E-3*temperature^2 ...
            +168.58*weight_fraction-1.1338*weight_fraction*temperature ...
            +2.3810*10^-3*weight_fraction*temperature^2-134.05*weight_fraction^2 ...
            +0.62288*weight_fraction^2*temperature-5.8641*10^-4*weight_fraction^2*temperature^2;
        LN_mu_DipropyleneGlicolMix=weight_fraction*LN_mu_DipropyleneGlycol+(1-weight_fraction)*LN_mu_Water...
            +(LN_mu_DipropyleneGlycol-LN_mu_Water)*weight_fraction*(1-weight_fraction)*(1.6832-5.1426*weight_fraction-2.8529*10^-3*temperature+3.4638*weight_fraction^2);
        mu_DipropyleneGlycolMix=exp(LN_mu_DipropyleneGlicolMix);
        lamda_DipropyleneGlycolMix=weight_fraction*lamda_DipropyleneGlycol+(1-weight_fraction)*lamda_Water...
            +(lamda_DipropyleneGlycol-lamda_Water)*weight_fraction*(1-weight_fraction)*(0.40207+1.9196*10^-2*weight_fraction-2.3796*10^-4*temperature);
        T_m = 64.24;
        X_m = 45.59;
        X = weight_fraction*100;
        T1 = temperature-T_m;
        X1 = X-X_m;
        A1 = 4.0419793600;
        A2 = -0.0100217627;
        A3 = -0.0002168138;
        A4 = 1.461539E-06;
        A5 =  3.410252E-08;
        A6 =  -6.763198E-10;
        A7 =  0.0009645226;
        A8 =  3.689494E-05;
        A9 =  1.240151E-06;
        A10 = -3.761298E-10;
        A11 = -3.041988E-10;
        A12 = -6.402469E-07;
        A13 =  3.419809E-08;
        A14 =  1.505691E-09;
        A15 = -3.158085E-10;
        A16 =  3.413751E-08;
        A17 =  1.732337E-09;
        A18 = -4.340647E-11;
        cp_DipropyleneGlycolMix= A1 + A2*X1 + A3*X1^2 + A4*X1^3 + A5*X1^4 + A6*X1^5 + A7*T1 + A8*T1*X1 + A9*T1*X1^2 ...         % Adjusted from http://dowac.custhelp.com/app/answers/detail/a_id/7470/~/propylene-glycols---specific-heat-values
            + A10*T1*X1^3 + A11*T1*X1^4 +  A12*T1^2 + A13*T1^2*X1 + A14*T1^2*X1^2 + A15*T1^2*X1^3 +  A16*T1^3 ...
            + A17*T1^3*X1 + A18*T1^3*X1^2;

        rho =rho_DipropyleneGlycolMix;
        mu = mu_DipropyleneGlycolMix;
        lamda =lamda_DipropyleneGlycolMix;
        cp = cp_DipropyleneGlycolMix*4.1868*1000;
    
end

end

