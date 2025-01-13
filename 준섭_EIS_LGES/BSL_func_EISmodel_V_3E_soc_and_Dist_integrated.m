function [output,paras] = BSL_func_EISmodel_V_3E_soc_and_Dist_integrated(f_vector,factors,multi_soc_range,T,type_acf,cell_type,type_dist)

% Notes: 			This is an analytical model to predict the EIS for a full cell with intercalation based electrodes, separated by an ion conducting separator. The numerical equivalent was developed in COMSOL and compared to
% 					the analytical model in J. Electrochem. Soc. 2007 volume 154, issue 1,A43-A54 
% [corrected] 2008 papaer, Eqn 31: the first large parenthesis should have a negtive sign


% This version is based on JS_EIS_model_V6

omegagen= f_vector*(2*pi()); % frequency vector in [rad]
N = length(omegagen);% [v6]-taking from the input

A_coat = 12.6*2/10000; % % [m2] % LGES 2024 02

% type_anode = 0;
type_anode = evalin('base','type_anode');
 
addpath('C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\1_standalone\interpolations')

%% Thernodynamic Configs

% SOC and stoic 
    x_1 = 0.8781; % anode stoic when soc =0
    x_0 = 0.0216; % anode stoic when soc =1
    y_0 = 0.9319; % cathode stoic when soc =0
    y_1 = 0.3532; % cathode stoic when soc =1
    
   
    plateau1_start = 0.07; %for dUdc plateau setting
    plateau1_end = 0.15;
    
    plateau2_start = 0.18;
    plateau2_end = 0.50;
    
    plateau3_start = 0.51;
    plateau3_end = 0.98;

    %initiallized matrix
    c_imp = zeros(N,length(multi_soc_range));
    a_imp = zeros(N,length(multi_soc_range));
    s_imp = zeros(N,length(multi_soc_range));
    fc_imp = zeros(N,length(multi_soc_range));  

    for i = 1:length(multi_soc_range)
        soc = multi_soc_range(i)*0.01;

        x = x_0 + (x_1 - x_0)*soc; % anode stoic
        y = y_0 + (y_1 - y_0)*soc; % cathode stoic

%% Kinetics Parameters
    
    R_itsca(i) = factors(1,i)*0.0755; % [Ohm] % LGES 2024 02
    R_itscc(i) = factors(1,i+length(multi_soc_range))*0.0755; % [Ohm] % LGES 2024 02

    % Particle parameters
    Rfa = 0;                            Rfc = 0;        % {modified} [v6b - anode film *+*] [Ohm.m2]
    C_filma = 0;                        C_filmc = 0;    % {modified}  [v6b - anode film *+*] [F/m2]
    if type_anode == 0
    Rpa =  (17.8)*1e-6;              Rpc = (10)*1e-6;     % [um] Radius of particles  % LGES 2024 02
    elseif type_anode == 1
    Rpa =  (16.5)*1e-6;              Rpc = (10)*1e-6;     % [um] Radius of particles  % LGES 2024 02, Blending
    elseif type_anode == 2
    Rpa =  (17.8)*1e-6;              Rpc = (10)*1e-6;     % [um] Radius of particles  % LGES 2024 02, Natural
    end
    Dsa = factors(4,i)*Dsa_function(x,T);      Dsc = factors(4,i+length(multi_soc_range))*Dsc_function(y,T); % {modified} [m^2/sec]      
    Cdla =  factors(3,i)*0.2;             Cdlc = factors(3,i+length(multi_soc_range))*0.2;             % [F/m2]     
    cta = 29626;                        ctc = 48786;            % [mol/m3]

    % Exchange current density
    ka = factors(2,i)*ka_function(x,T,cta); % {modified}
    kc = factors(2,i+length(multi_soc_range))*kc_function(y,T,ctc); % {modified}
    c_e_ref = 1000; % {modified} [mol/m3] reference concetration of electrolyte

    % DRT Distribution parameters % LGES V2 2024 05 % HERE - ADD FACTORS
    mean_drt_ra = 1;                        mean_drt_rc = 1;
    std_drt_ra = factors(6,i);                       std_drt_rc = factors(6,i+length(multi_soc_range));      % [*+*] thse are parameters for observed lognormal distribution, dimensionless
    [mu_drt_ra,sig_drt_ra] = normal_para(mean_drt_ra,std_drt_ra); [mu_drt_rc,sig_drt_rc] = normal_para(mean_drt_rc,std_drt_rc); % converting them to parameters of underlying normal distribution

    % DDT Distribution parameters % LGES V2 2024 05 % HERE - ADD FACTORS
    mean_ddt_ra = 1;                        mean_ddt_rc = 1;
    std_ddt_ra = factors(7,i);                       std_ddt_rc = factors(7,i+length(multi_soc_range));      % [*+*] thse are parameters for observed lognormal distribution, dimensionless
    [mu_ddt_ra,sig_ddt_ra] = normal_para(mean_ddt_ra,std_ddt_ra); [mu_ddt_rc,sig_ddt_rc] = normal_para(mean_ddt_rc,std_ddt_rc); % converting them to parameters of underlying normal distributio

    % Porous electrode
    if type_anode == 0
    La = 79e-6;                         Lc = 60.0e-6;           % [m] LGES 2024 02
    elseif type_anode == 1
    La = 77e-6;                         Lc = 60.0e-6;           % [m] LGES 2024 02
    elseif type_anode == 2
    La = 79e-6;                         Lc = 60.0e-6;           % [m] LGES 2024 02
    end
    bruga = 1.44;                       brugc = 1.44;   % {modified}

    % epsla =   0.237;                    epslc = 0.234;         % {modified} void fraction LGES 2024 02

    if factors(5,2*length(multi_soc_range)+1) == 0 %for 3E_calc_el_multi_soc_V1
        epsla =   0.237;                    epslc = 0.234;         % {modified} void fraction LGES 2024 02
    elseif factors(5,2*length(multi_soc_range)+1) ~= 0 %for 3E_calc_el_multi_soc_V2 
        epsla =   factors(5,2*length(multi_soc_range)+1);     epslc = 0.234;          % {modified} void fraction LGES 2024 02
    end 

    n_am1_vf =    0.977;                p_am1_vf = 0.9792;     % {modified}LGES 2024 02
    epssa =   (1-epsla)*n_am1_vf;       epssc = (1-epslc)*p_am1_vf;  % {modified}
    taua = epsla^(-bruga);              tauc = epslc^(-brugc);    %{modifed} tortuosity of electrodes.
    if type_anode == 0
    sigmaa = 5.6022e+03;                     sigmac = 24.2718;             % [S/m] this is the solid phase conductivity
    elseif type_anode == 1
    sigmaa = 4.2088e+03;                     sigmac = 24.2718;            % [S/m] this is the solid phase conductivity
    elseif type_anode == 2
    sigmaa = 5.6022e+03;                     sigmac = 24.2718;            % [S/m] this is the solid phase conductivity
    end
    cs0a = x*cta;                       cs0c = y*ctc;     % OK
    alphaa = 0.5;                       % same alphaa           % OK                     
    alphac = 0.5;                       % same alphac           % OK


    % Electrolyte and Separator
    c_e = 1120;                 % {modified} [mol/m3] Initial electrolyte concentration
    Di0 = factors(2,2*length(multi_soc_range)+1)*De_function(c_e/1000,T);       % {modified} [m2/s] c_e input concentration in [mol/liter]
    epsls = 0.5;               % OK
    Lsep = 12.5e-6;              % OK % LGES 2024 02
    F = 96487;                  % OK
    R = 8.314;                  % OK
    iapp = 1;                   % OK - should not matter in impedance
    tplus = 0.363;              % OK
    nu=1;                       % OK
    % brugs = 3.0;              % not used anymore
    taus = 1.8;                % {modified} [1] tortuosity in separator
    % fc = 1.32;                
    % dfdx =1.7e-3;
    dlnfdlnc = (0.601-0.24*(c_e/(1000))^0.5+0.982*(1-0.0052*(T-298.15))*(c_e/(1000))^1.5)/(1-0.363)-1; % {modified} replacing f and dfdc
    kappa= factors(1,2*length(multi_soc_range)+1)*kappae_function(c_e/1000,T);                 % {modified} c_e input in [mol/liter] 


    %% 
% Parameter Expressions


    i0a = F*ka*((c_e/c_e_ref)^alphaa)*((cta-cs0a)^alphaa)*cs0a^alphac;                    % {modified} c_e_ref
    i0c = F*kc*((c_e/c_e_ref)^alphaa)*((ctc-cs0c)^alphaa)*cs0c^alphac;                    % {modified} c_e_ref
    
    a_n1 = factors(1,2*length(multi_soc_range)+2); %a_n is stand for anode expansion gradient
    a_n2 = factors(2,2*length(multi_soc_range)+2);
    a_n3 = factors(3,2*length(multi_soc_range)+2);
    a_n4 = factors(4,2*length(multi_soc_range)+2);
    a_n5 = factors(5,2*length(multi_soc_range)+2);

    a_p = factors(6,2*length(multi_soc_range)+2); %a_p is stand for cathode expansion gradient

    % aa =factors(5,i)*3*epssa/Rpa;   % {modifed} [m2/m3] this is specific area per a thickness % *+*
    % ac =factors(5,i+length(multi_soc_range))*3*epssc/Rpc;
    
    if factors(3,2*length(multi_soc_range)+1) == 0 %3E_calc_el_multi_soc_V1
       aa =factors(5,i)*3*epssa/Rpa;   % {modifed} [m2/m3] this is specific area per a thickness % *+*
       ac =factors(5,i+length(multi_soc_range))*3*epssc/Rpc;

    elseif factors(3,2*length(multi_soc_range)+1) ~= 0 %for 3E_calc_el_multi_soc_V2
       aa =factors(3,2*length(multi_soc_range)+1)*3*epssa/Rpa;   % {modifed} [m2/m3] this is specific area per a thickness % *+*
       ac =factors(4,2*length(multi_soc_range)+1)*3*epssc/Rpc;

       % if (soc < 0.12) % Anode specific sueface area decreases upon soc
       %    aa = aa*(1 - a_n1*soc); 
       % elseif (0.12 <= soc) && (soc < 0.18)
       %    aa = aa*(1 - a_n2*soc); 
       % elseif (0.18 <= soc) && (soc < 0.24)
       %    aa = aa*(1 - a_n3*soc); 
       % elseif (0.24 <= soc) && (soc < 0.50)
       %    aa = aa*(1 - a_n4*soc);  
       % elseif (0.50 <= soc)
       %    aa = aa*(1 - a_n5*soc);  
       % end 
       
       aa = aa*(1 - a_n5*soc);
       ac = ac*(1 + a_p*soc); % Cathode specific sueface area increases upon soc
    end

    sigmaeffa =(epssa/taua)*sigmaa; % {modified} all Bruggman relationships are modified.
    sigmaeffc =(epssc/tauc)*sigmac;
    
    kappaeffa = (epsla/taua)*kappa;
    kappaeffc = (epslc/tauc)*kappa;
    kappaeffs = (epsls/taus)*kappa;
    
    Dieffa = (epsla/taua)*Di0;
    Dieffc = (epslc/tauc)*Di0;
    Dieffs = (epsls/taus)*Di0;
    
    dx = 0.0001; % finite difference step size.
    chg = 0.5; % amount weighting on charging curve wrpt discharging.
    dUdcc = (1/ctc)*(Uc_function_v2(y+dx,chg) - Uc_function_v2(y-dx,chg))/(2*dx);   % {modified}
    dUdca = (1/cta)*(Ua_function_v2(x+dx,chg) - Ua_function_v2(x-dx,chg))/(2*dx);    % *+*
            
    if (plateau1_start <= soc) && (soc <= plateau1_end) %soc 0.07 - 0.15
       dUdca = (1/cta)*(Ua_function_v2(plateau1_end,chg) - Ua_function_v2(plateau1_start,chg))/(plateau1_end-plateau1_start);
    elseif (plateau2_start <= soc) && (soc <= plateau2_end) % soc 0.18 - 0.50
       dUdca = (1/cta)*(Ua_function_v2(plateau2_end,chg) - Ua_function_v2(plateau2_start,chg))/(plateau2_end-plateau2_start);
    elseif (plateau3_start <= soc) && (soc <= plateau3_end) % 0.51 - 00.98
       dUdca = (1/cta)*(Ua_function_v2(plateau3_end,chg) - Ua_function_v2(plateau3_start,chg))/(plateau3_end-plateau3_start);
    end   
 
    % if (0.07 <= soc) && (soc <= 0.1) %soc 0.07 - 0.1
    %    dUdca = (1/cta)*(Ua_function_v2(0.1,chg) - Ua_function_v2(0.07,chg))/(0.1-0.07);
    % elseif (0.22 <= soc) && (soc <= plateau2_end) % soc 0.18 - 0.50
    %    dUdca = (1/cta)*(Ua_function_v2(plateau2_end,chg) - Ua_function_v2(0.22,chg))/(plateau2_end-0.22);
    % elseif (plateau3_start <= soc) && (soc <= 0.99) % 0.51 - 00.98
    %    dUdca = (1/cta)*(Ua_function_v2(0.99,chg) - Ua_function_v2(plateau3_start,chg))/(0.99-plateau3_start);
    % end 
%% Calculation


for k = 1:N
%*************************************************************************
% There are three types of dimensionless laplace variables defined in each
% region

%sa/sc - corresponds to the anode/cathode%      
sa=1i*omegagen(k)*epsla*La^2/Dieffa;          %[refer to list of symbols]        
sc=1i*omegagen(k)*epslc*Lc^2/Dieffc;          %[refer to list of symbols]                    
ss=1i*omegagen(k)*epsls*Lsep^2/Dieffs;        %[refer to list of symbols]        

%s2a/s2c - corresponds to the anode/cathode particle
spa=1i*omegagen(k);      %   s2a                 %[refer to list of symbols]        
spc=1i*omegagen(k);      %   s2c                 %[refer to list of symbols]        

%*************************************************************************

Rcta=R*T/i0a/F/(alphaa+alphac);                 %[refer to list of symbols]        
Rctc=R*T/i0c/F/(alphaa+alphac);                 %[refer to list of symbols]        

paraa = [dUdca,Rpa,Dsa,F,spa,Cdla,Rcta,C_filma,Rfa,mu_drt_ra,sig_drt_ra,mu_ddt_ra,sig_ddt_ra]; % LGES V2 2024 05
parac = [dUdcc,Rpc,Dsc,F,spc,Cdlc,Rctc,C_filmc,Rfc,mu_drt_rc,sig_drt_rc,mu_ddt_rc,sig_ddt_rc]; % LGES V2 2024 05

betaa = beta_func(paraa,type_dist); % LGES V2 2024 05
betac = beta_func(parac,type_dist); % LGES V2 2024 05

% %for comparison beta & zloc of integral & equation
% betaa_sum(1:2,k) = betaa;
% betac_sum(1:2,k) = betac;
% 
% betaa = betaa(1,:);
% betac = betac(1,:);
% %-----------end

% local impedance [m2/ohm] times (1/F) 
%1/betaa(c)/F=Z_p,i in the manuscript (refer to eqn [9])

%*******************************************************
% For Meyers Case III, the simplification indicates that Rf should be
% replaced with another term( please see derivation in handwork stuff)

%Rf = Rf + (1/((s2*Cdl2)+(1/Rctc)))

%*******************************************************


B1a = (1-tplus)*La/F/Dieffa*(2*R*T*(1 - tplus)/F*(1/c_e)*(1+dlnfdlnc))/(1/La/aa/F/betaa);    % {modified} to use dlnf/dlnc
B1c = (1-tplus)*Lc/F/Dieffc*(2*R*T*(1 - tplus)/F*(1/c_e)*(1+dlnfdlnc))/(1/Lc/ac/F/betac);    % {modified} to use dlnf/dlnc       

B2a = La*(1/sigmaeffa+1/kappaeffa)/(1/La/aa/F/betaa);   %[refer to list of symbols]        
B2c = Lc*(1/sigmaeffc+1/kappaeffc)/(1/Lc/ac/F/betac);   %[refer to list of symbols]        

B3 = 2*kappaeffs*R*T*(1 - tplus)/Dieffs/F*(1/c_e)*(1+dlnfdlnc);  % {modified} to use dlnf/dlnc        

%EIGEN VALUES

lambda1a = 1/2*(sa + B1a + B2a + sqrt(sa^2 + 2*B1a*sa - 2*sa*B2a + B1a^2 + 2* B1a*B2a + B2a^2));    % [Eqn [5]in paper]
lambda1c = 1/2*(sc + B1c + B2c + sqrt(sc^2 + 2*B1c*sc - 2*sc*B2c + B1c^2 + 2* B1c*B2c + B2c^2));    % [Eqn [5]in paper]
 
lambda2a = 1/2*(sa + B1a + B2a - sqrt(sa^2 + 2*B1a*sa - 2*sa*B2a + B1a^2 + 2*B1a*B2a + B2a^2));     % [Eqn [5]in paper]
lambda2c = 1/2*(sc + B1c + B2c - sqrt(sc^2 + 2*B1c*sc - 2*sc*B2c + B1c^2 + 2*B1c*B2c + B2c^2));     % [Eqn [5]in paper]
%----------------------------------------------<<<<<<<<<<< SEPARATOR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>------------------------------------------------------------%


% CONSTANTS FROM THE EQUATIONS IN THE POROUS SEPARATOR REGION


GAMMAA = -La^3*aa*iapp*(1-tplus)*betaa/sigmaeffa/nu/Dieffa/(lambda1a-lambda2a)*((1/sqrt(lambda2a)/sinh(sqrt(lambda2a))-1/sqrt(lambda1a)/sinh(sqrt(lambda1a)))...
+sigmaeffa/kappaeffa*(1/sqrt(lambda2a)/tanh(sqrt(lambda2a))-1/sqrt(lambda1a)/tanh(sqrt(lambda1a))));    % [Eqn [20] in paper]

GAMMAC = -Lc^3*ac*iapp*(1-tplus)*betac/sigmaeffc/nu/Dieffc/(lambda1c-lambda2c)*((1/sqrt(lambda2c)/sinh(sqrt(lambda2c))-1/sqrt(lambda1c)/sinh(sqrt(lambda1c)))...
+sigmaeffc/kappaeffc*(1/sqrt(lambda2c)/tanh(sqrt(lambda2c))-1/sqrt(lambda1c)/tanh(sqrt(lambda1c))));    % [Eqn [20] in paper]

MU = Lsep/sqrt(ss)/Dieffs;       %   [Eqn [20] in paper]

BETAA = La/Dieffa/(lambda1a-lambda2a)*((sa-lambda2a+B1a)/(sqrt(lambda1a))/tanh(sqrt(lambda1a))-(sa-lambda1a+B1a)/(sqrt(lambda2a))/tanh(sqrt(lambda2a)))+MU/tanh(sqrt(ss)); %  [Eqn [20] in paper]
BETAC = Lc/Dieffc/(lambda1c-lambda2c)*((sc-lambda2c+B1c)/(sqrt(lambda1c))/tanh(sqrt(lambda1c))-(sc-lambda1c+B1c)/(sqrt(lambda2c))/tanh(sqrt(lambda2c)))+MU/tanh(sqrt(ss)); %  [Eqn [20] in paper]

ZETAA = (GAMMAA+GAMMAC*MU/sinh(sqrt(ss))/BETAC)/(BETAA-(MU/sinh(sqrt(ss)))^2/BETAC);  % [Eqn [20] in paper]
ZETAC = (GAMMAC+GAMMAA*MU/sinh(sqrt(ss))/BETAA)/(BETAC-(MU/sinh(sqrt(ss)))^2/BETAA);  % [Eqn [20] in paper]


%%%%okokokokokko%%%%%

C_5 = B3/sqrt(ss)*(ZETAC/sinh(sqrt(ss))-ZETAA/tanh(sqrt(ss)));    %[Eqn[14] in paper]
% C_6 = B3*ZETAA/sqrt(ss);                                           %[Eqn[14] in paper]

c_sep_xs_0 = C_5;
c_sep_xs_1 = B3/sqrt(ss)*(ZETAC/tanh(sqrt(ss))-ZETAA/sinh(sqrt(ss)));

phi2_sep_xs_0 = Lsep*iapp/kappaeffs*((c_sep_xs_0-c_sep_xs_1)/iapp+1);     %[Combine Eqn [27] and the line above]
phi2_sep_xs_1 = 0;      %   reference point

%---------------------------------<<<<<<<<<<<CATHODE CATHODE CATHODE CATHODE >>>>>>>>>>>>>>>>>>>>>>>>---------------------------------------------------------

%---------------------------------------------------------------------

C_2_c= B1c*iapp/sqrt(lambda1c)/(lambda1c-lambda2c)*(sigmaeffc/kappaeffc-(sc+B1c-lambda2c)*sigmaeffc*ZETAC/iapp/Lc^2/ac/(1-tplus)/betac);        %   Eqn [7] in the paper
C_1_c=-B1c*iapp/sqrt(lambda1c)/(lambda1c-lambda2c)/sinh(sqrt(lambda1c))-C_2_c/tanh(sqrt(lambda1c));                                              %   Eqn [7] in the paper
C_4_c=-B1c*iapp/sqrt(lambda2c)/(lambda1c-lambda2c)*(sigmaeffc/kappaeffc-(sc+B1c-lambda1c)*sigmaeffc*ZETAC/iapp/Lc^2/ac/(1-tplus)/betac);       %   Eqn [7] in the paper
C_3_c= B1c*iapp/sqrt(lambda2c)/(lambda1c-lambda2c)/sinh(sqrt(lambda2c))-C_4_c/tanh(sqrt(lambda2c));                                               %   Eqn [7] in the paper

C_7_c=-Lc*iapp/sigmaeffc*(1-sc*Lc^2*ac*F*betac/sigmaeffc/lambda1c/lambda2c);
C_8_c=-Lc/sigmaeffc*((sc-lambda1c)/B1c*C_1_c*(1-Lc^2*ac*F*betac/sigmaeffc/lambda1c)+(sc-lambda2c)/B1c*C_3_c*(1-Lc^2*ac*F*betac/sigmaeffc/lambda2c));

GAMMALAMBDAC1=-B1c*iapp/sqrt(lambda1c)/(lambda1c-lambda2c)/tanh(sqrt(lambda1c))-C_2_c/sinh(sqrt(lambda1c));
GAMMALAMBDAC2= B1c*iapp/sqrt(lambda2c)/(lambda1c-lambda2c)/tanh(sqrt(lambda2c))-C_4_c/sinh(sqrt(lambda2c));

phi1x1c=-Lc^3*ac*F*betac/sigmaeffc^2*((sc-lambda1c)/B1c/lambda1c*(GAMMALAMBDAC1)+(sc-lambda2c)/B1c/lambda2c*(GAMMALAMBDAC2))+C_7_c+C_8_c;



%-----------------------------------------<<<<<<<<<<<<<<<<<< ANODE ANODE ANODE ANODE ANODE ANODE >>>>>>>>>>>>>>>>>>>>>>------------------------

C_1_a =  B1a*iapp/sqrt(lambda1a)/(lambda1a-lambda2a)*(1/tanh(sqrt(lambda1a))+sigmaeffa/kappaeffa/sinh(sqrt(lambda1a))-(sa+B1a-lambda2a)*nu*sigmaeffa*ZETAA/La^2/aa/iapp/betaa/(1-tplus)/sinh(sqrt(lambda1a)));
C_2_a = -B1a*iapp/sqrt(lambda1a)/(lambda1a-lambda2a);
C_3_a = -B1a*iapp/sqrt(lambda2a)/(lambda1a-lambda2a)*(1/tanh(sqrt(lambda2a))+sigmaeffa/kappaeffa/sinh(sqrt(lambda2a))-(sa+B1a-lambda1a)*nu*sigmaeffa*ZETAA/La^2/aa/iapp/betaa/(1-tplus)/sinh(sqrt(lambda2a)));
C_4_a =  B1a*iapp/sqrt(lambda2a)/(lambda1a-lambda2a);


GAMMALAMBDAA1= B1a*iapp/sqrt(lambda1a)/(lambda1a-lambda2a)/sinh(sqrt(lambda1a))+B1a*iapp/sqrt(lambda1a)/(lambda1a-lambda2a)/tanh(sqrt(lambda1a))*(sigmaeffa/kappaeffa -(sa+B1a-lambda2a)*nu*sigmaeffa*ZETAA/La^2/aa/iapp/(1-tplus)/betaa);
GAMMALAMBDAA2=-B1a*iapp/sqrt(lambda2a)/(lambda1a-lambda2a)/sinh(sqrt(lambda2a))-B1a*iapp/sqrt(lambda2a)/(lambda1a-lambda2a)/tanh(sqrt(lambda2a))*(sigmaeffa/kappaeffa -(sa+B1a-lambda1a)*nu*sigmaeffa*ZETAA/La^2/aa/iapp/(1-tplus)/betaa);


C_7_a=-La*iapp/sigmaeffa*(1-La^2*aa*F*betaa/iapp/sigmaeffa*((sa-lambda1a)/B1a/sqrt(lambda1a)*C_2_a+(sa-lambda2a)/B1a/sqrt(lambda2a)*C_4_a));
C_8_a=-La/sigmaeffa*((sa-lambda1a)/B1a*GAMMALAMBDAA1*(1-La^2*aa*F*betaa/sigmaeffa/lambda1a)+(sa-lambda2a)/B1a*GAMMALAMBDAA2*(1-La^2*aa*F*betaa/sigmaeffa/lambda2a))+Lsep*iapp/kappaeffs*(1+(c_sep_xs_0-c_sep_xs_1)/iapp)-C_7_a;

phi1x1a = - La^3*aa*F*betaa/sigmaeffa^2*((sa-lambda1a)*C_1_a/B1a/lambda1a+(sa-lambda2a)*C_3_a/B1a/lambda2a) + C_8_a;

%*************************************************************************

%------------Overall Cell Potential Drop (Sandwich)------------------------

c_imp(k,i) = -(phi1x1c-phi2_sep_xs_1)/iapp;
a_imp(k,i) = -(phi2_sep_xs_0-phi1x1a)/iapp;
s_imp(k,i) = -(phi2_sep_xs_1-phi2_sep_xs_0)/iapp;
fc_imp(k,i) = -(phi1x1c-phi1x1a)/iapp;

%cell_potentiala (k) = Zc;
%cell_potentialb(k) = -phi1x1a;

%lambda1set(k)=lambda1;
%lambda2set(k)=lambda2;
%C_1set(k)=C_1;C_2set(k)=C_2;
%C_3set(k)=C_3;C_4set(k)=C_4;
%B1_set(k)=B1;
%B2_set(k)=B2;
%sset(k)=s;
%Capprox (k) = (s-lambda1+B1);
%Meyers(k)=1/beta/F;
end
parasa(:,i) = [R_itsca(i), i0a, Cdla, Dsa, aa]';
parasc(:,i) = [R_itscc(i), i0c, Cdlc, Dsc, ac, std_drt_ra,std_ddt_ra, std_drt_rc,std_ddt_rc]';
    end 

Av_grad_dist = [a_n1 a_n2 a_n3 a_n4 a_n5 0 a_p]';
assignin('base',"Av_grad_dist",Av_grad_dist)
% toc;
% w_vector = omegagen/(2*pi()); % into [Hz] from [rad] [v6] - taking from input
%cathode = cathode_impedance(1:N);
%anode = anode_impedance(1:N);
%separator = separator_impedance(1:N);
%fullcell = full_cell_impedance(1:N);

%end
%****************NYQUIST PLOTS*********************************************
 
% data = [omegagen(1:N)',real(c_imp)'*1e4,imag(c_imp)'*1e4]
 
% simulated_data = [real(cell_potential(start:ende))'*1e4,-imag(cell_potential(start:ende))'*1e4]; 
%{
% These are in unit of [Ohm/m2]
plot(real(c_imp(1:N)),-imag(c_imp(1:N)),'r-','Linewidth',2)
hold on 
plot(real(s_imp(1:N)),-imag(s_imp(1:N)),'g-','Linewidth',2)
hold on
plot(real(a_imp(1:N)),-imag(a_imp(1:N)),'b-','Linewidth',2)
hold on
plot(real(fc_imp(1:N)),-imag(fc_imp(1:N)),'m-','Linewidth',2)
%hold on
%plot(real(c_imp(1:N)+a_imp(1:N)+s_imp(1:N)),-imag(c_imp(1:N)+a_imp(1:N)+s_imp(1:N)),'bo')

% Bode plots

figure(2); hold on;
subplot(2,1,1)
loglog(w_vector,abs(c_imp))
subplot(2,1,2)
semilogx(w_vector,angle(c_imp)/pi()*180)
%}

%{
% These are in unit of [Ohm/cm2]
plot(real(c_imp(1:N))*1e4,-imag(c_imp(1:N))*1e4,'r-')
hold on 
plot(real(s_imp(1:N))*1e4,-imag(s_imp(1:N))*1e4,'go')
hold on
plot(real(a_imp(1:N))*1e4,-imag(a_imp(1:N))*1e4,'b-')
hold on
plot(real(fc_imp(1:N))*1e4,-imag(fc_imp(1:N))*1e4,'m-')
hold on
plot(real(c_imp(1:N)+a_imp(1:N)+s_imp(1:N))*1e4,-imag(c_imp(1:N)+a_imp(1:N)+s_imp(1:N))*1e4,'bo')
%}


%% Output data - changed for fitting

% output = [w_vector.', fc_imp.', c_imp.', a_imp.', s_imp.'];
% output = [R_itsc+(1/A_coat)*real(fc_imp.'),(1/A_coat)*imag(fc_imp.')]; % unit is [Ohm] format of real matrix

if type_acf ==1 % anode
    for i = 1:length(multi_soc_range)
        output(:,2*i-1) = R_itsc(i)+(1/A_coat)*real(a_imp(:,i));
        output(:,2*i) = (1/A_coat)*imag(a_imp(:,i));
    end
    
    paras =[kappa, Di0]';
elseif type_acf ==2 % cathode
    for i = 1:length(multi_soc_range)
        output(:,2*i-1) = R_itsc(i)+(1/A_coat)*real(c_imp(:,i));
        output(:,2*i) = (1/A_coat)*imag(c_imp(:,i));
        
    end
    
    paras =[kappa, Di0]';
elseif type_acf == 3 % 3E_Sum
        
        for i = 1:length(multi_soc_range)
            output(:,2*i-1) = R_itsca(i)+(1/A_coat)*real(a_imp(:,i)) + (1/A_coat)*real(c_imp(:,i)) +2.*real(s_imp(:,i));

            output(:,2*i) = (1/A_coat)*imag(a_imp(:,i))+(1/A_coat)*imag(c_imp(:,i))+2.*imag(s_imp(:,i)) ;
            paras(:,i) = [parasa(:,i);parasc(:,i); kappa; Di0];
        end

elseif type_acf == 4 % 3E_Simul
       
        for i = 1:length(multi_soc_range)
            output(:,2*i-1) = R_itsca(i)+(1/A_coat)*real(a_imp(:,i));
            output(:,2*(i+length(multi_soc_range))-1) = R_itscc(i)+(1/A_coat)*real(c_imp(:,i));

            output(:,2*i) = (1/A_coat)*imag(a_imp(:,i));
            output(:,2*(i+length(multi_soc_range))) = (1/A_coat)*imag(c_imp(:,i));

            paras(:,i) = [parasa(:,i);parasc(:,i); kappa; Di0; epsla;]; 
        end
else
    error('select anode (1), cathode (2), 3E_Sum (3), 3E_Simul (4)')
end
 

end

function beta = beta_func(para,type_dist) % LGES V2 2024 05
% para = [dUdc,Rp,Ds,F,sp,Cdl,Rct,C_film,Rf,mu_drt,sig_drt,mu_ddt,sig_ddt]
    if para(11) > 0.001
        if type_dist == 0 % drt
            upper_func = @(t)lognpdf(t,para(10),para(11)).*zloc_func_drt(t,para,type_dist); upper = integral(upper_func,0.001,100);
            down_func = @(t)lognpdf(t,para(10),para(11)); down = integral(down_func,0.001,100);
            beta = (upper./down)^-1/para(4);
        elseif type_dist == 1 % ddt
            upper_func = @(t)lognpdf(t,para(12),para(13)).*yloc_func_ddt(t,para,type_dist); upper = integral(upper_func,0.001,100);
            down_func = @(t)lognpdf(t,para(12),para(13)); down = integral(down_func,0.001,100);
            beta = upper./down/para(4);
        elseif type_dist == 2 % Integrated
            upper_func_drt = @(t)lognpdf(t,para(10),para(11)).*zloc_func_drt(t,para,type_dist); upper_drt = integral(upper_func_drt,0.001,100);
            down_func_drt = @(t)lognpdf(t,para(10),para(11)); down_drt = integral(down_func_drt,0.001,100);

            upper_func_ddt = @(t)lognpdf(t,para(12),para(13)).*yloc_func_ddt(t,para,type_dist); upper_ddt = integral(upper_func_ddt,0.001,100);
            down_func_ddt = @(t)lognpdf(t,para(12),para(13)); down_ddt = integral(down_func_ddt,0.001,100);

          
            
          

                %(10-16 combining; impedance sum (not addmitance sum) sum in beta^-1)
                % beta is an addmitance (with unit adjusted.)
            % beta_drt = (upper_drt./down_drt)^-1/para(4);
            % beta_ddt = upper_ddt./down_ddt/para(4); 
            % beta_p2d = zloc_func_p2d(para)^-1/para(4);
            % 
            % zloc = (beta_drt)^-1 + (beta_ddt)^-1 - (beta_p2d)^-1; 
            % 
            % 
            % beta = zloc^-1; % equation

            


            %% beta separation 
            beta_drt = (upper_drt./down_drt)^-1/para(4);
            beta_ddt = upper_ddt./down_ddt/para(4); 

            beta = (beta_drt^-1 + beta_ddt^-1)^-1; 

            % beta = [beta1; beta2];
        end
    else
        beta = yloc_func_ddt(1,para)/para(4);
    end
end

function yloc = yloc_func_ddt(t,para,type_dist) % LGES V2 2024 05
% para = [dUdc,Rp,Ds,F,sp,Cdl,Rct,C_film,Rf,mu_drt,sig_drt,mu_ddt,sig_ddt]; 

% Rdif=-dUdc*Rp/Ds/F;
Rdif = - para(1)*para(2)/para(3)/para(4);
% Ys=(sqrt(sp*Rp^2/Ds)-tanh(sqrt(sp*Rp^2/Ds)))/tanh(sqrt(sp*Rp^2/Ds))
Ys=(sqrt(para(5)*t*(para(2)).^2/para(3))-tanh(sqrt(para(5)*t*(para(2)).^2/para(3))))./tanh(sqrt(para(5)*t*(para(2)).^2/para(3)));  
% zeta = sp*Cdl + 1/(Rct+Rdif/Ys);          %[Eqn [9] in the paper]


if type_dist == 1
zeta = para(5)*para(6) + 1./(para(7)+Rdif./Ys); %[Eqn [9] in the paper]
elseif type_dist == 2
zeta = 1./(Rdif./Ys); 
% zeta = para(5)*para(6) + 1./(para(7)+Rdif./Ys);
end

%New Beta based on Meyers Case 2
yloc = para(5)*para(8) + 1./(1./zeta+para(9));  % edit -GS  [Eqn[8] in the paper]


end

function zloc = zloc_func_drt(t,para,type_dist) % LGES V2 2024 05
% para = [dUdc,Rp,Ds,F,sp,Cdl,Rct,C_film,Rf,mu_drt,sig_drt,mu_ddt,sig_ddt]; 

% Rdif=-dUdc*Rp/Ds/F;

Rdif = - para(1)*para(2)/para(3)/para(4);
% Ys=(sqrt(sp*Rp^2/Ds)-tanh(sqrt(sp*Rp^2/Ds)))/tanh(sqrt(sp*Rp^2/Ds))
Ys=(sqrt(para(5)*(para(2)).^2/para(3))-tanh(sqrt(para(5)*(para(2)).^2/para(3))))./tanh(sqrt(para(5)*(para(2)).^2/para(3)));  
% zeta = sp*Cdl + 1/(Rct+Rdif/Ys);          %[Eqn [9] in the paper]

if type_dist == 0
zeta = para(5)*para(6)*t + 1./(para(7)+Rdif./Ys); %[Eqn [9] in the paper]
elseif type_dist == 2
zeta = para(5)*para(6)*t + 1./para(7);          %[Eqn [9] in the paper] omit Rdif
% zeta = para(5)*para(6)*t + 1./(para(7)+Rdif./Ys);
end 

yloc = para(5)*para(8) + 1./(1./zeta+para(9));  % edit -GS  [Eqn[8] in the paper]

zloc = yloc.^-1;

end 


function zloc =  zloc_func_p2d(para)
% Rdif = -dUdc*Rp/Ds/F;
Rdif = -para(1)*para(2)/para(3)/para(4);

% Ys =(sqrt(sp*Rp^2/Ds)-tanh(sqrt(sp*Rp^2/Ds)))/tanh(sqrt(sp*Rp^2/Ds))
Ys = (sqrt(para(5)*para(2)^2/para(3))-tanh(sqrt(para(5)*para(2)^2/para(3))))/tanh(sqrt(para(5)*para(2)^2/para(3)));

% zeta = sp*Cdl + 1/(Rct+Rdif/Ys);
zeta = para(5)*para(6)+1/(para(7)+Rdif/Ys);

% admittance
yloc = (para(5)*para(8)+1/(1/zeta+para(9)));

zloc = yloc^-1;
end

function [mu,sig]=normal_para(mean,std) % LGES V2 2024 05
% converting parameters of observed lognormal distribution
% to parameters of underlying normal distribution
m=mean;
v=std^2;
mu=log(m^2/(v+m^2)^0.5);
sig=(log(v/m^2+1))^0.5;


end

