% AEROSPACE PROPULSION 4 - NEW EFFICIENT PASSENGER AIRCRAFT (NEPA)
% THERMODYNAMIC CYCLE DESIGN
% Written by S. Messina 2477336M
% James Watt School of Engineering
% University of Glasgow
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

clear all
clc
close all

% Display user message
fprintf('THERMODYNAMIC CYCLE DESIGN OF A TURBOFAN ENGINE \n\n');
% Define options 
fprintf('Please choose an option \n');
fprintf ('1. Analysis varying Overall Compressor ratio \n');
fprintf ('2. Analysis varying Turbine Entry Temperature (TET) \n');
fprintf ('3. Close the application \n');
option = 0;

%Make sure the option entered by the user is valid
while (option~=1 || option~=2 || option~=3)
    option = input('Choose an option: ', 's');
    fprintf('\n');
    % Check validity of the option entered by the user
    while ( not(isnumeric(str2num(option))) || isempty(str2num(option)))
        fprintf('You have entered a blank or a non-numeric value. Please, try again.\n');
        option = input('Choose an option: ','s');
        fprintf('\n');
    end %end while
    option = str2num(option);
    switch option
        case 1
            clear
            option = 0;
            % Data available
            ec = 0.915; % Polytropic Efficiency of IP and HP compressors
            et_h = 0.93; % HP Turbine polytropic efficiency
            et_l = 0.92; % LP Turbine polytropic efficiency
            ni = 0.98; % Intake Isentropic efficiency
            nb = 0.99; % Combustion burner efficiency
            nf = 0.91; % Fan Isentropic efficiency
            nn = 0.985; % Exit nozzle (core and bypass) isentropic efficiency
            p_drop = 0.04; % Tot pressure drop in combustion chamber
            nm = 0.985; % Mechanical Efficiency
            
            Qr = 43100; % kJ/kgK Calorific Value of Fuel
            
            % Air Properties 
            c_p = 1005; % J/kgK specific heat 
            gamma = 1.4;  % ratio of specific heats
            R = 287; % J/kgK gas constant
            
            % Gas Properties after Combustion
            c_pg = 1.148 * 10^3; % J/kgK specific heat 
            gamma_g = 1.333;
            
            % Ambient Conditions (Cruise)
            Mf = 0.85; % Flight Mach Number
            T_a= 218.9; % Ambient Temperature (K) at 35,000 ft
            p_a= 23.84; % Ambient (static) pressure (kPa) at 35,000 ft
            
            % Design Point Data
            F_T = 102.9136514; % Total Thrust (kN)
            sfc = 0.017671556; % Specific Fuel Consumption (kg/kNs)
            %*******************    OPTIMISATION PARAMETERS    ************************
            %pi_c_hp = 20; % HP Compressor Pressure Ratio (p_o3/p_o25)
            %pi_c_ip = 4; % IP Compressor Pressure Ratio (p_o25/p_o21)
            T_o4 = 1760; % TET (Turbine Entry Temperature)
            %pi_fc = 1.8; % Fan pressure Ratio for core (p_o21/p_o2)
            pi_fb = 2; % Fan pressure Ratio for bypass (p_o13/p_o2)
            
            pi_c_hp = 6 : (12-6)/40 : 14;
            pi_c_ip = 2.5;
            pi_fc = 1.6;
            n = size(pi_c_hp,2);

            B = 6:0.5:9;
            s = size(B,2);
            
            %**************************************************************************
            for j = 1:1:s
                for i = 1:1:n
                pi_c(i) = pi_c_hp(i) .* pi_c_ip .* pi_fc;
            
            Vf = Mf*sqrt(gamma*R*T_a);
            
            %% INTAKE
            
            % Stations Ambient to 2
            p_oa = p_a * (1 + ((gamma-1)/2)*Mf^2)^(gamma/(gamma-1)); % ambient stagnation pressure
            T_o2 = T_a * (1 + ((gamma-1)/2)*Mf^2); % Stagnation temperature at fan inlet
            T_o2s = ni * (T_o2-T_a) + T_a; % Ideal stag. temp. att fan inlet
            p_o2 = p_a * (T_o2s/T_a)^(gamma/(gamma-1)); % Tot pressure at fan inlet
            
            
            %% BYPASS DUCT
            
            % Stations 2 to 13 FAN BYPASS
            p_o13 = pi_fb * p_o2;
            T_o13 = T_o2 * (1 + (1/nf)*((p_o13/p_o2)^((gamma-1)/gamma) - 1));
            
            
            % Stations 13 to 19 NOZZLE BYPASS
            % Assuming no stagnation pressure loss through the bypass 
            p_o19 = p_o13;
            T_o19 = T_o13;
            r_critical_bp = 1/(1 -(1/nn)*((gamma - 1)/(gamma + 1)))^(gamma/(gamma-1));
            p_ratio_b = p_o19/p_a;
            
            if p_ratio_b < r_critical_bp 
                T_19s = T_o13 * (p_a / p_o13)^((gamma-1)/gamma);
                T_19 = T_o13 + (nn * (T_19s - T_o13));
                p_19 = p_a;
                V_19(i,j) = sqrt(2 * c_p * (T_o13-T_19));
                M_19 = V_19/sqrt(gamma * R * T_19);
            else % bypass nozzle is choked
                p_19 = p_o19/r_critical_bp;
                T_19 = T_o19 * (p_19/p_o19)^((gamma-1)/gamma);
                M_19 = 1;
                V_19(i,j) = sqrt(gamma * R * T_19);
            end
            
            
            %% FAN CORE
            
            % Stations 2 to 2.1
            p_o21 = p_o2 * pi_fc;
            T_o21 = T_o2 * (1 + (1/nf)*((p_o21/p_o2)^((gamma-1)/gamma) -1 ));
            
            %% IP COMPRESSOR
            
            % Stations 2.1 to 2.5
            p_o25 = p_o21 * pi_c_ip;
            T_o25 = T_o21 * pi_c_ip^((gamma-1)/(ec*gamma));
            
            %% HP COMPRESSOR
            
            % Stations 2.5 to 3
            p_o3 = p_o25 * pi_c_hp(i);
            T_o3 = T_o25 * pi_c_hp(i)^((gamma-1)/(ec*gamma));
            
            % Overall Compressor Pressure ratio
            %pi_c = pi_c_hp .* pi_c_ip; % p_o3/p_o21
            
            %% COMBUSTION
            
            % Stations 3 to 4
            p_o4 = (1-p_drop)*p_o3;
            f = ((c_pg * T_o4) - (c_p * T_o3))/((Qr*10^3)-(c_pg * T_o4));
            f_actual = f /nb ;
            
            %% HP TURBINE
            
            % Stations 4 to 4.5
            T_o45 = T_o4 - (c_p * (T_o3 - T_o25))/((1+f_actual)*nm*c_pg);
            p_o45 = p_o4 * (T_o45./T_o4)^(gamma_g/(et_h*(gamma_g - 1)));
            
            %% LP TURBINE
            
            % Stations 4.5 to 5
            T_o5 = T_o45 - (B(j)*(T_o13-T_o2) *c_p + c_p * (T_o25-T_o21) + ...
                c_p *(T_o21 - T_o2))/(nm * c_pg * (1+f_actual));
            p_o5 = p_o45 * (T_o5./T_o45)^(gamma_g / (et_l*(gamma_g - 1)));
            
            
            %Evaluate Critical Pressure Ratio
            r_critical = 1/(1-(1/nn)*((gamma_g - 1)/(gamma_g + 1)))^...
                (gamma_g/(gamma_g-1));
            p_ratio = p_o5/p_a;
            
            %% CORE NOZZLE
            
            if p_ratio < r_critical
                T_9s = T_o5 * (p_a / p_o5)^((gamma_g-1)/gamma_g);
                T_9 = T_o5 + (nn * (T_9s - T_o5));
                p_9 = p_a;
                V_9(i,j) = sqrt(2 * c_pg * (T_o5-T_9));
                M_9 = V_9/sqrt(gamma_g * R * T_9);
            else 
                % Flow is choked, i.e. M_9 = 1, exit pressure is p*
                p_9 = p_o5/r_critical;
                T_9 = T_o5 * (p_9./p_o5)^((gamma_g -1)/gamma_g);
                M_9 = 1;
                V_9(i,j) = sqrt(gamma_g * R * T_9);
            end
            
            %% ENGINE PERFORMANCE
            
            % Specific Thrust Core
            F_sc = (1/(B(j)+1)) * (((1+f_actual)*V_9(i,j)) - Vf) + ((1+f_actual)*R*T_9*...
                (p_9-p_a))/(V_9(i,j)*p_9*(B(j)+1));
            
            % Specific Thrust Bypass
            F_sb = (B(j)/(B(j)+1)) * (V_19(i,j) - Vf) + (R*T_19*B(j)*(p_19-p_a))/...
                (V_19(i,j)*p_19*(B(j)+1));
            
            % Total Specific Thrust
%             F_S(i,j) = (1/(B(j)+1))*((1+f_actual)*V_9-Vf)+(((1+f_actual)*R*T_9)/(V_9*p_9*(B(j)+1)))...
%                 *(p_9-p_a)+(B(j)/(B(j)+1))*(V_19-Vf)+((R*T_19*B(j))/(V_19*p_19*(B(j)+1)))*(p_19-p_a);

            % Total Specific Thrust
            F_S(i,j) = F_sc + F_sb;
            
            
            % Specific Fuel Consumption
            tsfc(i,j) =10^3 * f_actual/((B(j)+1)*F_S(i,j));
            
            
            % Propulsive Efficiency
            n_p(i,j) = (F_S(i,j)*Vf)/((F_S(i,j)*Vf) + 0.5*((1/(B(j)+1)*...
            (V_9(i,j)-Vf)^2 + (B(j)/(B(j)+1)*(V_19(i,j)-Vf)^2))));
    
            n_th(i,j) =((F_S(i,j)*Vf )+(0.5*((1/ (B(j)+1))*( (V_9(i,j)-Vf)^2)+...
                (B(j)/(B(j)+1) )*((V_19(i,j)-Vf)^2))))/((f_actual/(B(j)+1))*...
            Qr*1000); % Thermal Efficiency

            n_o(i,j)=n_p(i,j)*n_th(i,j); % Overall Efficiency 
            
            
            
            %% ENGINE MASS FLOW AND SIZE
            m_in = (F_T*10^3)/F_S(i,j); % mass flow at inlet
            m_c = (1/(B(j)+1)) * m_in; % mass flow in core
            m_b = (B(j)/(B(j)+1)) * m_in; % mass flow in bypass
            
            A_9 = ((1+f_actual)*m_c)/((p_9/(0.287*T_9)) * V_9(i,j));
            A_19 = m_b/((p_19/(0.287*T_19)) * V_19(i,j));
            
            % ADDITIONAL DESIGN PARAMETER
            tip_fan = 0.35; % Ratio of hub to tip radius
            
            % METHOD 1 
            %T_choked = T_o2 * (2/(gamma + 1));
            %p_choked = p_o2 * (2/(gamma + 1))^(gamma/(gamma - 1));
            %V_choked = sqrt(gamma * R * T_choked);
            %m_a = (p_choked * V_choked)/(0.287 * T_choked);
            %A_in = m_in/(0.88 * m_a);
            
            % METHOD 2
            M_2 = 0.6;
            T_2 = T_o2/(1 + ((gamma-1)/2)*M_2^2);
            p_2 = p_o2/(1 + ((gamma-1)/2)*M_2^2)^(gamma/(gamma-1));
            V_2 = sqrt(2 * c_p *(T_o2-T_2));
            A_in = m_in/((p_2/(0.287*T_2))*V_2);
            r_t = sqrt(A_in / (pi * (1 - tip_fan^2)));
            d_i(i,j) = 2 * r_t;
            d_h = tip_fan * 2 * r_t;
            
            
            v_exit_ratio(i,j) = V_19(i,j)/V_9(i,j);

            if B(j)==8 && i==76
                a = pi_fc;
                b = pi_c_ip;
                c = pi_c_hp(i);
                d = F_S(i,j)
                e = V_9(i,j)
                f = V_19(i,j)
            end


            
                end
            
            end
            
            figure
            nexttile
            plot(pi_c, tsfc, '-x')
            xlabel('Compressor Pressure Ratio, \pi_c')
            ylabel('Specific Fuel Consumption (kg/kNs)')
            yline(sfc, '-.', 'LineWidth',1.2);
            legendStrings = "B = " + string(B);
            legend(legendStrings)
            hold off
            
            nexttile
            plot(pi_c, F_S, '-*')
            xlabel('Compressor Pressure Ratio, \pi_c')
            ylabel('Specific Thrust (Ns/kg)')
            legendStrings = "B = " + string(B);
            legend(legendStrings)
            hold off
            
            nexttile
            plot(pi_c, d_i, '-x')
            xlabel('Compressor Pressure Ratio, \pi_c')
            ylabel('Fan Inlet Diameter, d_{\it i} (m)')
            yline(3.2, '-.', 'LineWidth',1.2)
            legendStrings = "B = " + string(B);
            legend(legendStrings)
            hold off
            
            nexttile
            plot(pi_c, v_exit_ratio, '-*')
            xlabel('Compressor Pressure Ratio, \pi_c')
            ylabel('Ratio of bypass to core exit velicities (V_{19}/V_{9})')
            %yline(0.85, '-.', 'LineWidth',1.2)
            yline(0.7, '-.', 'LineWidth',1.2)
            legendStrings = "B = " + string(B);
            legend(legendStrings)
            hold off
            
            figure
            nexttile
            plot(pi_c, n_p*100, '-*')
            xlabel('Compressor Pressure Ratio, \pi_c')
            ylabel('Propulsive Efficiency, \eta_p (%)')
            legendStrings = "B = " + string(B);
            legend(legendStrings)
            hold off
    
            nexttile
            plot(pi_c, n_th*100, '-*')
            xlabel('Compressor Pressure Ratio, \pi_c')
            ylabel('Thermal Efficiency, \eta_{th} (%)')
            legendStrings = "B = " + string(B);
            legend(legendStrings)
            hold off
            
            nexttile
            plot(pi_c, n_o*100, '-*')
            xlabel('Compressor Pressure Ratio, \pi_c')
            ylabel('Overall Efficiency, \eta_o (%)')
            legendStrings = "B = " + string(B);
            legend(legendStrings)
            hold off

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

        case 2
            clear
            option = 0;
            % Data available
            ec = 0.915; % Polytropic Efficiency of IP and HP compressors
            et_h = 0.93; % HP Turbine polytropic efficiency
            et_l = 0.92; % LP Turbine polytropic efficiency
            ni = 0.98; % Intake Isentropic efficiency
            nb = 0.99; % Combustion burner efficiency
            nf = 0.91; % Fan Isentropic efficiency
            nn = 0.985; % Exit nozzle (core and bypass) isentropic efficiency
            p_drop = 0.04; % Tot pressure drop in combustion chamber
            nm = 0.985; % Mechanical Efficiency
            
            Qr = 43100; % kJ/kgK Calorific Value of Fuel
            
            % Air Properties
            c_p = 1005; % J/kgK specific heat 
            gamma = 1.4;  % ratio of specific heats
            R = 287; % J/kgK gas constant
            
            % Gas Properties after Combustion
            c_pg = 1.148 * 10^3; % J/kgK specific heat 
            gamma_g = 1.333;
            
            % Ambient Conditions (Cruise)
        Mf = 0.85; % Flight Mach Number
        T_a= 218.9; % Ambient Temperature (K) at 35,000 ft
        p_a= 23.84; % Ambient (static) pressure (kPa) at 35,000 ft
        
        % Design Point Data
        F_T = 102.9136514; % Total Thrust (kN)
        sfc = 0.017671556; % Specific Fuel Consumption (kg/kNs)
        pi_fb = 1.8; % Fan pressure Ratio for bypass (p_o13/p_o2)
            
        pi_c_hp = 10.5;
        pi_c_ip = 2.5;
        pi_fc = 1.6;
        
        pi_c = pi_fc * pi_c_hp * pi_c_ip;
        T_o4 = 1600:5:1800;
        c = size(T_o4,2);
         
        B = 6:0.5:9;
        s = size(B,2);
        
        %**************************************************************************
        for j = 1:1:s
            for i = 1:1:c
        Vf = Mf*sqrt(gamma*R*T_a);
        
        %* 
        

        %% INTAKE
        
        % Stations Ambient to 2
        p_oa = p_a * (1 + ((gamma-1)/2)*Mf^2)^(gamma/(gamma-1)); % ambient stagnation pressure
        T_o2 = T_a * (1 + ((gamma-1)/2)*Mf^2); % Stagnation temperature at fan inlet
        T_o2s = ni * (T_o2-T_a) + T_a; % Ideal stag. temp. att fan inlet
        p_o2 = p_a * (T_o2s/T_a)^(gamma/(gamma-1)); % Tot pressure at fan inlet
        
        
        %% BYPASS DUCT
        
        % Stations 2 to 13 FAN BYPASS
        p_o13 = pi_fb * p_o2;
        T_o13 = T_o2 * (1 + (1/nf)*((p_o13/p_o2)^((gamma-1)/gamma) - 1));
        
        
        % Stations 13 to 19 NOZZLE BYPASS
        % Assuming no stagnation pressure loss through the bypass 
        p_o19 = p_o13;
        T_o19 = T_o13;
        r_critical_bp = 1/(1 -(1/nn)*((gamma - 1)/(gamma + 1)))^(gamma/(gamma-1));
        p_ratio_b = p_o19/p_a;
        
        if p_ratio_b < r_critical_bp 
            T_19s = T_o13 * (p_a / p_o13)^((gamma-1)/gamma);
            T_19 = T_o13 + (nn * (T_19s - T_o13));
            p_19 = p_a;
            V_19 = sqrt(2 * c_p * (T_o13-T_19));
            M_19 = V_19/sqrt(gamma * R * T_19);
        else % bypass nozzle is choked
            p_19 = p_o19/r_critical_bp;
            T_19 = T_o19 * (p_19/p_o19)^((gamma-1)/gamma);
            M_19 = 1;
            V_19 = sqrt(gamma * R * T_19);
        end
        
        
        %% FAN CORE
        
        % Stations 2 to 2.1
        p_o21 = p_o2 * pi_fc;
        T_o21 = T_o2 * (1 + (1/nf)*((p_o21/p_o2)^((gamma-1)/gamma) -1 ));
        
        %% IP COMPRESSOR
        
        % Stations 2.1 to 2.5
        p_o25 = p_o21 * pi_c_ip;
        T_o25 = T_o21 * pi_c_ip^((gamma-1)/(ec*gamma));
        
        %% HP COMPRESSOR
        
        % Stations 2.5 to 3
        p_o3 = p_o25 * pi_c_hp;
        T_o3 = T_o25 * pi_c_hp^((gamma-1)/(ec*gamma));
        
        % Overall Compressor Pressure ratio
        %pi_c = pi_c_hp .* pi_c_ip; % p_o3/p_o21
        
        %% COMBUSTION
        
        % Stations 3 to 4
        p_o4 = (1-p_drop)*p_o3;
        f = ((c_pg * T_o4(i)) - (c_p * T_o3))/((Qr*10^3)-(c_pg * T_o4(i)));
        f_actual = f /nb ;
        
        %% HP TURBINE
        
        % Stations 4 to 4.5
        T_o45 = T_o4(i) - (c_p * (T_o3 - T_o25))/((1+f_actual)*nm*c_pg);
        p_o45 = p_o4 * (T_o45./T_o4(i))^(gamma_g/(et_h*(gamma_g - 1)));
        
        %% LP TURBINE
        
        % Stations 4.5 to 5
        T_o5 = T_o45 - (B(j)*(T_o13-T_o2) *c_p + c_p * (T_o25-T_o21) + ...
            c_p *(T_o21 - T_o2))/(nm * c_pg * (1+f_actual));
        p_o5 = p_o45 * (T_o5./T_o45)^(gamma_g / (et_l*(gamma_g - 1)));
        
        
        %Evaluate Critical Pressure Ratio
        r_critical = 1/(1-(1/nn)*((gamma_g - 1)/(gamma_g + 1)))^...
            (gamma_g/(gamma_g-1));
        p_ratio = p_o5/p_a;
        
        %% CORE NOZZLE
        
        if p_ratio < r_critical
            T_9s = T_o5 * (p_a / p_o5)^((gamma_g-1)/gamma_g);
            T_9 = T_o5 + (nn * (T_9s - T_o5));
            p_9 = p_a;
            V_9 = sqrt(2 * c_pg * (T_o5-T_9));
            M_9 = V_9/sqrt(gamma_g * R * T_9);
        else 
            % Flow is choked, i.e. M_9 = 1, exit pressure is p*
            p_9 = p_o5/r_critical;
            T_9 = T_o5 * (p_9./p_o5)^((gamma_g -1)/gamma_g);
            M_9 = 1;
            V_9 = sqrt(gamma_g * R * T_9);
        end
        
        %% ENGINE PERFORMANCE
        
        % Specific Thrust Core
        F_sc = (1/(B(j)+1)) * (((1+f_actual)*V_9) - Vf) + ((1+f_actual)*R*T_9*...
            (p_9-p_a))/(V_9*p_9*(B(j)+1));
        
        % Specific Thrust Bypass
        F_sb = (B(j)/(B(j)+1)) * (V_19 - Vf) + (R*T_19*B(j)*(p_19-p_a))/...
            (V_19*p_19*(B(j)+1));
        
        % Total Specific Thrust
        F_S(i,j) = F_sc + F_sb;
        
        
        % Specific Fuel Consumption
        tsfc(i,j) =10^3 * f_actual/((B(j)+1)*F_S(i,j));
        
        
        % Propulsive Efficiency
        n_p(i,j) = 100*(F_S(i,j)*Vf)/((F_S(i,j)*Vf) + 0.5*((1/(B(j)+1) *...
        (V_9-Vf)^2 + (B(j)/(B(j)+1)*(V_19-Vf)^2))));

        n_th(i,j) =((F_S(i,j)*Vf )+(0.5*((1/ (B(j)+1))*( (V_9-Vf)^2)+...
            (B(j)/(B(j)+1) )*((V_19-Vf)^2))))/((f_actual/(B(j)+1))*...
        Qr*1000); % Thermal Efficiency
        n_o(i,j)=(F_S(i,j)*Vf)/((f/(B(j)+1))*Qr*1000); % Overall Efficiency 
        
        
        
        %% ENGINE MASS FLOW AND SIZE
        m_in = (F_T*10^3)/F_S(i,j); % mass flow at inlet
        m_c = (1/(B(j)+1)) * m_in; % mass flow in core
        m_b = (B(j)/(B(j)+1)) * m_in; % mass flow in bypass
        
        A_9 = ((1+f_actual)*m_c)/((p_9/(0.287*T_9)) * V_9);
        A_19 = m_b/((p_19/(0.287*T_19)) * V_19);
        d_9 = 2*sqrt(A_9)/pi;

        % Consider exit area as a single nozzle
        A_exit = A_19 + A_9;
        d_exit = 2*sqrt((A_9+A_19)/pi);
        
        % ADDITIONAL DESIGN PARAMETER
        tip_fan = 0.35; % Ratio of hub to tip radius
        
        % METHOD 2
        M_2 = 0.6;
        T_2 = T_o2/(1 + ((gamma-1)/2)*M_2^2);
        p_2 = p_o2/(1 + ((gamma-1)/2)*M_2^2)^(gamma/(gamma-1));
        V_2 = sqrt(2 * c_p *(T_o2-T_2));
        A_in = m_in/((p_2/(0.287*T_2))*V_2);
        r_t = sqrt(A_in / (pi * (1 - tip_fan^2)));
        d_i(i,j) = 2 * r_t;
        d_h = tip_fan * 2 * r_t;
        
        
        v_exit_ratio(i,j) = V_19/V_9;

        %*
        if B(j)==9 & T_o4(i)==1755 % Extract data for HP Turbine
                mass_t = m_c;
                p_in_hp_turbine = p_o4;
                p_out_hp_turbine = p_o45;
                T_out_hp_turbine = T_o45;
                target_no = n_o(i,j);
                target_nth = n_th(i,j);
                target_np = n_p(i,j);
                target_d_in = d_i(i,j);
                target_d_hub = d_h;
                inlet_area = A_in;
                inlet_diameter = d_i(i,j);
                %break
        end
        
        
            end
        
        end
        
        figure
        nexttile
        plot(T_o4, tsfc, '-x')
        xlabel('TET, T_{o4} (K)')
        ylabel('Specific Fuel Consumption (kg/kNs)')
        yline(sfc, '-.', 'LineWidth',1.2);
        legendStrings = "B = " + string(B);
        legend(legendStrings)
        hold off
        
        nexttile
        plot(T_o4, F_S, '-*')
        xlabel('TET, T_{o4} (K)')
        ylabel('Specific Thrust (Ns/kg)')
        legendStrings = "B = " + string(B);
        legend(legendStrings)
        hold off
        
        nexttile
        plot(T_o4, d_i, '-x')
        xlabel('TET, T_{o4} (K)')
        ylabel('Fan Inlet Diameter, d_{\it i} (m)')
        yline(3.2, '-.', 'LineWidth',1.2)
        legendStrings = "B = " + string(B);
        legend(legendStrings)
        hold off

        nexttile
        plot(T_o4, v_exit_ratio, '-*')
        xlabel('TET, T_{o4} (K)')
        ylabel('Ratio of bypass to core exit velicities (V_{19}/V_{9})')
        %yline(0.85, '-.', 'LineWidth',1.2)
        yline(0.7, '-.', 'LineWidth',1.2)
        legendStrings = "B = " + string(B);
        legend(legendStrings)
        hold off


        figure
        nexttile
        plot(T_o4, n_p, '-*')
        xlabel('TET, T_{o4} (K)')
        ylabel('Propulsive Efficiency, \eta_p (%)')
        legendStrings = "B = " + string(B);
        legend(legendStrings)
        hold off

        nexttile
        plot(T_o4, n_th, '-*')
        xlabel('TET, T_{o4} (K)')
        ylabel('Thermal Efficiency, \eta_{th} (%)')
        legendStrings = "B = " + string(B);
        legend(legendStrings)
        hold off

        nexttile
        plot(T_o4, n_o*100, '-*')
        xlabel('TET, T_{o4} (K)')
        ylabel('Overall Efficiency, \eta_o (%)')
        legendStrings = "B = " + string(B);
        legend(legendStrings)
        hold off

        case 3
            fprintf('Exiting the code... \n')
            break

        otherwise
            fprintf('You chose a wrong option. Please choose again. \n\n');
     end %end switch
end %end while loop

