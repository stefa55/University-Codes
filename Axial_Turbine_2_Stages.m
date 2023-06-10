% AEROSPACE PROPULSION 4 - NEW EFFICIENT PASSENGER AIRCRAFT (NEPA)
% COMPONENT DESGIN - AXIAL TURBINE
% Written by S. Messina 2477336M
% James Watt School of Engineering
% University of Glasgow
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

clear all
clc
close all

%% Data & Assumptions

% DATA FROM THERMODYNAMIC CYCLE ANALYSIS
core_flow = 55.496877740179585; % Core mass flow
p_in = 1.528054871397713e+03; % Pressure at the inlet of HPT
p_exit = 5.665337120441504e+02; % Pressure at HPT outlet
T_exit = 1.393686082024633e+03; % Temperature at the exit of HPT
TET = 1755; % Turbine Entry Temperature

poly = 0.93; % HP Turbine polytropic efficiency
c_pg = 1148; % c_p of gases after combustion (J/kgK)
gamma_g = 1.333;
R = 287; 

% Define HPT data
Delta_T_stage1 = (TET - T_exit)*0.58;
Delta_T_stage2 = (TET-T_exit)*0.42;

T_exit1 = TET-Delta_T_stage1;
p_in1 = p_in;


p_exit1 = p_in1 * (T_exit1/TET)^(gamma_g/(poly*(gamma_g - 1)));


% Determine
% no. of HPT stages required
% stage loading psi, flow coefficient phi
% Reactions Lambda, Blade Angles, Flow Velocities and Blade Speeds at
% mid-height, hub and tip for each stage (free vortex blading)
% Annulus Dimensions: Flow Area, Blade Radii (hub, mid-height and tip),
% blade heights


% ASSUMPTIONS AND CONTROL PARAMETERS
alpha_1 = 0; % Flow is axial at turbine entry
% C_1 = C_3; % Absolute inlet and exit velocities to the stage are equal
% C_a2 = C_a3; % Axial velocity at rotor inlet and outlet is the same
flow_coeff = [0.4 0.6 0.8 1.0 1.2 1.4 1.6]; % Assume flow coefficient based on typical values
exit_swirl_alpha = 40;  % Assume exit swirl angle

shaft_speed = (8500/60); % Shaft speed in rev/s, has to be between 7500-9000rpm
blade_speed = 250:5:450; % Mean blade speed


%% STAGE 1
[ca_1,psi_1,am1,ar1,at1,bm1,br1,bt1,ur1,ut1,root_r1,tip_r1,Cm1,Cr1,Ct1,...
    Vm1,Vr1,Vt1,f1,l1,h1,M1, phi_chosen, U_chosen, hub_to_tip, mean_radius] = turbine_stage1(core_flow,...
    Delta_T_stage1, TET, T_exit1, exit_swirl_alpha,flow_coeff...
    ,blade_speed,shaft_speed,poly,p_in1,p_exit1);

%% STAGE 2

% Need to feed:
% alpha1= alpha 3 stage 1 at mid-height, root and tip
% beta1 = beta 3 stage 1 at mid-height, root and tip
% C1 = C3 stage 1 at mid-height, root and tip
% V1 = V3 stage 1 at mid-height, root and tip
% p_o1 = p_o3 stage 1
% p_o3 = p_exit

alpha_1m_s2 = am1(2);
alpha_1r_s2 = ar1(2);
alpha_1t_s2 = at1(2);
beta_1m_s2 = bm1(2);
beta_1r_s2 = br1(2);
beta_1t_s2 = bt1(2);
C_1m_s2 = Cm1(3);
C_1r_s2 = Cr1(2);
C_1t_s2 = Ct1(2);
V_1m_s2 = Vm1(2);
V_1r_s2 = Vr1(2);
V_1t_s2 = Vt1(2);

T_in2 = T_exit1;
T_exit2 = T_exit;
p_in2 = p_exit1;
p_exit2 = p_exit;

% Control Parameters
exit_swirl_alpha2 = 14;


[ca_2,psi_2,am2,ar2,at2,bm2,br2,bt2,ur2,ut2,root_r2,tip_r2,Cm2,Cr2,Ct2,...
    Vm2,Vr2,Vt2,f2,l2,h2,M2] = turbine(core_flow, Delta_T_stage2, T_in2, T_exit2, exit_swirl_alpha2,phi_chosen...
    ,U_chosen,shaft_speed,poly,p_in2,p_exit2, alpha_1m_s2, beta_1m_s2,...
    C_1m_s2, C_1r_s2, C_1t_s2, V_1m_s2, V_1r_s2, V_1t_s2);



%% FUNCTIONS
function [C_a_s1,psi_s1,alpha_m_s1,alpha_r_s1,alpha_t_s1,beta_m_s1,...
    beta_r_s1,beta_t_s1,Ur_s1,Ut_s1,rr_s1,rt_s1,C_m_s1,C_r_s1,...
    C_t_s1,V_m_s1,V_r_s1,V_t_s1,phi_s1,Lambda_s1,h_s1,M_2a_s1, phi_ext, U_ext, h_to_t, r_m_s1]...
    = turbine_stage1(m_in, T_drop, T_o1, T_o3, alpha_3m, phi, U, N, et_h ,p_o1, p_o3)

c_pg = 1148; % c_p of gases after combustion (J/kgK)
gamma_g = 1.333;
R = 287; 
s = size(U,2);

for j = 1:1:7

for i=1:1:s

%% ANGLES DEFINITION 
psi(i) = (c_pg*T_drop)/U(i)^2; % Stage loading

beta_3m(j) = atand((1/phi(j)) + tand(alpha_3m)); % Calculate relative angle at rotor exit
Lambda(i,j) = 0.5*((2*phi(j)*tand(beta_3m(j))) - psi(i)); % Reaction at mid-blade height

beta_2m(i,j) = atand(tand(beta_3m(j)) - ((2*Lambda(i,j))/phi(j)));
alpha_2m(i,j) = atand((1/phi(j)) + tand(beta_2m(i,j)));

%% VELOCITY TRIANGLES
C_a2(i,j) = phi(j)*U(i); % Axial Velocity at Interstage (2)
C_2(i,j) = C_a2(i,j)/cosd(alpha_2m(i,j)); % Absolute Velocity at Interstage (2)
T_2(i,j) = T_o1 - C_2(i,j)^2/(2*c_pg); % Static Temperature at NGV Exit

NGV_pratio = (T_o1/T_2(i,j))^(gamma_g / (et_h*(gamma_g - 1)));

r_crit = ((gamma_g+1)/2)^(gamma_g/(et_h*(gamma_g-1)));

if NGV_pratio > r_crit
    % Check for choked flow
    p_2(i,j) = p_o1/r_crit;
else 
    p_2(i,j) = p_o1/NGV_pratio;
end

A_2(i,j) = ((R*T_2(i,j))/(p_2(i,j)*1000)) * m_in/C_a2(i,j);
C_a3(i,j) = C_a2(i,j); % Axial Velocity at Rotor Exit equals the value at Rotor Inlet
C_3(i,j) = C_a3(i,j)/(cosd(alpha_3m));
C_1 (i,j)= C_3(i,j);
C_a1(i,j) = C_1(i,j);

% At turbine Inlet
T_1(i,j) = T_o1 - C_1(i,j)^2/(2*c_pg);
p_1(i,j) = p_o1 * (T_1(i,j)/T_o1)^(gamma_g/(gamma_g - 1));
A_1(i,j) = ((R*T_1(i,j))/(p_1(i,j)*1000)) * m_in/C_a1(i,j);


% At Stage Outlet
T_3(i,j) = T_o3 - C_3(i,j)^2/(2*c_pg);
p_3(i,j) = p_o3*(T_3(i,j)/T_o3)^(gamma_g/(gamma_g-1));
A_3(i,j) = ((R*T_3(i,j))/(p_3(i,j)*1000)) * m_in/C_a3(i,j);


% Relative Velocities at mid-height
V_2(i,j) = C_a2(i,j)/cosd(beta_2m(i,j));
V_3(i,j) = C_a3(i,j)/cosd(beta_3m(j));



%% BLADE STUDY
   
r_m(i) = U(i)/(2*pi*N); % Mean Radius

% Station 1
h_1(i,j) = A_1(i,j)*(N/U(i)); % Blade height
rt_1(i,j) = r_m(i) + (h_1(i,j)/2); % Tip radius
rr_1(i,j) = r_m(i) - (h_1(i,j)/2); % Root radius
ratiort_1(i,j) = rt_1(i,j)/rr_1(i,j); % Ratio of Tip to Root radius

ratior_1(i,j) = r_m(i)/rr_1(i,j) ; % Mean Radius to Root radius ratio
ratiot_1(i,j) = r_m(i)/rt_1(i,j) ; % Mean Radius to Tip radius ratio


% Station 2
h_2(i,j) = A_2(i,j)*(N/U(i)); % Blade height
rt_2(i,j) = r_m(i) + (h_2(i,j)/2); % Tip radius
rr_2(i,j) = r_m(i) - (h_2(i,j)/2); % Root radius
ratiort_2(i,j) = rt_2(i,j)/rr_2(i,j); % Ratio of Tip to Root radius

ratior_2(i,j) = r_m(i)/rr_2(i,j) ; % Mean Radius to Root radius ratio
ratiot_2(i,j) = r_m(i)/rt_2(i,j) ; % Mean Radius to Tip radius ratio

% Station 3
h_3(i,j) = A_3(i,j)*(N/U(i)); % Blade height
rt_3(i,j) = r_m(i) + (h_3(i,j)/2); % Tip radius
rr_3(i,j) = r_m(i) - (h_3(i,j)/2); % Root radius
ratiort_3(i,j) = rt_3(i,j)/rr_3(i,j); % Ratio of Tip to Root radius

ratior_3 (i,j)= r_m(i)/rr_3(i,j) ; % Mean Radius to Root radius ratio
ratiot_3(i,j) = r_m(i)/rt_3(i,j) ; % Mean Radius to Tip radius ratio


%% RADIAL VARIATION OF ANGLES (FREE VORTEX THEORY)

% Station 2
alpha_2r(i,j) = atand(ratior_2(i,j)*tand(alpha_2m(i,j))); % Absolute Angle at Blade Root
alpha_2t(i,j) = atand(ratiot_2(i,j)*tand(alpha_2m(i,j))); % Absolute Angle at Blade Tip

beta_2r(i,j) = atand((ratior_2(i,j)*tand(alpha_2m(i,j))) - (1/ratior_2(i,j))*U(i)/C_a2(i,j)); % Absolute Angle at Blade Root
beta_2t(i,j) = atand((ratiot_2(i,j)*tand(alpha_2m(i,j))) - (1/ratiot_2(i,j))*U(i)/C_a2(i,j)); % Absolute Angle at Blade Tip


% Station 3
alpha_3r(i,j) = atand(ratior_3(i,j)*tand(alpha_3m)); % Absolute Angle at Blade Root
alpha_3t(i,j) = atand(ratiot_3(i,j)*tand(alpha_3m)); % Absolute Angle at Blade Tip

beta_3r(i,j) = atand((ratior_3(i,j)*tand(alpha_3m)) + (1/ratior_3(i,j))*U(i)/C_a3(i,j)); % Relative Angle at Blade Root
beta_3t(i,j) = atand((ratiot_3(i,j)*tand(alpha_3m)) + (1/ratiot_3(i,j))*U(i)/C_a3(i,j)); % Relative Angle at Blade Tip


U_1r(i,j) = (2*pi*N * rr_1(i,j)); % Blade Speed at root at Station 1
U_1t(i,j) = (2*pi*N * rt_1(i,j)); % Blade Speed at tip at Station 1

U_2r(i,j) = (2*pi*N * rr_2(i,j)); % Blade Speed at root at Station 2
U_2t(i,j) = (2*pi*N * rt_2(i,j)); % Blade Speed at tip at Station 2

U_3r(i,j) = (2*pi*N * rr_3(i,j)); % Blade Speed at root at Station 3
U_3t(i,j) = (2*pi*N * rt_3(i,j)); % Blade Speed at tip at Station 3



phi_r(i,j) = 1/(tand(beta_3r(i,j)) - tand(alpha_3r(i,j))); % Flow coefficient at Root
Lambda_r(i,j) = (phi_r(i,j)/2)*(tand(beta_3r(i,j)) - tand(beta_2r(i,j))); % Stage loading at Root
V_2r(i,j) = C_a2(i,j)/cosd(beta_2r(i,j)); % Relative Velocity at Blade Root
C_2r(i,j) = C_a2(i,j)/cosd(alpha_2r(i,j)); % Absolute Velocity at Blade Tip
T_2r = T_o1 - (C_2r(i,j)^2/(2*c_pg)); % Static Temperature at blade root 
M_2ar(i,j) = C_a2(i,j)/sqrt(gamma_g * R * T_2r);

V_2t(i,j) = C_a2(i,j)/cosd(beta_2t(i,j)); % Relative Velocity at Blade Root
C_2t(i,j) = C_a2(i,j)/cosd(alpha_2t(i,j)); % Absolute Velocity at Blade Tip



phi_t(i,j) = 1/(tand(beta_3t(i,j)) - tand(alpha_3t(i,j))); % Flow coefficient at Root
Lambda_t(i,j) = (phi_t(i,j)/2)*(tand(beta_3t(i,j)) - tand(beta_2t(i,j))); % Reaction at Root
V_3r(i,j) = C_a3(i,j)/cosd(beta_3r(i,j)); % Relative Velocity at Blade Root
C_3r(i,j) = C_a3(i,j)/cosd(alpha_3r(i,j)); % Absolute Velocity at Blade Tip

V_3t(i,j) = C_a3(i,j)/cosd(beta_3t(i,j)); % Relative Velocity at Blade Root
C_3t(i,j) = C_a3(i,j)/cosd(alpha_3t(i,j)); % Absolute Velocity at Blade Tip


end

end

figure
nexttile
plot(U, Lambda, '-x')
xlabel('Mean Blade Speed, U (m/s)')
ylabel('Degree of Reaction, \Lambda')
yline(0.5, '-.', 'LineWidth',1.2);
legendStrings = "\phi = " + string(phi);
legend(legendStrings)
hold off


nexttile
plot(U, h_2, '-x')
xlabel('Mean Blade Speed, U (m/s)')
ylabel('Blade Heeight at Rotor Inlet, h_2 (m)')
yline(0.04, '-.', 'LineWidth',1.2);
legendStrings = "\phi = " + string(phi);
legend(legendStrings)
hold off


nexttile
plot(U, M_2ar, '-x')
xlabel('Mean Blade Speed, U (m/s)')
ylabel('Mach Number (Absolute) rotor inlet, M_{2a}')
yline(0.4, '-.', 'LineWidth',1.2);
yline(0.32, '-.', 'LineWidth',1.2);
legendStrings = "\phi = " + string(phi);
legend(legendStrings)
hold off


nexttile
plot(U, psi, '-x')
xlabel('Mean Blade Speed, U (m/s)')
ylabel('Stage Loading, \psi')
yline(2.45, '-.', 'LineWidth',1.2);
yline(1.1, '-.', 'LineWidth',1.2);
legendStrings = "\psi vs U";
legend(legendStrings)
hold off


nexttile
plot(U, Lambda_r, '-x')
xlabel('Mean Blade Speed, U (m/s)')
ylabel('Reaction of Stage at Root, \Lambda_r')
legendStrings = "\phi = " + string(phi);
legend(legendStrings)
hold off


nexttile
plot(U, r_m, '-x')
xlabel('Mean Blade Speed, U (m/s)')
ylabel('Mid-Height radius, r_m')
yline(0.3, '-.', 'LineWidth',1.2);
yline(0.4, '-.', 'LineWidth',1.2);
legendStrings = "r_m vs U";
legend(legendStrings)
hold off


% EXTRACT VALUES FOR PHI AND U FROM FIRST STAGE
phi_ext = input('Choose a value for the Flow Coefficient \Phi: ');
    fprintf('\n');
    % Check validity of the option entered by the user
    while ( not(isnumeric(phi_ext)) || isempty(phi_ext))
        fprintf('You have entered a blank or a non-numeric value. Please, try again.\n');
        phi_ext = input('Choose an option: ');
        fprintf('\n');
    end %end while

U_ext = input('Choose a value for the Mean Blade Speed U: ');
    fprintf('\n');
    % Check validity of the option entered by the user
    while ( not(isnumeric(U_ext)) || isempty(U_ext))
        fprintf('You have entered a blank or a non-numeric value. Please, try again.\n');
        U_ext = input('Choose an option: ');
        fprintf('\n');
    end %end while
    
    
f = find(phi==phi_ext);
v = find(U==U_ext);

C_a_s1 = [C_a1(v,f) C_a2(v,f) C_a3(v,f)];
psi_s1 = psi(v);

alpha_m_s1 = [alpha_2m(v,f) alpha_3m]
alpha_r_s1 = [alpha_2r(v,f) alpha_3r(v,f)];
alpha_t_s1 = [alpha_2t(v,f) alpha_3t(v,f)];

beta_m_s1 = [beta_2m(v,f) beta_3m(f)];
beta_r_s1 = [beta_2r(v,f) beta_3r(v,f)];
beta_t_s1 = [beta_2t(v,f) beta_3t(v,f)];

Ur_s1 = [U_2r(v,f) U_3r(v,f)];
Ut_s1 = [U_2t(v,f) U_3t(v,f)];

rr_s1 = [rr_1(v,f) rr_2(v,f) rr_3(v,f)];
rt_s1 = [rt_1(v,f) rt_2(v,f) rt_3(v,f)];

C_m_s1 = [C_1(v,f) C_2(v,f) C_3(v,f)];
C_r_s1 = [C_2r(v,f) C_3r(v,f)];
C_t_s1 = [C_2t(v,f) C_3t(v,f)];

V_m_s1 = [V_2(v,f) V_3(v,f)];
V_r_s1 = [V_2r(v,f) V_3r(v,f)];
V_t_s1 = [V_2t(v,f) V_3t(v,f)];

phi_s1 = [phi_r(v,f) phi_t(v,f)];
Lambda_s1 = [Lambda(v,f) Lambda_r(v,f) Lambda_t(v,f)]

h_s1 = [h_1(v,f) h_2(v,f) h_3(v,f)];
M_2a_s1 = M_2ar(v,f);

h_to_t = rr_1(v,f)/rt_1(v,f);
r_m_s1 = r_m(v);


end


function [C_a,psi,alpha_m,alpha_r,alpha_t,beta_m,beta_r,beta_t,U_r,U_t,rr...
    ,rt,C_m,C_r,C_t,V_m,V_r,V_t,phi,DOR,h,M_2ar]...
    = turbine(m_in, T_drop, T_o1, T_o3, alpha_3m, phi, U, N, et_h ,p_o1,...
    p_o3, alpha_1m, beta_1m, C_1, C_1r,...
    C_1t, V_1m, V_1r, V_1t)

c_pg = 1148; % c_p of gases after combustion (J/kgK)
gamma_g = 1.333;
R = 287; 
%% ANGLES DEFINITION 
psi = (c_pg*T_drop)/U^2; % Stage loading

beta_3m = atand((1/phi) + tand(alpha_3m)); % Calculate relative angle at rotor exit
Lambda = 0.5*((2*phi*tand(beta_3m)) - psi); % Reaction at mid-blade height

beta_2m = atand(tand(beta_3m) - ((2*Lambda)/phi));
alpha_2m = atand((1/phi) + tand(beta_2m));

%% VELOCITY TRIANGLES
C_a2 = phi*U; % Axial Velocity at Interstage (2)
C_2 = C_a2/cosd(alpha_2m); % Absolute Velocity at Interstage (2)
T_2 = T_o1 - C_2^2/(2*c_pg); % Static Temperature at NGV Exit

NGV_pratio = (T_o1/T_2)^(gamma_g / (et_h*(gamma_g - 1)));

r_crit = ((gamma_g+1)/2)^(gamma_g/(et_h*(gamma_g-1)));

if NGV_pratio > r_crit
    % Check for choked flow
    p_2 = p_o1/r_crit;
else 
    p_2 = p_o1/NGV_pratio;
end

A_2 = ((R*T_2)/(p_2*1000)) * m_in/C_a2;
C_a3 = C_a2; % Axial Velocity at Rotor Exit equals the value at Rotor Inlet
C_3 = C_1;
C_a1 = C_a2; % Keep C_a constant throughout stages

% At turbine Inlet
T_1 = T_o1 - C_1^2/(2*c_pg);
p_1 = p_o1 * (T_1/T_o1)^(gamma_g/(gamma_g - 1));
A_1 = ((R*T_1)/(p_1*1000)) * m_in/C_a1;


% At Stage Outlet
T_3 = T_o3 - C_3^2/(2*c_pg);
p_3 = p_o3*(T_3/T_o3)^(gamma_g/(gamma_g-1));
A_3 = ((R*T_3)/(p_3*1000)) * m_in/C_a3;


% Relative Velocities at mid-height
V_2 = C_a2/cosd(beta_2m);
V_3 = C_a3/cosd(beta_3m);

V_m = [V_1m V_2 V_3];

%% BLADE STUDY
   
r_m = U/(2*pi*N) % Mean Radius

A = [A_1 A_2 A_3];

for n = 1:3
h(n) = A(n)*(N/U); % Blade height
rt(n) = r_m + (h(n)/2); % Tip radius
rr(n) = r_m - (h(n)/2); % Root radius
ratiort(n) = rt(n)/rr(n); % Ratio of Tip to Root radius

ratior(n) = r_m/rr(n) ; % Mean Radius to Root radius ratio
ratiot(n) = r_m/rt(n); % Mean Radius to Tip radius ratio
end


%% RADIAL VARIATION OF ANGLES (FREE VORTEX THEORY)

alpha_m = [alpha_1m alpha_2m alpha_3m];
beta_m = [beta_1m beta_2m beta_3m];
C_a = [C_a1 C_a2 C_a3];
sign = [1 1 -1];
for k = 1:3
alpha_r(k)= atand(ratior(k)*tand(alpha_m(k))); % Absolute Angle at Blade Root
alpha_t(k) = atand(ratiot(k)*tand(alpha_m(k))); % Absolute Angle at Blade Tip

beta_r(k) = atand((ratior(k)*tand(alpha_m(k))) - sign(k)*(1/ratior(k))*U/C_a(k)); % Absolute Angle at Blade Root
beta_t(k) = atand((ratiot(k)*tand(alpha_m(k))) - sign(k)*(1/ratiot(k))*U/C_a(k)); % Absolute Angle at Blade Tip

end

U_1r = (2*pi*N * rr(1)); % Blade Speed at root at Station 1
U_1t = (2*pi*N * rt(1)); % Blade Speed at tip at Station 1


U_2r = (2*pi*N * rr(2)); % Blade Speed at root at Station 2
U_2t = (2*pi*N * rt(2)); % Blade Speed at tip at Station 2

U_3r = (2*pi*N * rr(3)); % Blade Speed at root at Station 3
U_3t = (2*pi*N * rt(3)); % Blade Speed at tip at Station 3


U_r = [U_1r U_2r U_3r]; % Store blade speed at root for all stations
U_t = [U_1t U_2t U_3t]; % Store blade speed at root for all stations
C_m = [C_1 C_2 C_3]; % Store absolute speed at mid-height for all stations


phi_r = 1/(tand(beta_r(3)) - tand(alpha_r(3))); % Flow coefficient at Root

Lambda_r = (phi_r/2)*(tand(beta_r(3)) - tand(beta_r(2))); % Stage loading at Root
V_2r = C_a(2)/cosd(beta_r(2)); % Relative Velocity at Blade Root
C_2r = C_a(2)/cosd(alpha_r(2)); % Absolute Velocity at Blade Tip
T_2r = T_o1 - (C_2r^2/(2*c_pg)); % Static Temperature at blade root 
M_2ar = C_a(2)/sqrt(gamma_g * R * T_2r);

V_2t = C_a(2)/cosd(beta_t(2)); % Relative Velocity at Blade Root
C_2t = C_a(2)/cosd(alpha_t(2)); % Absolute Velocity at Blade Tip

phi_t = 1/(tand(beta_t(3)) - tand(alpha_t(3))); % Flow coefficient at Root

Lambda_t = (phi_t/2)*(tand(beta_t(3)) - tand(beta_t(2))); % Reaction at Root
V_3r = C_a(3)/cosd(beta_r(3)); % Relative Velocity at Blade Root
C_3r = C_a(3)/cosd(alpha_r(3)); % Absolute Velocity at Blade Tip

V_3t = C_a(3)/cosd(beta_t(3)); % Relative Velocity at Blade Root
C_3t = C_a(3)/cosd(alpha_t(3)); % Absolute Velocity at Blade Tip

V_r = [V_1r V_2r V_3r]; % Store relative velocities at root for all stations
V_t = [V_1t V_2t V_3t]; % Store relative velocities at tip for all stations
C_r = [C_1r C_2r C_3r]; % Store absolute velocities at root for all stations
C_t = [C_1t C_2t C_3t]; % Store absolute velocities at tip for all stations

DOR = [Lambda Lambda_r Lambda_t]

end