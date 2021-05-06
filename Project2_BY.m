clc
clear
close all

%% Initial Comments 
%{

%=========================================================================================================================================%
AEM-368-001, Spring 2021
Project 2
Kyra Bryan, Kristan Young
%=========================================================================================================================================%

Description:
MATLAB program that calculates the aerodynamic efficiency and longitudinal
stability curves versus angle of attack for a prototype UAV based on the
principles discussed in AEM 368.
Valid only for the input bounds given in the prompt

Nomenclature and variables:

h_g = flight altitude, geometric
M_flight = flight Mach number
c_d0_fuselage = parasitic drag coefficient of the fuselage, relative to the wing area

S_w = wing area
AR_w = wing aspect ratio
taper_w = wing taper ratio
sweep_w = wing quarter-chord sweep
thickness_w = wing thickness ratio
x_cg = x-axis location of the center of gravity over the root chord LE of the wing; measured from the root chord leading edge to the center of gravity, positive aft
z_cg = z-axis location of the center of gravity over the root chord of the wing; measured from the center of gravity to the wing aerodynamic center, positive down
c_M_AC_W = moment of the wing-body aerodynamic center

S_t = tail area
AR_t = tail aspect ratio
taper_t = tail taper ratio
sweep_t = tail quarter-chord sweep
thickness_t = tail thickness ratio
l2_t_b_w = double length of tail over the span of wing; measured from the center of gravity to the tail aerodynamic center, positive aft
z_t_c_r_w = z-axis location of tail center of gravity; measured from the center of gravity to the tail aerodynamic center, positive down
i_t = tail incidence angle, relative to the wing angle of attack, positive LE down
dw_partial = partial derivative of downwash angle over angle of attack 

r_earth = radius of earth
g0 = gravitational constant
T_SSL = temperature at SSL
P_SSL = pressure at SSL
rho_SSL = density at SSL
R = gas constant 
a0 = lapse rate for altitudes under 30000ft 
gamma = heat capacity ratio
Re_transition = we assume transition from laminar to turbulent at this Re value 
a_0_liftslope = initial lift curve slope
eta = percent lift from wing
epsilon0 = downwash
 
h = geopotential altitude 
T_alt_R = temperature for a given flight altitude in Rankine
T_alt_R = temperature for a given flight altitude in Kelvin
mu_SI = dynamic viscosity in SI units at flight altitude
mu_alt = dynamic viscosity in US units at flight altitude
P_alt = pressure for a given flight altitude
rho_alt = density for a given flight altitude
T_alt_F = temperature for a given flight altitude in Fahrenheit 

b_w = wing span
c_r_w = wing root chord
MAC_w = wing mean aerodynamic chord
b_t = tail span
c_r_t = tail root chord
MAC_t = tail mean aerodynamic chord
v_SOS = speed of sound at altitude
v_flight = flight speed
q = dynamic pressure
Re_mac = Reynolds number based on wing MAC
Re_mac_t = Reynolds number based on tail MAC
x_AC = location of aerodynamic center
h_AC = ratio of location of aerodynamic center to wing MAC
x_cg = location of center of gravity 
h_CG = ratio of location of center of gravity to wing MAC
CRminusACratio = ratio of the distance between root chord leading edge and wing aerodynamic center over the root chord
tau = slope parameter
e_wing = spanwise efficiency factor
e_oswald = oswald efficiency factor 

a_w_rad = wing lift slope in radians
a_t_rad = tail lift slope in radians
a_w = wing lift slope in 1/deg
a_t = tail lift slope in 1/deg
V_HT = horizontal tail volume ratio
alpha = angle of attacks (assume both wing and absolute)
c_m_cg = coefficient of moment about the center of gravity 
alpha_trim = trim angle of attack
SM = static margin
c_l_trim = coefficient of lift at trim angle of attack
W_trim = weight at trim angle of attack
wingloading_trim = wing loading at trim angle of attack
c_l_w = coefficient of lift of wing
c_l_t = coefficient of lift of tail
c_l_total = total coefficient of lift
formfactor_w = form factor for wing
formfactor_t = form factor for tail
C_f_w = coefficient of form drag for wing
C_f_t = coefficient of form drag for tail
Swet_w = wetted area of wing
D_f_w = form drag from wing
Swet_t = wetted area of tail
D_f_t = form drag from tail
D_total_parasitic = total form drag
c_d_parasitic = total coeffcient of parasitic drag
c_d_induced = induced drag coefficient
c_d_total = total coeffcient of drag
c_d_induced_trim = induced drag coefficient at trim angle of attack
c_d_total_trim = total coeffcient of drag at trim angle of attack
%}

%% Read Input File, Check If Within Bounds

% Read Input Files
data_array1 = load('specs.txt'); % reads (loads) the data as a column array
h_g = data_array1(1); %ft, sets var to the value on the first line
M_flight = data_array1(2); %sets var to the value on the second line
c_d0_fuselage = data_array1(3); %sets var to the value on the third line

data_array2 = load('wing.txt'); 
S_w = data_array2(1); %ft^2
AR_w = data_array2(2); 
taper_w = data_array2(3); 
sweep_w = data_array2(4); %deg
thickness_w = data_array2(5); 
x_cg = data_array2(6); 
z_cg = data_array2(7); 
c_M_AC_W = data_array2(8); %lb*ft

data_array3 = load('tail.txt'); 
S_t = data_array3(1); %ft^2
AR_t = data_array3(2); 
taper_t = data_array3(3); 
sweep_t = data_array3(4); %deg
thickness_t = data_array3(5); 
l2_t_b_w = data_array3(6);
z_t_c_r_w = data_array3(7); 
i_t = data_array3(8); %deg
partial_dw = data_array3(9); 

% Check If Within Bounds, and Modify as Needed 
if h_g < 0
    h_g = 0;
elseif h_g > 30000
    h_g = 30000;
end 
if M_flight < 0.05
    M_flight = 0.05;
elseif M_flight > 0.4
    M_flight = 0.4;
end 
if c_d0_fuselage < 0.01
    c_d0_fuselage = 0.01;
elseif c_d0_fuselage > 0.05
    c_d0_fuselage = 0.05;
end  

if S_w < 25
    S_w = 25;
elseif S_w > 250
    S_w = 250;
end 
if AR_w < 4
    AR_w = 4;
elseif AR_w > 16
    AR_w = 16;
end 
if taper_w < 0.2
    taper_w = 0.2;
elseif taper_w > 1
    taper_w = 1;
end 
if sweep_w < 0
    sweep_w = 0;
elseif sweep_w > 20
    sweep_w = 20;
end 
if thickness_w < 0.05
    thickness_w = 0.05;
elseif thickness_w > 0.15
    thickness_w = 0.15;
end 
if x_cg < -2
    x_cg = -2;
elseif x_cg > 3
    x_cg = 3;
end 
if z_cg < -1
    z_cg = -1;
elseif z_cg > 1
    z_cg = 1;
end 
if c_M_AC_W < -0.2
    c_M_AC_W = -0.2;
elseif c_M_AC_W > 0
    c_M_AC_W = 0;
end 

if S_t < (0.1*S_w)
    S_t = (0.1*S_w);
elseif S_t > (0.4*S_w)
    S_t = (0.4*S_w);
end 
if AR_t < 3
    AR_t = 3;
elseif AR_t > 6
    AR_t = 6;
end 
if taper_t < 0.4
    taper_t = 0.4;
elseif taper_t > 1
    taper_t = 1;
end 
if sweep_t < 0
    sweep_t = 0;
elseif sweep_t > 10
    sweep_t = 10;
end 
if thickness_t < 0.05
    thickness_t = 0.05;
elseif thickness_t > 0.15
    thickness_t = 0.15;
end 
if l2_t_b_w < 0.5
    l2_t_b_w = 0.5;
elseif l2_t_b_w > 2
    l2_t_b_w = 2;
end 
if z_t_c_r_w < -1
    z_t_c_r_w = -1;
elseif z_t_c_r_w > 0.5
    z_t_c_r_w = 0.5;
end 
if i_t < -5
    i_t = -5;
elseif i_t > 5
    i_t = 5;
end 
if partial_dw < 0
    partial_dw = 0;
elseif partial_dw > 1
    partial_dw = 1;
end 

%% Other Inputs and Constants

r_earth = 2.09e7; %ft
g0 = 32.2; %ft/s^2
T_SSL = 518.67; %R
P_SSL = 2116.8; %lb/ft^2
rho_SSL = 0.002377; %slug/ft^3
R = 1716; %ft*lb/(slug*R)
a0 = -0.00357; %R/ft
gamma = 1.4;
Re_transition = 500000; 

a_0_liftslope = 2*pi;
eta = 1; %L7.4.2 example 
epsilon0 = 0; 

%% Standard Atmosphere

h = (h_g*r_earth)/(h_g+r_earth); %ft
T_alt_R = T_SSL + a0*h; %R
T_alt_K = T_alt_R*(5/9); %K
mu_SI = 1.458e-6*((T_alt_K^1.5)/(T_alt_K+110.4)); %kg/(m*s), at altitude
mu_alt = (mu_SI*0.0685)/(3.281); %slug/(ft*s)
P_alt = P_SSL*(T_alt_R/T_SSL)^(-g0/(a0*R)); %lb/ft^2
rho_alt = P_alt/(T_alt_R*R); %slug/ft^3
T_alt_F = T_alt_R - 459.67; %F

%% Aerodynamics

b_w = sqrt(S_w*AR_w);  %ft
c_r_w = (2*S_w)/(b_w*(1+taper_w)); %ft
MAC_w = ((2/3)*c_r_w)*((taper_w^2+taper_w+1)/(taper_w+1)); %ft

b_t = sqrt(S_t*AR_t);  %ft
c_r_t = (2*S_t)/(b_t*(1+taper_t)); %ft
MAC_t = ((2/3)*c_r_t)*((taper_t^2+taper_t+1)/(taper_t+1)); %ft

v_SOS = sqrt(gamma*T_alt_R*R); %ft/s
v_flight = M_flight*v_SOS; %ft/s
q = 0.5*rho_alt*v_flight^2; %slug/(ft*s^2)
Re_mac = (rho_alt*MAC_w*v_flight)/(mu_alt);
Re_mac_t = (rho_alt*MAC_t*v_flight)/(mu_alt);

% Start with leading edge

x_AC = c_r_w/4 + MAC_w*tand(sweep_w); %ft
h_AC = x_AC/MAC_w; 
x_cg = x_cg*c_r_w; %ft
h_CG = x_cg/MAC_w; 
CRminusACratio = abs(x_cg-x_AC)/c_r_w;

%to interpolate tau and e_wing
AR_sample = [16 12 8 4 2;
        16 12 8 4 2;
        16 12 8 4 2;
        16 12 8 4 2;
        16 12 8 4 2;
        16 12 8 4 2];
taper_sample = [1 1 1 1 1;
    0.8 0.8 0.8 0.8 0.8;
    0.6 0.6 0.6 0.6 0.6;
    0.4 0.4 0.4 0.4 0.4;
    0.2 0.2 0.2 0.2 0.2;
    0 0 0 0 0];
tau_sample = [0.2956 0.2511 0.1952 0.1195 0.0682;
              0.2252 0.1908 0.1478 0.0901 0.0513;
              0.1527 0.1288 0.0991 0.0598 0.0338;
              0.0933 0.0789 0.0610 0.0369 0.0208;
              0.0939 0.0846 0.0715 0.0500 0.0318;
              0.3236 0.3091 0.2850 0.2343 0.1778];
e_wing_sample = [0.882 0.907 0.937 0.972 0.990;
                0.909 .929 0.952 0.979 0.993;
                0.938 0.953 0.969 0.987 0.995;
                0.960 0.970 0.980 0.992 0.997;
                0.942 0.950 0.962 0.980 0.991;
                0.783 0.797 0.820 0.865 0.913];
tau = interp2(AR_sample,taper_sample,tau_sample,AR_w,taper_w);
e_wing = interp2(AR_sample,taper_sample,e_wing_sample,AR_w,taper_w);

e_oswald = e_wing*0.75;

%% Stability Theory and More Aerodyanmics 

a_w_rad = (a_0_liftslope)/(1+((a_0_liftslope*(1+tau))/(pi*AR_w))); %L2.4.2, L7.4.2
a_t_rad = (a_0_liftslope)/(1+((a_0_liftslope*(1+tau))/(pi*AR_t))); %L2.4.2, L7.4.2
a_w = a_w_rad/(180/pi);
a_t = a_t_rad/(180/pi);

V_HT = (S_t*((l2_t_b_w*b_w)/2))/(S_w*MAC_w); %L7.4.2

%some values are arrays based off of AoA
alpha = linspace(-6,10); %deg
c_m_cg = c_M_AC_W + a_w*alpha*(h_CG-h_AC-(a_t/a_w)*eta*V_HT*(1-partial_dw)) + a_t*eta*V_HT*(i_t+epsilon0); %L7.4.1

p1 = polyfit(c_m_cg,alpha,1);
alpha_trim = polyval(p1,0); %deg

SM = h_AC-h_CG+((a_t/a_w)*eta*V_HT*(1-partial_dw)); %L7.4.2, p4; Should be between 0.05 and 0.25
c_l_trim = a_w*alpha_trim; %L7.4.2 example
W_trim = c_l_trim*0.5*rho_alt*v_flight^2*S_w; %lb
wingloading_trim = W_trim/S_w; %lb/ft^s

%lift calcs
c_l_w = a_w*alpha; %L7.4.2 example
c_l_t = a_t*(alpha-i_t-(partial_dw*alpha));
c_l_total = c_l_w + c_l_t*eta*(S_t/S_w);

%parasitic drag calcs
formfactor_w = 5.46*(thickness_w)^2 + (1.55-sind(sweep_w))*(thickness_w) + 1;
formfactor_t = 5.46*(thickness_t)^2 + (1.55-sind(sweep_t))*(thickness_t) + 1;
if Re_mac < Re_transition
    C_f_w = 1.33/sqrt(Re_mac);
    C_f_t = 1.33/sqrt(Re_mac_t);
else 
    C_f_w = 0.455/((log(Re_mac))^2.58); % Assume totally turbulent, transitional
    C_f_t = 0.455/((log(Re_mac_t))^2.58); % Assume totally turbulent, transitional
end 
Swet_w = S_w*2*(1+(0.2*thickness_w));
D_f_w = C_f_w*formfactor_w*q*Swet_w; %lb
Swet_t = S_t*2*(1+(0.2*thickness_t));
D_f_t = C_f_t*formfactor_t*q*Swet_t;%lb
D_total_parasitic = D_f_w + D_f_t;%lb
c_d_parasitic = (D_total_parasitic/(q*S_w)) + c_d0_fuselage;

%other drag
c_d_induced = ((c_l_total).^2/(pi*AR_w*e_oswald));
c_d_total = c_d_induced + c_d_parasitic;

c_d_induced_trim = ((c_l_trim).^2/(pi*AR_w*e_oswald));
c_d_total_trim = c_d_induced_trim + c_d_parasitic;

%% Output

fprintf('Temperature at an altitude of %0.0fft = %0.2f °F\n',h_g, T_alt_F);
fprintf('Pressure at an altitude of %0.0fft = %0.2f psf\n',h_g, P_alt);
fprintf('Density at an altitude of %0.0fft = %0.6f sl/ft^3\n',h_g, rho_alt);
fprintf('Flight speed = %0.2f ft/s\n',v_flight);
fprintf('Reynolds number based on wing MAC = %0.2f\n',Re_mac);
fprintf('Ratio of wing MAC over root chord = %0.4f\n',MAC_w/c_r_w);
fprintf('Ratio of the distance between root chord LE and wing AC over the root chord = %0.4f\n',CRminusACratio);
fprintf('Wing lift curve slope = %0.2f 1/rad\n',a_w_rad);
fprintf('Tail lift curve slope = %0.2f 1/rad\n',a_t_rad);
fprintf('Horizontal tail volume ratio from MAC = %0.4f\n',V_HT);
fprintf('Trim angle of attack = %0.2f deg\n',alpha_trim);
fprintf('Static margin relative to the MAC = %0.2f\n',SM);
fprintf('Lift coefficient at trim angle of attack = %0.2f\n',c_l_trim);
fprintf('Aircraft weight at trim angle of attack = %0.2f lb\n',W_trim);
fprintf('Aircraft wing loading at trim angle of attack = %0.2f lb/ft^s\n',wingloading_trim);

%% Plotting, Final Output

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]) 

% Plots the wing and tail as line segments on a graph
x1 = 0; % Plots leading edge of wing
y1 = 0;
x2 = b_w/2*sweep_w/100;
y2 = b_w/2;
z1 = [x1 x2];
z1y = [y1 y2];
y3 = -b_w/2;
z2 = [x1 x2];
z2y = [y1 y2];
x4 = x2+ c_r_w*taper_w; % Plots tips of edges of wing
z3 = [x2 x4];
z3y = [y2 y2];
z4 = [x2 x4];
z4y = [y3 y3];
x5 = c_r_w; % Plots trailing edges of wing
z5 = [x5 x4];
z5y = [y1 y2];
z6 = [x5 x4];
z6y = [y1 y3];
x6 = x_cg*c_r_w;
x7 = x6 + i_t*(b_w/2);
x8 = x7 +(b_t/2)*(sweep_t/100); % Plots leading edge of tail
y4 = b_t/2;
y5 = -b_t/2;
z7 = [x7 x8];
z7y = [y1 y4];
z8 = [x7 x8];
z8y = [y1 y5];
x9 = x8 + c_r_w*taper_w; % Plots tip edges of tail
z9 = [x8 x9];
z9y = [y4 y4];
z10 = [x8 x9];
z10y = [y5 y5];
x10 = x8+ c_r_t;
z11 = [x10 x9]; % Plots trailing edges of tail
z11y = [y1 y4];
z12 = [x10 x9];
z12y = [y1 y5];

line(z1,z1y)
hold on
line(z2,z2y)
line(z3,z3y)
line(z4,z4y)
line(z5,z5y)
line(z6,z6y)
line(z7,z7y)
line(z8,z8y)
line(z9,z9y)
line(z10,z10y)
line(z11,z11y)
line(z12,z12y)
line(z1,z6y)
grid on
title('Figure 1: Plot of Wing and Tail View');
xlabel('x - Direction (ft)');
ylabel('y - Direction (ft)');
axis equal
hold off

figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2])
plot(c_d_total,c_l_total)
hold on
grid on
title("Figure 2: Lift-drag polar for wing-tail configuration");
xlabel("Coefficient of drag, C_D");
ylabel("Coefficient of lift, C_L");
hold off

figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2])
plot(alpha,c_m_cg/MAC_w)
hold on
grid on
title("Figure 3: Coefficient of moment about the CG vs. angle of attack for wing-tail configuration");
xlabel("Angle of attack (deg)");
ylabel("Coefficient of moment about the CG with MAC, C_M_,_C_G");
hold off

figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2])
plot(alpha,c_l_total./c_d_total,'DisplayName','Lift-Drag Curve')
hold on
grid on
scatter(alpha_trim,c_l_trim/c_d_total_trim,75,'filled','DisplayName','Point of trim')
legend('Location','bestoutside')
title("Figure 4: Lift/drag ratio vs. angle of attack for wing-tail configuration");
xlabel("Angle of attack (deg)");
ylabel("Lift/drag ratio, c_L/c_D");
hold off

fprintf('\nThe four requested plots have been generated.\n'); 
fprintf('AEM368 Project 2 script complete. --------------------------------------------------\n\n');