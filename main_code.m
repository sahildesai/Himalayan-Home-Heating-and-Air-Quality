% To find the temperature dynamics of the house
% Composite wall is considered

disp('***** HOUSE TEMPERATURE AND AIR QUALITY DYNAMICS *****')
% Import the variables from the excel file:

%%  IMPORTING VARIABLES    
% Importing the variables to the MATLAB workspace
% Importing and Converting imported data to individual variables
filename = 'input.xlsx';
len_sim = xlsread(filename,'Sheet1','C1');%;input('How many hours is the simulation run for ? : ');
%% STOVE SIGNAL
ss = xlsread(filename,'D4:DS4');   % ss = stove signal

ss1 = zeros(3600*120,1);
ss1(1:3600) = ss(1);
for i = 2:120      
    ss1((i-1)*3600+1:i*3600) = ss(i);
end
stove_signal = timeseries(ss1);


%% STOVE FIREPOWER
fp = xlsread(filename,'D5:DS5');   % ot = outside temperature

fp1 = zeros(3600*120,1);
fp1(1:3600) = fp(1);
for i = 2:120      
    fp1((i-1)*3600+1:i*3600) = fp(i);
end
fire_power = timeseries(fp1);

%% OUTSIDE TEMPERATURE
ot = xlsread(filename,'D6:DS6');   % ot = outside temperature

ot1 = zeros(3600*120,1);
ot1(1:3600) = ot(1);
for i = 2:120      
    ot1((i-1)*3600+1:i*3600) = ot(i);
end
outside_temp = timeseries(ot1);

%% NUMBER OF AIR CHANGES
ac = xlsread(filename,'D7:DS7');   % ac = air changes

ac1 = zeros(3600*120,1);
ac1(1:3600) = ac(1);
for i = 2:120      
    ac1((i-1)*3600+1:i*3600) = ac(i);
end
air_changes = timeseries(ac1);

%% WALLS
wall_a = xlsread(filename,'K21:AI21');
wall_b = xlsread(filename,'K22:AI22');
wall_c = xlsread(filename,'K23:AI23');
wall_d = xlsread(filename,'K24:AI24');

[wall_a_1, wall_a_2, wall_a_3, wall_a_4, wall_a_5] = ...
    wall_properties(wall_a);
[wall_b_1, wall_b_2, wall_b_3, wall_b_4, wall_b_5] = ...
    wall_properties(wall_b);
[wall_c_1, wall_c_2, wall_c_3, wall_c_4, wall_c_5] = ...
    wall_properties(wall_c);
[wall_d_1, wall_d_2, wall_d_3, wall_d_4, wall_d_5] = ...
    wall_properties(wall_d);

k_wall_a_1 = wall_a_1(1); Cp_wall_a_1 = wall_a_1(2); rho_wall_a_1 = wall_a_1(3); alpha_wall_a_1 = wall_a_1(4); t_wall_a_1 = wall_a_1(5);
k_wall_a_2 = wall_a_2(1); Cp_wall_a_2 = wall_a_2(2); rho_wall_a_2 = wall_a_2(3); alpha_wall_a_2 = wall_a_2(4); t_wall_a_2 = wall_a_2(5);
k_wall_a_3 = wall_a_3(1); Cp_wall_a_3 = wall_a_3(2); rho_wall_a_3 = wall_a_3(3); alpha_wall_a_3 = wall_a_3(4); t_wall_a_3 = wall_a_3(5);
k_wall_a_4 = wall_a_4(1); Cp_wall_a_4 = wall_a_4(2); rho_wall_a_4 = wall_a_4(3); alpha_wall_a_4 = wall_a_4(4); t_wall_a_4 = wall_a_4(5);
k_wall_a_5 = wall_a_5(1); Cp_wall_a_5 = wall_a_5(2); rho_wall_a_5 = wall_a_5(3); alpha_wall_a_5 = wall_a_5(4); t_wall_a_5 = wall_a_5(5);

k_wall_b_1 = wall_b_1(1); Cp_wall_b_1 = wall_b_1(2); rho_wall_b_1 = wall_b_1(3); alpha_wall_b_1 = wall_b_1(4); t_wall_b_1 = wall_b_1(5);
k_wall_b_2 = wall_b_2(1); Cp_wall_b_2 = wall_b_2(2); rho_wall_b_2 = wall_b_2(3); alpha_wall_b_2 = wall_b_2(4); t_wall_b_2 = wall_b_2(5);
k_wall_b_3 = wall_b_3(1); Cp_wall_b_3 = wall_b_3(2); rho_wall_b_3 = wall_b_3(3); alpha_wall_b_3 = wall_b_3(4); t_wall_b_3 = wall_b_3(5);
k_wall_b_4 = wall_b_4(1); Cp_wall_b_4 = wall_b_4(2); rho_wall_b_4 = wall_b_4(3); alpha_wall_b_4 = wall_b_4(4); t_wall_b_4 = wall_b_4(5);
k_wall_b_5 = wall_b_5(1); Cp_wall_b_5 = wall_b_5(2); rho_wall_b_5 = wall_b_5(3); alpha_wall_b_5 = wall_b_5(4); t_wall_b_5 = wall_b_5(5);

k_wall_c_1 = wall_c_1(1); Cp_wall_c_1 = wall_c_1(2); rho_wall_c_1 = wall_c_1(3); alpha_wall_c_1 = wall_c_1(4); t_wall_c_1 = wall_c_1(5);
k_wall_c_2 = wall_c_2(1); Cp_wall_c_2 = wall_c_2(2); rho_wall_c_2 = wall_c_2(3); alpha_wall_c_2 = wall_c_2(4); t_wall_c_2 = wall_c_2(5);
k_wall_c_3 = wall_c_3(1); Cp_wall_c_3 = wall_c_3(2); rho_wall_c_3 = wall_c_3(3); alpha_wall_c_3 = wall_c_3(4); t_wall_c_3 = wall_c_3(5);
k_wall_c_4 = wall_c_4(1); Cp_wall_c_4 = wall_c_4(2); rho_wall_c_4 = wall_c_4(3); alpha_wall_c_4 = wall_c_4(4); t_wall_c_4 = wall_c_4(5);
k_wall_c_5 = wall_c_5(1); Cp_wall_c_5 = wall_c_5(2); rho_wall_c_5 = wall_c_5(3); alpha_wall_c_5 = wall_c_5(4); t_wall_c_5 = wall_c_5(5);

k_wall_d_1 = wall_d_1(1); Cp_wall_d_1 = wall_d_1(2); rho_wall_d_1 = wall_d_1(3); alpha_wall_d_1 = wall_d_1(4); t_wall_d_1 = wall_d_1(5);
k_wall_d_2 = wall_d_2(1); Cp_wall_d_2 = wall_d_2(2); rho_wall_d_2 = wall_d_2(3); alpha_wall_d_2 = wall_d_2(4); t_wall_d_2 = wall_d_2(5);
k_wall_d_3 = wall_d_3(1); Cp_wall_d_3 = wall_d_3(2); rho_wall_d_3 = wall_d_3(3); alpha_wall_d_3 = wall_d_3(4); t_wall_d_3 = wall_d_3(5);
k_wall_d_4 = wall_d_4(1); Cp_wall_d_4 = wall_d_4(2); rho_wall_d_4 = wall_d_4(3); alpha_wall_d_4 = wall_d_4(4); t_wall_d_4 = wall_d_4(5);
k_wall_d_5 = wall_d_5(1); Cp_wall_d_5 = wall_d_5(2); rho_wall_d_5 = wall_d_5(3); alpha_wall_d_5 = wall_d_5(4); t_wall_d_5 = wall_d_5(5);

w_wall = xlsread(filename,'Sheet1','I28'); l_wall = xlsread(filename,'Sheet1','I29'); h_wall = xlsread(filename,'Sheet1','I30');

%% WINDOW

num_window = xlsread(filename,'Sheet1','I33'); w_window = xlsread(filename,'Sheet1','I34'); h_window = xlsread(filename,'Sheet1','I35');
window = xlsread(filename,'N30:R30');
k_window = window(1); Cp_window = window(2); rho_window = window(3); alpha_window = window(4); t_window = window(5); 

ws = xlsread(filename,'D9:DS9');
ws1 = zeros(3600*120,1);
ws1(1:3600) = ws(1);
for i = 2:120      
    ws1((i-1)*3600+1:i*3600) = ws(i);
end
window_signal = timeseries(ws1);

%% DOOR
w_door = xlsread(filename,'Sheet1','I31'); h_door = xlsread(filename,'Sheet1','I32');
door = xlsread(filename,'N29:R29');
k_door = door(1); Cp_door = door(2); rho_door = door(3); alpha_door = door(4); t_door = door(5);

ds = xlsread(filename,'D8:DS8');
ds1 = zeros(3600*120,1);
ds1(1:3600) = ds(1);
for i = 2:120      
    ds1((i-1)*3600+1:i*3600) = ds(i);
end
door_signal = timeseries(ds1);

%% ROOF
roof = xlsread(filename,'K25:AI25');
[roof_1, roof_2, roof_3, roof_4, roof_5] = ...
    wall_properties(roof);

k_roof_1 = roof_1(1); Cp_roof_1 = roof_1(2); rho_roof_1 = roof_1(3);...
    alpha_roof_1 = roof_1(4); t_roof_1 = roof_1(5);

k_roof_2 = roof_2(1); Cp_roof_2 = roof_2(2); rho_roof_2 = roof_2(3); alpha_roof_2 = roof_2(4); t_roof_2 = roof_2(5);
k_roof_3 = roof_3(1); Cp_roof_3 = roof_3(2); rho_roof_3 = roof_3(3); alpha_roof_3 = roof_3(4); t_roof_3 = roof_3(5);
k_roof_4 = roof_4(1); Cp_roof_4 = roof_4(2); rho_roof_4 = roof_4(3); alpha_roof_4 = roof_4(4); t_roof_4 = roof_4(5);
k_roof_5 = roof_5(1); Cp_roof_5 = roof_5(2); rho_roof_5 = roof_5(3); alpha_roof_5 = roof_5(4); t_roof_5 = roof_5(5);

%% ADDITIONAL PROPERTIES
initial_room_temp = xlsread(filename,'Sheet1','X28'); 
Cp_air = xlsread(filename,'Sheet1','X29');
h_inf = xlsread(filename,'Sheet1','X31');
h_room = xlsread(filename,'Sheet1','X32');


amb_CO2 = xlsread(filename,'Sheet1','X32');
amb_CO = xlsread(filename,'Sheet1','X33');
amb_PM = xlsread(filename,'Sheet1','X34');
house_dim = [w_wall l_wall h_wall num_window w_window h_window w_door h_door];
To = initial_room_temp;
floor = xlsread(filename,'N31:R31');
k_floor = floor(1); Cp_floor = floor(2); rho_floor = floor(3); alpha_floor = floor(4); t_floor = floor(5);

% Layer by Layer Properties
hb = 1.2; hinf = 1.2;
roof_c1 = alpha_roof_1/(t_roof_1*t_roof_1);
roof_c2 = alpha_roof_2/(t_roof_2*t_roof_2);
roof_c3 = alpha_roof_3/(t_roof_3*t_roof_3);
roof_c4 = alpha_roof_4/(t_roof_4*t_roof_4);
roof_c5 = alpha_roof_5/(t_roof_5*t_roof_5);
roof_p1 = 4*t_roof_1*hb/k_roof_1;
roof_p2 = 4*t_roof_5*hinf/k_roof_5;

%% RUNNING THE SIMULATION

sim('house_dynamics',len_sim*3600);


%% PLOTTING

% Room Temperature
l_T = length(plot_RT);
t_T = linspace(0,len_sim,l_T);
figure (1)
plot(t_T,plot_RT)
hold on
plot(t_T,plot_OT)
grid on
grid minor
xlim([0 len_sim])
xlabel('Time (hr)')
ylabel('Temperature (deg C)')
title('Room Temperature variation with reference to outside temperature variation')
legend('Rooom Temperature','Outside Temperature')

%Room Air Density
l_rho = length(plot_rho);
t_rho = linspace(0,len_sim,l_rho);
figure (2)
plot(t_rho,plot_rho)
grid on
grid minor
xlim([0 len_sim])
xlabel('Time (hr)')
ylabel('Densiity (kg/m3)')
title('Room Density variation with time')
legend('Room Density')

% Air Quality
l_AQ = length(plot_CO2);
t_AQ = linspace(0,len_sim,l_AQ);
figure (3)
plot(t_AQ,plot_CO2)
hold on
plot(t_AQ,plot_CO)
plot(t_AQ,plot_PM)
grid on
grid minor
xlim([0 len_sim])
xlabel('Time (hr)')
ylabel('Concentration (ppm)')
title('CO_2 CO and PM2.5 concentration variation with time')
legend('CO_2','CO','PM2.5')

% Stove Heat
l_stove = length(plot_stoveheat);
t_stove = linspace(0,len_sim,l_stove);
figure (4)
plot(t_stove,plot_stoveheat/1000)
grid on
grid minor
xlim([0 len_sim])
xlabel('Time (hr)')
ylabel('Stove Heat Output (kW)')
title('Stove Heat Output Variation')
legend('Stove Heat')