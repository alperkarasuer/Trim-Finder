clear
close
clc 

%% Preliminaries and Inputs

% If the data file doesn't exist then
% run the script, load the resulting data for initial values as a struct
if(isfile('data.mat'))
    inits = load('data.mat');
else
    run('data.m')
    inits = load('data.mat');
end

% Input trim conditions (leave empty for defaults)
trimConds.mach = input('Mach Number: ');
trimConds.alpha = deg2rad(input('Angle of Attack: '));
trimConds.beta = deg2rad(input('Sideslip Angle: '));
trimConds.altitude = input('Altitude: ');

% Defaults --> M = 0.9, alpha = 5 deg, beta = 5 deg, h = 3000 m
if isempty(trimConds.mach)
    trimConds.mach = 0.9;
end
if isempty(trimConds.alpha)
    trimConds.alpha = deg2rad(5);
end
if isempty(trimConds.beta)  
    trimConds.beta = deg2rad(5);
end
if isempty(trimConds.altitude)
    trimConds.altitude = 3000;
end

% Physical constraints
physConstr.finLim = 30*pi/180;
physConstr.delim = 30*pi/180;
physConstr.drlim = 30*pi/180;
physConstr.dalim = 30*pi/180;
physConstr.pLim = deg2rad(0.001);
physConstr.qLim = deg2rad(200);
physConstr.rLim = deg2rad(200);
physConstr.phiLim = deg2rad(0.01);
physConstr.thetaLim = deg2rad(90);

% Calculate density, temperature and speed of sound at the given altitude
[trimConds.rho, trimConds.T, trimConds.speedOfSound] = altitudeProp(trimConds.altitude, inits);

%% Solve the problem
% Find the trim condition --> trimValues = [def_dr, def_de, def_da, p, q, r, theta, phi]
[trimValues, fval] = trim_search(inits, physConstr, trimConds);

%% Function to find atmospheric properties 
function [rho, T, speedOfSound] = altitudeProp(h, consts)
    if((h) <= 10000)
        rho = consts.rho0*(1 - 0.00002256*(h))^4.256;
        T = consts.T0*(1 - 0.00002256*(h));
    else
        rho = 0.412*exp(-0.000151*(h-10000));
        T = 0.7744*consts.T0;
    end
    speedOfSound = sqrt(consts.k*consts.R*T); 
end

