clear
close
clc 

% If the data file doesn't exist then
% run the script, load the resulting data for initial values as a struct
if(isfile('data.mat'))
    inits = load('data.mat');
else
    run('data.m')
    inits = load('data.mat');
end


Lbw = @(alpha, beta) [cos(beta)*cos(alpha), sin(beta), -cos(beta)*sin(alpha);...
    -sin(beta)*cos(alpha), cos(beta), sin(alpha)*sin(beta);...
    sin(alpha), 0, cos(alpha)];

Lbe = @(theta,phi) [cos(theta), 0, -sin(theta);...
    sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta);...
    cos(phi)*sin(theta)*cos(phi) -sin(phi), cos(phi)*cos(theta)];

input_trim = [0.9 5 5 3000];

states.mach = input_trim(1);
states.alpha = input_trim(2)*pi/180;
states.beta = input_trim(3)*pi/180;
states.height = input_trim(4);

constr.finLim = 30*pi/180;
constr.delim = 30*pi/180;
constr.drlim = 30*pi/180;
constr.dalim = 30*pi/180;
constr.pLim = deg2rad(0.001);
constr.qLim = deg2rad(200);
constr.rLim = deg2rad(200);
constr.phiLim = deg2rad(0.01);
constr.thetaLim = deg2rad(90);

%%

[altitude.rho, altitude.T, altitude.speedOfSound] = altitudeProp(states.height, inits);

bodyVel = Lbw(states.alpha, states.beta)*[states.mach*altitude.speedOfSound; 0; 0];
states.u = bodyVel(1);
states.v = bodyVel(2);
states.w = bodyVel(3);

%%

Cd = interp1(inits.Machpoints, inits.Cd_data, states.mach);

Cza = interp1(inits.Machpoints, inits.Cza_data, states.mach);

Czq = interp1(inits.Machpoints, inits.Czq_data, states.mach);

Cma = interp1(inits.Machpoints, inits.Cma_data, states.mach);

Cmq = interp1(inits.Machpoints, inits.Cmq_data, states.mach);

Clp = interp1(inits.Machpoints, inits.Clp_data, states.mach);

Czd = interp1(inits.Machpoints, inits.Czd_data, states.mach);

Cmd = interp1(inits.Machpoints, inits.Cmd_data, states.mach);

Cld = interp1(inits.Machpoints, inits.Cld_data, states.mach);

%%

Cx = Cd;

Cy = @(def_dr, r) Cza*states.beta + Czd*def_dr - Czq*r*inits.d/(2*norm(bodyVel));

Cz = @(def_de, q) Cza*states.alpha + Czd*def_de + Czq*q*inits.d/(2*norm(bodyVel));

Cl = @(def_da, p) Cld*def_da + Clp*p*inits.d/(2*norm(bodyVel));

Cm = @(def_de, q) Cma*states.alpha + Cmd*def_de + Cmq*q*inits.d/(2*norm(bodyVel));

Cn = @(def_dr, r) -Cma*states.beta - Cmd*def_dr + Cmq*r*inits.d/(2*norm(bodyVel));

%%

Qd = 0.5*altitudeProp(states.height, inits)*norm(bodyVel)^2;

forces = @(def_dr, def_de, def_da, p, q, r, theta, phi) ...
    Qd*inits.d*transpose([Cx, Cy(def_dr, r), Cz(def_de, q)]) + inits.mass*Lbe(theta,phi)*transpose([inits.g 0 0]);

moments = @(def_dr, def_de, def_da, p, q, r) Qd*inits.A*inits.d*transpose([Cl(def_da,p) Cm(def_de, q) Cn(def_dr, r)]);

%%

idx = @(expr, index) expr(index); % helper func

omega = @(p, q, r) [0 -r q; r 0 -p; -q p 0];

velDots = @(def_dr, def_de, def_da, p, q, r, theta, phi) -omega(p,q,r)*[states.u; states.v; states.w] + ...
    (1/inits.mass)*forces(def_dr, def_de, def_da, p, q, r, theta, phi);

angRateDots = @(def_dr, def_de, def_da, p, q, r) inits.I\moments(def_dr, def_de, def_da, p, q, r) - ...
    inits.I\(cross([p; q; r], inits.I*[p; q; r]));


vDot = @(def_dr, def_de, def_da, p, q, r, theta, phi) idx(velDots(def_dr, def_de, def_da, p, q, r, theta, phi), 2);
wDot =  @(def_dr, def_de, def_da, p, q, r, theta, phi) idx(velDots(def_dr, def_de, def_da, p, q, r, theta, phi), 3);
pDot = @(def_dr, def_de, def_da, p, q, r) idx(angRateDots(def_dr, def_de, def_da, p, q, r), 1);
qDot = @(def_dr, def_de, def_da, p, q, r) idx(angRateDots(def_dr, def_de, def_da, p, q, r), 2);
rDot = @(def_dr, def_de, def_da, p, q, r) idx(angRateDots(def_dr, def_de, def_da, p, q, r), 3);
phiDot = @(p, q, r, theta, phi) p + (q*sin(phi) + r*cos(phi))*tan(theta);
%%


J = @(x) (vDot(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8)))^2 ...
+ (wDot(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8)))^2 ... 
+ (pDot(x(1), x(2), x(3), x(4), x(5), x(6)))^2 + (qDot(x(1), x(2), x(3), x(4), x(5), x(6)))^2 ... 
+ (rDot(x(1), x(2), x(3), x(4), x(5), x(6)))^2 + (phiDot(x(4), x(5), x(6), x(7), x(8)))^2;

% J = @(def_dr, def_de, def_da, p, q, r, theta, phi) (vDot(def_dr, def_de, def_da, p, q, r, theta, phi))^2 ...
% + (wDot(def_dr, def_de, def_da, p, q, r, theta, phi))^2 ... 
% + (pDot(def_dr, def_de, def_da, p, q, r))^2 + (qDot(def_dr, def_de, def_da, p, q, r))^2 ... 
% + (rDot(def_dr, def_de, def_da, p, q, r))^2 + (phiDot(p, q, r, theta, phi))^2;

lb = -[constr.drlim, constr.delim, constr.dalim, constr.pLim, constr.qLim, constr.rLim, constr.thetaLim, constr.phiLim];
ub = [constr.drlim, constr.delim, constr.dalim, constr.pLim, constr.qLim, constr.rLim, constr.thetaLim, constr.phiLim];

x0 = (lb+ub)/2;
A = [1, 1, 1, 0, 0, 0, 0, 0;...
    -1, 1, 1, 0, 0, 0, 0, 0;...
    -1, -1, 1, 0, 0, 0, 0, 0;...
    1, -1, 1, 0, 0, 0, 0, 0];

b = 30*pi/180*ones(size(A,1),1);

Aeq = [];
beq = [];

    
[x, fval] = fmincon(J, x0, A, b, Aeq, beq, lb, ub)
%%
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

