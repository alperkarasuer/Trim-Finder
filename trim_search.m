function [trimVals, fval] = trim_search(initCond, constr, trimConds)
    % Wind axis to body axis transformation matrix
    Lbw = @(alpha, beta) [cos(alpha), 0, -sin(alpha); 0 1 0; sin(alpha) 0 cos(alpha)]*...
        [cos(-beta), sin(-beta), 0; -sin(-beta), cos(-beta), 0; 0 0 1];

    % NED axis to body axis transformation matrix
    Lbe = @(theta,phi) [cos(theta), 0, -sin(theta);...
        sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta);...
        cos(phi)*sin(theta)*cos(phi) -sin(phi), cos(phi)*cos(theta)];
    
    % Find speeds in body axis
    bodyVel = Lbw(trimConds.alpha, trimConds.beta)*[trimConds.mach*trimConds.speedOfSound; 0; 0];
    [u, v, w] = deal(bodyVel(1), bodyVel(2), bodyVel(3));
    
    % Interpolation of aerodynamic data
    Cd = interp1(initCond.Machpoints, initCond.Cd_data, trimConds.mach);
    Cza = interp1(initCond.Machpoints, initCond.Cza_data, trimConds.mach);
    Czq = interp1(initCond.Machpoints, initCond.Czq_data, trimConds.mach);
    Cma = interp1(initCond.Machpoints, initCond.Cma_data, trimConds.mach);
    Cmq = interp1(initCond.Machpoints, initCond.Cmq_data, trimConds.mach);
    Clp = interp1(initCond.Machpoints, initCond.Clp_data, trimConds.mach);
    Czd = interp1(initCond.Machpoints, initCond.Czd_data, trimConds.mach);
    Cmd = interp1(initCond.Machpoints, initCond.Cmd_data, trimConds.mach);
    Cld = interp1(initCond.Machpoints, initCond.Cld_data, trimConds.mach);
    
    % Functions to calculate aerodynamic coefficients
    Cx = Cd;
    Cy = @(def_dr, r) Cza*trimConds.beta + Czd*def_dr - Czq*r*initCond.d/(2*norm(bodyVel));
    Cz = @(def_de, q) Cza*trimConds.alpha + Czd*def_de + Czq*q*initCond.d/(2*norm(bodyVel));
    Cl = @(def_da, p) Cld*def_da + Clp*p*initCond.d/(2*norm(bodyVel));
    Cm = @(def_de, q) Cma*trimConds.alpha + Cmd*def_de + Cmq*q*initCond.d/(2*norm(bodyVel));
    Cn = @(def_dr, r) -Cma*trimConds.beta - Cmd*def_dr + Cmq*r*initCond.d/(2*norm(bodyVel));
    
    % Dynamic pressure
    Qd = 0.5*trimConds.rho*norm(bodyVel)^2;

    % Functions to create vectors for forces and moments
    forces = @(def_dr, def_de, def_da, p, q, r, theta, phi) ...
        Qd*initCond.A*transpose([Cx, Cy(def_dr, r), Cz(def_de, q)]) + initCond.mass*Lbe(theta,phi)*transpose([initCond.g 0 0]);
    moments = @(def_dr, def_de, def_da, p, q, r) Qd*initCond.A*initCond.d*transpose([Cl(def_da,p) Cm(def_de, q) Cn(def_dr, r)]);
    
    % Functions to create vectors of body velocities, angular rates and roll angle

    velDots = @(def_dr, def_de, def_da, p, q, r, theta, phi) ((1/initCond.mass)*...
        forces(def_dr, def_de, def_da, p, q, r, theta, phi)) - cross([p; q; r], [u; v; w]);

    angRateDots = @(def_dr, def_de, def_da, p, q, r) initCond.I\moments(def_dr, def_de, def_da, p, q, r) - ...
        initCond.I\(cross([p; q; r], initCond.I*[p; q; r]));    
    
    phiDot = @(p, q, r, theta, phi) p + (q*sin(phi) + r*cos(phi))*tan(theta);    
    
    % Helper function to access individual elements of output of anonymous
    % function
    idx = @(expr, index) expr(index);

    % Distribute vectors to individual elements
    vDot = @(def_dr, def_de, def_da, p, q, r, theta, phi) idx(velDots(def_dr, def_de, def_da, p, q, r, theta, phi), 2);
    wDot =  @(def_dr, def_de, def_da, p, q, r, theta, phi) idx(velDots(def_dr, def_de, def_da, p, q, r, theta, phi), 3);
    pDot = @(def_dr, def_de, def_da, p, q, r) idx(angRateDots(def_dr, def_de, def_da, p, q, r), 1);
    qDot = @(def_dr, def_de, def_da, p, q, r) idx(angRateDots(def_dr, def_de, def_da, p, q, r), 2);
    rDot = @(def_dr, def_de, def_da, p, q, r) idx(angRateDots(def_dr, def_de, def_da, p, q, r), 3);

    % Cost function
    % x = [def_dr, def_de, def_da, p, q, r, theta, phi]
    J = @(x) (vDot(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8)))^2 ...
    + (wDot(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8)))^2 ... 
    + (pDot(x(1), x(2), x(3), x(4), x(5), x(6)))^2 + (qDot(x(1), x(2), x(3), x(4), x(5), x(6)))^2 ... 
    + (rDot(x(1), x(2), x(3), x(4), x(5), x(6)))^2 + (phiDot(x(4), x(5), x(6), x(7), x(8)))^2;

    % Upper and lower bounds of deflections
    lb = -[constr.drlim, constr.delim, constr.dalim, constr.pLim, constr.qLim, constr.rLim, constr.thetaLim, constr.phiLim];
    ub = [constr.drlim, constr.delim, constr.dalim, constr.pLim, constr.qLim, constr.rLim, constr.thetaLim, constr.phiLim];
    
    % Initialize optimization from the mean of upper and lower boundaries
    x0 = (lb+ub)/2;
    
    % Linear inequality constraints
    A = [1, 1, 1, 0, 0, 0, 0, 0;...
        -1, 1, 1, 0, 0, 0, 0, 0;...
        -1, -1, 1, 0, 0, 0, 0, 0;...
        1, -1, 1, 0, 0, 0, 0, 0]; % di = x*dr + y*de + z*da relation
    
    A = [A; -A];
    
    b = 30*pi/180*ones(size(A,1),1); % finLim constraints

    % There are no equality constraints and nonlinear constraints
    Aeq = [];
    beq = [];

    % Find minimum of the cost function J given the constraints
    [trimVals, fval] = fmincon(J, x0, A, b, Aeq, beq, lb, ub);
    
end
    
    
