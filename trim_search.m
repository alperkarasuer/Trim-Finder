 function [trimVals, fval, deriv] = trim_search(initConds, constr, trimConds)
    % Wind axis to body axis transformation matrix
    Lbw = @(alpha, beta) [cos(alpha), 0, -sin(alpha); 0 1 0; sin(alpha) 0 cos(alpha)]*...
        [cos(-beta), sin(-beta), 0; -sin(-beta), cos(-beta), 0; 0 0 1];

    % NED axis to body axis transformation matrix
    Lbe = @(phi, theta) [cos(theta), 0, -sin(theta);...
        sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta);...
        cos(phi)*sin(theta)*cos(phi) -sin(phi), cos(phi)*cos(theta)];
    
    % Find speeds in body axis
    bodyVel = Lbw(trimConds.alpha, trimConds.beta)*[trimConds.mach*trimConds.speedOfSound; 0; 0];
    [u, v, w] = deal(bodyVel(1), bodyVel(2), bodyVel(3));
    
    % Interpolation of aerodynamic data
    Cd = interp1(initConds.Machpoints, initConds.Cd_data, trimConds.mach);
    Cza = interp1(initConds.Machpoints, initConds.Cza_data, trimConds.mach);
    Czq = interp1(initConds.Machpoints, initConds.Czq_data, trimConds.mach);
    Cma = interp1(initConds.Machpoints, initConds.Cma_data, trimConds.mach);
    Cmq = interp1(initConds.Machpoints, initConds.Cmq_data, trimConds.mach);
    Clp = interp1(initConds.Machpoints, initConds.Clp_data, trimConds.mach);
    Czd = interp1(initConds.Machpoints, initConds.Czd_data, trimConds.mach);
    Cmd = interp1(initConds.Machpoints, initConds.Cmd_data, trimConds.mach);
    Cld = interp1(initConds.Machpoints, initConds.Cld_data, trimConds.mach);
    
    
    % Functions to calculate aerodynamic coefficients
    Cx = Cd;
    Cy = @(def_dr, r, beta) Cza*beta + Czd*def_dr - Czq*r*initConds.d/(2*norm(bodyVel));
    Cz = @(def_de, q, alpha) Cza*alpha + Czd*def_de + Czq*q*initConds.d/(2*norm(bodyVel));
    Cl = @(def_da, p) Cld*def_da + Clp*p*initConds.d/(2*norm(bodyVel));
    Cm = @(def_de, q, alpha) Cma*alpha + Cmd*def_de + Cmq*q*initConds.d/(2*norm(bodyVel));
    Cn = @(def_dr, r, beta) -Cma*beta - Cmd*def_dr + Cmq*r*initConds.d/(2*norm(bodyVel));
    
    % Dynamic pressure
    Qd = 0.5*trimConds.rho*norm(bodyVel)^2;

    % Functions to create vectors for forces and moments
    forces = @(def_de, def_dr, def_da, p, q, r, phi, theta, alpha, beta) ...
        Qd*initConds.A*transpose([Cx, Cy(def_dr, r, beta), Cz(def_de, q, alpha)]) + initConds.mass*Lbe(phi, theta)*transpose([0 0 initConds.g]);
    
    moments = @(def_de, def_dr, def_da, p, q, r, alpha, beta) Qd*initConds.A*initConds.d*transpose([Cl(def_da,p) Cm(def_de, q, alpha) Cn(def_dr, r, beta)]);
    
    % Functions to create vectors of body velocities, angular rates and roll angle

    velDots = @(def_de, def_dr, def_da, p, q, r, phi, theta, alpha, beta) ((1/initConds.mass)*...
        forces(def_de, def_dr, def_da, p, q, r, phi, theta, alpha, beta)) - cross([p; q; r], [u; v; w]);

    angRateDots = @(def_de, def_dr, def_da, p, q, r, alpha, beta) initConds.I\moments(def_de, def_dr, def_da, p, q, r, alpha, beta) - ...
        initConds.I\(cross([p; q; r], initConds.I*[p; q; r]));    
    
    phiDot = @(p, q, r, phi, theta) p + (q*sin(phi) + r*cos(phi))*tan(theta);    
    
    % Helper function to access individual elements of output of anonymous
    % function
    idx = @(expr, index) expr(index);

    % Distribute vectors to individual elements
    vDot = @(def_de, def_dr, def_da, p, q, r, phi, theta, alpha, beta) idx(velDots(def_de, def_dr, def_da, p, q, r, phi, theta, alpha, beta), 2);
    wDot = @(def_de, def_dr, def_da, p, q, r, phi, theta, alpha, beta) idx(velDots(def_de, def_dr, def_da, p, q, r, phi, theta, alpha, beta), 3);
    pDot = @(def_de, def_dr, def_da, p, q, r, alpha, beta) idx(angRateDots(def_de, def_dr, def_da, p, q, r, alpha, beta), 1);
    qDot = @(def_de, def_dr, def_da, p, q, r, alpha, beta) idx(angRateDots(def_de, def_dr, def_da, p, q, r, alpha, beta), 2);
    rDot = @(def_de, def_dr, def_da, p, q, r, alpha, beta) idx(angRateDots(def_de, def_dr, def_da, p, q, r, alpha, beta), 3);

    % Cost function
    % x = [def_de, def_dr, def_da, p, q, r, phi, theta]
    J = @(x) (vDot(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), trimConds.alpha, trimConds.beta))^2 ...
    + (wDot(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), trimConds.alpha, trimConds.beta))^2 ... 
    + (pDot(x(1), x(2), x(3), x(4), x(5), x(6), trimConds.alpha, trimConds.beta))^2 ...
    + (qDot(x(1), x(2), x(3), x(4), x(5), x(6), trimConds.alpha, trimConds.beta))^2 ... 
    + (rDot(x(1), x(2), x(3), x(4), x(5), x(6), trimConds.alpha, trimConds.beta))^2 ...
    + (phiDot(x(4), x(5), x(6), x(7), x(8)))^2;

    % Upper and lower bounds of deflections
    lb = -[constr.delim, constr.drlim, constr.dalim, constr.pLim, constr.qLim, constr.rLim, constr.phiLim, constr.thetaLim];
    ub = [constr.delim, constr.drlim, constr.dalim, constr.pLim, constr.qLim, constr.rLim, constr.phiLim, constr.thetaLim];
    
    % Initialize optimization from the mean of upper and lower boundaries
    x0 = (lb+ub)/2;
    
    % Linear inequality constraints
    A = [1, 1, 1, 0, 0, 0, 0, 0;...
        1, -1, 1, 0, 0, 0, 0, 0;...
        -1, -1, 1, 0, 0, 0, 0, 0;...
        -1, 1, 1, 0, 0, 0, 0, 0]; % di = x*de + y*dr + z*da relation
    
    A = [A; -A];
    
    b = 30*pi/180*ones(size(A,1),1); % finLim constraints

    % There are no equality constraints and nonlinear constraints
    Aeq = [];
    beq = [];

    % Find minimum of the cost function J given the constraints
    % trimVals = [def_de, def_dr, def_da, p, q, r, phi, theta]
    options = optimoptions('fmincon','Display','notify');
    [trimVals, fval] = fmincon(J, x0, A, b, Aeq, beq, lb, ub, [], options);
    
    [def_de, def_dr, def_da, p, q, r] = deal(trimVals(1), trimVals(2), ...
        trimVals(3), trimVals(4), trimVals(5), trimVals(6));
    
    % Rows are for Cza, Czde, Cma, Cmde, Cyb, Cydr, Clb, Clda, Cnb, Cndr
    % First column is positive, second column is negative perturbation
    deltaCoef = zeros(10,2);
    
    for i = 1:2
        deltaCoef(1,i) = Cz(def_de, q, trimConds.alpha + ((-1)^(i+1))*initConds.delAlpha);
        deltaCoef(2,i) = Cz(def_de + ((-1)^(i+1))*initConds.delDef, q, trimConds.alpha);
        deltaCoef(3,i) = Cm(def_de, q, trimConds.alpha + ((-1)^(i+1))*initConds.delAlpha);
        deltaCoef(4,i) = Cm(def_de + ((-1)^(i+1))*initConds.delDef, q, trimConds.alpha);
        deltaCoef(5,i) = Cy(def_dr, r, trimConds.beta + ((-1)^(i+1))*initConds.delBeta);
        deltaCoef(6,i) = Cy(def_dr + ((-1)^(i+1))*initConds.delDef, q, trimConds.beta);
        deltaCoef(7,i) = Cl(def_da, p);
        deltaCoef(8,i) = Cl(def_da + ((-1)^(i+1))*initConds.delDef, p);
        deltaCoef(9,i) = Cn(def_dr, r, trimConds.beta + ((-1)^(i+1))*initConds.delBeta);
        deltaCoef(10,i) = Cn(def_dr + ((-1)^(i+1))*initConds.delDef, r, trimConds.beta);
    end
    
    % Derivatives
    Cza = (deltaCoef(1,1)-deltaCoef(1,2))/(2*initConds.delAlpha);
    Czde = (deltaCoef(2,1)-deltaCoef(2,2))/(2*initConds.delDef);
    Cma = (deltaCoef(3,1)-deltaCoef(3,2))/(2*initConds.delAlpha);
    Cmde = (deltaCoef(4,1)-deltaCoef(4,2))/(2*initConds.delDef);
    Cyb = (deltaCoef(5,1)-deltaCoef(5,2))/(2*initConds.delBeta);
    Cydr = (deltaCoef(6,1)-deltaCoef(6,2))/(2*initConds.delDef);
    Clb = (deltaCoef(7,1)-deltaCoef(7,2))/(2*initConds.delBeta);
    Clda = (deltaCoef(8,1)-deltaCoef(8,2))/(2*initConds.delAlpha);
    Cnb = (deltaCoef(9,1)-deltaCoef(9,2))/(2*initConds.delBeta);
    Cndr = (deltaCoef(10,1)-deltaCoef(10,2))/(2*initConds.delDef);
    
    deriv = [Cza Czde Cma Cmde Cyb Cydr Clb Clda Cnb Cndr];
    
end
    
    
