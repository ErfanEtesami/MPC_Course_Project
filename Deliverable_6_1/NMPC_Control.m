function opti_eval = NMPC_Control(rocket, H)

import casadi.*
opti = casadi.Opti(); % Optimization problem

N = ceil(H/rocket.Ts); % MPC horizon
nx = 12; % Number of states
nu = 4;  % Number of inputs

% Decision variables (symbolic)
X_sym = opti.variable(nx, N);       % state trajectory
U_sym = opti.variable(nu, N-1);     % control trajectory)

% Parameters (symbolic)
x0_sym  = opti.parameter(nx, 1);  % initial state
ref_sym = opti.parameter(4, 1);   % target position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

% Decision varaibles for steady state (symbolic)
X_ref = opti.variable(nx, 1);   % state trajectory
U_ref = opti.variable(nu, 1);   % control trajectory

f_discrete = @(x,u) RK4(x, u, rocket);

% Constraint Matrices
% u in U = { u | Mu <= m }
M = [1 0 0 0;
     -1 0 0 0;
     0 1 0 0;
     0 -1 0 0;
     0 0 1 0;
     0 0 -1 0;
     0 0 0 1;
     0 0 0 -1];
m = [0.26;
     0.26; 
     0.26; 
     0.26; 
     80; 
     -50; 
     20; 
     20];
% x in X = { x | Fx <= f }
F = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0]; 
f = [deg2rad(85); deg2rad(85)];

% Cost matrix for steady state
Rs = diag([1, 1, 1, 1]);

obj = 0;

% Constraints and Objective for steady state
opti.subject_to(X_ref == f_discrete(X_ref, U_ref));
opti.subject_to(F*X_ref <= f);
opti.subject_to(M*U_ref <= m);
opti.subject_to(X_ref([10, 11, 12, 6]) == ref_sym);
obj = obj + (U_ref)'*Rs*(U_ref);

% Cost matrices for tracking
Q = diag([100, 100, 1, 1, 1, 1000, 1, 1, 1, 1000, 1000, 1000]);
R = diag([1, 1, 1, 1]);

% Constraints and Objective for tracking
opti.subject_to(X_sym(:, 2) == f_discrete(X_sym(:, 1), U_sym(:, 1)));
opti.subject_to(M*U_sym(:, 1) <= m);
obj = obj + (U_sym(:, 1)-U_ref)'*R*(U_sym(:, 1)-U_ref);
for i = 2:1:N-1
    opti.subject_to(X_sym(:, i+1) == f_discrete(X_sym(:, i), U_sym(:, i)));
    opti.subject_to(F*X_sym(:, i) <= f);
    opti.subject_to(M*U_sym(:, i) <= m);
    obj = obj + (X_sym(:, i)-X_ref)'*Q*(X_sym(:, i)-X_ref) + (U_sym(:, i)-U_ref)'*R*(U_sym(:, i)-U_ref);
end

% Initial conditions
opti.subject_to(X_sym(:, 1) == x0_sym);

% Compute LQR controller for unconstrained system
[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
sys_d = c2d(sys, rocket.Ts);
[A, B, ~, ~] = ssdata(sys_d);
[~, Qf, ~] = dlqr(A, B, Q, R);

% Terminal Cost
obj = obj + (X_sym(:, N)-X_ref)'*Qf*(X_sym(:, N)-X_ref);

% Objective
opti.minimize(obj);

% YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Setup solver ------
ops = struct('ipopt', struct('print_level', 0, 'tol', 1e-3), 'print_time', false);
opti.solver('ipopt', ops);

% Create function to solve and evaluate opti
opti_eval = @(x0_, ref_) solve(x0_, ref_, opti, x0_sym, ref_sym, U_sym);
end

function u = solve(x0, ref, opti, x0_sym, ref_sym, U_sym)

% ---- Set the initial state and reference ----
opti.set_value(x0_sym, x0);
opti.set_value(ref_sym, ref);

% ---- Solve the optimization problem ----
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');

u = opti.value(U_sym(:,1));

% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end