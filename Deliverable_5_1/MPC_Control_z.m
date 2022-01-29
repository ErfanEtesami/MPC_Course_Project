classdef MPC_Control_z < MPC_Control
    properties
        A_bar, B_bar, C_bar % Augmented system for disturbance rejection
        L                   % Estimator gain for disturbance rejection
    end
    
    methods
        function mpc = MPC_Control_z(sys, Ts, H)
            mpc = mpc@MPC_Control(sys, Ts, H);
            
            [mpc.A_bar, mpc.B_bar, mpc.C_bar, mpc.L] = mpc.setup_estimator();
        end
        
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   d_est        - disturbance estimate
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = ceil(H/Ts); % Horizon steps
            
            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.3)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar(1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            
            % Cost matrices
            Q = diag([1, 100]);     % nx = 2
            R = 1;                  % nu = 1
            
            Ts = 1/20; % Sample time
            rocket = Rocket(Ts);
            [xs, us] = rocket.trim();
            
            % Constraints
            % u in U = { u | Mu <= m }
            M = [1; -1]; 
            m = [80; -50] - M*us(3);

            % Compute LQR controller for unconstrained system
            [~, Qf, ~] = dlqr(mpc.A, mpc.B, Q, R);
            
            % Objective and constraints YALMIP
            obj = 0;
            con = [];
            con = con + (X(:, 2) == mpc.A*X(:, 1) + mpc.B*U(:, 1) + mpc.B*d_est);
            con = con + (M*U(:, 1) <= m);
            obj = obj + (U(:, 1)-u_ref)'*R*(U(:, 1)-u_ref); 
            for i = 2:1:N-1
                con = con + (X(:, i+1) == mpc.A*X(:, i) + mpc.B*U(:, i) + mpc.B*d_est);
                con = con + (M*U(:, i) <= m);
                obj = obj + (X(:, i)-x_ref)'*Q*(X(:, i)-x_ref) + (U(:, i)-u_ref)'*R*(U(:, i)-u_ref);
            end
            obj = obj + (X(:, N)-x_ref)'*Qf*(X(:, N)-x_ref);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref, d_est}, U(:,1));
        end
        
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx = size(mpc.A, 1);
            
            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.3)
            ref = sdpvar;
            
            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            
            % Cost matrix
            Rs = 1;
            
            Ts = 1/20; % Sample time
            rocket = Rocket(Ts);
            [xs_l, us_l] = rocket.trim();
            xs_l = xs_l([9, 12]);
            us_l = us_l(3);
            
            % Constraints
            % u in U = { u | Mu <= m }
            M = [1; -1]; 
            m = [80; -50] - M*us_l;
            
            obj = 0;
            con = [];
            con = con + (xs == mpc.A*xs + mpc.B*us + mpc.B*d_est);
            con = con + (M*us <= m);
            con = con + (mpc.C*(xs+xs_l) == ref);
            obj = obj + (us+us_l)'*Rs*(us+us_l);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), {ref, d_est}, {xs, us});
        end
        
        
        % Compute augmented system and estimator gain for input disturbance rejection
        function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)
            
            %%% Design the matrices A_bar, B_bar, L, and C_bar
            %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
            %%% converges to the correct state and constant input disturbance
            %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            
            A_bar = [mpc.A, mpc.B; zeros(1, 2), 1];
            B_bar = [mpc.B; 0];
            C_bar = [mpc.C, 0];
            
            L = -place(A_bar', C_bar', [0.1, 0.2, 0.3])';
            
            Q = diag([1, 100]);
            R = 1;
            [Klqr, ~, ~] = dlqr(mpc.A, mpc.B, Q, R);
            Acl = mpc.A - mpc.B * Klqr;
            disp("Closed-loop poles of the unconstrained infinite-horizon system without observer");
            disp(eig(Acl));
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end