classdef MPC_Control_y < MPC_Control
    properties
        flag = 0;
    end
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = ceil(H/Ts); % Horizon steps

            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.2)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            
            % Cost matrices
            Q = diag([1, 1, 1, 100]);       % nx = 4
            R = 1;                          % nu = 1
            
            Ts = 1/20; % Sample time
            rocket = Rocket(Ts);
            [xs, us] = rocket.trim();
            
            % Constraints
            % u in U = { u | Mu <= m }
            M = [1; -1]; 
            m = [0.26; 0.26] - M*us(1);
            % x in X = { x | Fx <= f }
            F = [0, 1, 0, 0; 0, -1, 0, 0]; 
            f = [0.0873; 0.0873] - F*xs([1, 4, 8, 11]);

            % Compute LQR controller for unconstrained system
            [K, Qf, ~] = dlqr(mpc.A, mpc.B, Q, R);
            % MATLAB defines K as -K, so invert its signal
            K = -K; 

            % Compute maximal invariant set
            Xf = polytope([F; M*K], [f; m]);
            Acl = mpc.A + mpc.B*K;
            while 1
                prevXf = Xf;
                [T, t] = double(Xf);
                preXf = polytope(T*Acl, t);
                Xf = intersect(Xf, preXf);
                if isequal(prevXf, Xf)
                    break;
                end
            end
            [Ff, ff] = double(Xf);
            
            if mpc.flag == 0 
                % Plot Xf
                figure;
                for i = 1:1:3
                    subplot(1, 3, i);
                    Xf.projection(i: i+1).plot();
                    title("y controller : Xf for dimension " + num2str(i) + ", " + num2str(i+1));
                end
                mpc.flag = 1;
            end
            
            % Objective and constraints YALMIP
            obj = 0;
            con = [];
            con = con + (X(:, 2) == mpc.A*X(:, 1) + mpc.B*U(:, 1));
            con = con + (M*U(:, 1) <= m);
            obj = obj + U(:, 1)'*R*U(:, 1); 
            for i = 2:1:N-1
                con = con + (X(:, i+1) == mpc.A*X(:, i) + mpc.B*U(:, i));
                con = con + (F*X(:, i) <= f);
                con = con + (M*U(:, i) <= m);
                obj = obj + X(:, i)'*Q*X(:, i) + U(:, i)'*R*U(:, i);
            end
            con = con + (Ff*X(:, N) <= ff);
            obj = obj + X(:, N)'*Qf*X(:, N);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref}, U(:,1));
        end
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx = size(mpc.A, 1);

            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            obj = 0;
            con = [xs == 0, us == 0];
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end
    end
end
