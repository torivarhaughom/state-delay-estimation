clear all
close all
clc

global n_x n_u n_y A B C d1_min d1_max d1_star d2_min d2_max d2_star u1dot_min u1dot_max ... 
    u2dot_min u2dot_max epsilon L_R

% Define system dimensions
n_x = 6;
n_u = 2;
n_y = 2;

% Define system matrices
A = [-2  1  -1  1  0  0;
      1 -3  2  1  0  0;
      1  2 -2  1  3  0;
      0  2  -2 -5  1  0;
      1  3  -2  0 -6  1;
      0  -2  0  1  1 -7];

B = eye(n_x,n_u);
C = eye(n_y,n_x);

% Define delay and input derivative bounds
d1_min = 0.2; d1_max = 0.8; d1_mid = (d1_min + d1_max)/2; d1_star = d1_mid; 
d2_min = 0.2; d2_max = 0.8; d2_mid = (d2_min + d2_max)/2; d2_star = d2_mid; 
u1dot_min = -10; u1dot_max = 10;
u2dot_min = -10; u2dot_max = 10;

% Defining constants
epsilon = 0.01;
lambda_plus = 1.1;
lambda_minus = 1;
mu = 0.5;

alpha = 1;
beta = 50000;

% Switching law constants
lambda = 0.25*lambda_minus;
lambda_star = 0.75*lambda_minus;

% Augmented system (symbolic)
syms u1dot u2dot d1_star_sym d2_star_sym gamma_1 gamma_2 d1_dot d2_dot real 

U = diag([u1dot u2dot]); 
D_star = diag([d1_star_sym d2_star_sym]);
gamma = [gamma_1 ; gamma_2];
d_dot = [d1_dot ; d2_dot];

A_ = [A  -B*U ;
     zeros(n_u,n_x) zeros(n_u,n_u)];
B_ = [B B*D_star ;
      zeros(n_u,n_u) zeros(n_u,n_u)];
C_ = [C zeros(n_y,n_u)];

% Vertex system matrices for each region
A_region1_Q = zeros(n_x+n_u, n_x+n_u, n_u^2);
A_region2_Q = zeros(n_x+n_u, n_x+n_u, n_u^2);
A_region3_Q = zeros(n_x+n_u, n_x+n_u, n_u^2);
A_region4_Q = zeros(n_x+n_u, n_x+n_u, n_u^2);

A_region1_R = zeros(n_x+n_u, n_x+n_u, n_u^2);
A_region2_R = zeros(n_x+n_u, n_x+n_u, n_u^2);
A_region3_R = zeros(n_x+n_u, n_x+n_u, n_u^2);
A_region4_R = zeros(n_x+n_u, n_x+n_u, n_u^2);

% Region 1
k = 1;
for temp_u1 = [epsilon u1dot_max]
    for temp_u2 = [epsilon u2dot_max]
        A_numeric = double(subs(A_,[u1dot, u2dot], [temp_u1, temp_u2]));
        A_region1_Q(:,:,k) = A_numeric;
        k = k+1;
    end
end
k = 1;
for temp_u1 = [-epsilon u1dot_max]
    for temp_u2 = [-epsilon u2dot_max]
        A_numeric = double(subs(A_,[u1dot, u2dot], [temp_u1, temp_u2]));
        A_region1_R(:,:,k) = A_numeric;
        k = k+1;
    end
end

% Region 2
k = 1;
for temp_u1 = [u1dot_min -epsilon]
    for temp_u2 = [epsilon u2dot_max]
        A_numeric = double(subs(A_,[u1dot, u2dot], [temp_u1, temp_u2]));
        A_region2_Q(:,:,k) = A_numeric;
        k = k+1;
    end
end
k = 1;
for temp_u1 = [u1dot_min epsilon]
    for temp_u2 = [-epsilon u2dot_max]
        A_numeric = double(subs(A_,[u1dot, u2dot], [temp_u1, temp_u2]));
        A_region2_R(:,:,k) = A_numeric;
        k = k+1;
    end
end

% Region 3
k = 1;
for temp_u1 = [u1dot_min -epsilon]
    for temp_u2 = [u2dot_min -epsilon]
        A_numeric = double(subs(A_,[u1dot, u2dot], [temp_u1, temp_u2]));
        A_region3_Q(:,:,k) = A_numeric;
        k = k+1;
    end
end
k = 1;
for temp_u1 = [u1dot_min epsilon]
    for temp_u2 = [u2dot_min epsilon]
        A_numeric = double(subs(A_,[u1dot, u2dot], [temp_u1, temp_u2]));
        A_region3_R(:,:,k) = A_numeric;
        k = k+1;
    end
end

% Region 4
k = 1;
for temp_u1 = [epsilon u1dot_max]
    for temp_u2 = [u2dot_min -epsilon]
        A_numeric = double(subs(A_,[u1dot, u2dot], [temp_u1, temp_u2]));
        A_region4_Q(:,:,k) = A_numeric;
        k = k+1;
    end
end
k = 1;
for temp_u1 = [-epsilon u1dot_max]
    for temp_u2 = [u2dot_min epsilon]
        A_numeric = double(subs(A_,[u1dot, u2dot], [temp_u1, temp_u2]));
        A_region4_R(:,:,k) = A_numeric;
        k = k+1;
    end
end

% Group into cell arrays
A_region_Q = {A_region1_Q, A_region2_Q, A_region3_Q, A_region4_Q};
A_region_R = {A_region1_R, A_region2_R, A_region3_R, A_region4_R};

% Clear YALMIP's internal database 
yalmip('clear')

% Define variables at each vertex of each R-region
for k = 1:2^n_u
    X_region1_R{k} = sdpvar(n_x+n_u,n_y,'full');
    X_region2_R{k} = sdpvar(n_x+n_u,n_y,'full');
    X_region3_R{k} = sdpvar(n_x+n_u,n_y,'full');
    X_region4_R{k} = sdpvar(n_x+n_u,n_y,'full');
end
X_region_R = {X_region1_R, X_region2_R, X_region3_R, X_region4_R};

% Define variables for symmetric Lyapunov matrices for each region
P1 = sdpvar(n_x+n_u);
P2 = sdpvar(n_x+n_u);
P3 = sdpvar(n_x+n_u);
P4 = sdpvar(n_x+n_u);
P = {P1, P2, P3, P4};

% Initialize variables at each vertex of each Q-region
for k = 1:2^n_u
    X_region1_Q{k} = 0;
    X_region2_Q{k} = 0;
    X_region3_Q{k} = 0;
    X_region4_Q{k} = 0;
end

% Bilinear interpolation of observer gain vertex matrices in Q-regions from vertex gains in R-regions
for i = 1:2^n_u % Over all regions
    j=1;
    k=1;
    if i == 1 % Region 1
        u1dot_min_R = -epsilon;
        u1dot_max_R = u1dot_max;
        u2dot_min_R = -epsilon;
        u2dot_max_R = u2dot_max;
        for u1dot_Qv = [epsilon u1dot_max] 
            for u2dot_Qv = [epsilon u2dot_max]
                for temp_udot1 = [(u1dot_max_R-u1dot_Qv)/(u1dot_max_R-u1dot_min_R) (u1dot_Qv-u1dot_min_R)/(u1dot_max_R-u1dot_min_R)]
                    for temp_udot2 = [(u2dot_max_R-u2dot_Qv)/(u2dot_max_R-u2dot_min_R) (u2dot_Qv-u2dot_min_R)/(u2dot_max_R-u2dot_min_R)]
                        X_region1_Q{j} = X_region1_Q{j}+temp_udot1*temp_udot2*X_region1_R{k};
                        k = k+1;
                    end
                end
                j = j+1;
                k=1;
            end
        end
    elseif i == 2 % Region 2
        u1dot_min_R = u1dot_min;
        u1dot_max_R = epsilon;
        u2dot_min_R = -epsilon;
        u2dot_max_R = u2dot_max;
        for u1dot_Qv = [u1dot_min -epsilon]
            for u2dot_Qv = [epsilon u2dot_max]
                for temp_udot1 = [(u1dot_max_R-u1dot_Qv)/(u1dot_max_R-u1dot_min_R) (u1dot_Qv-u1dot_min_R)/(u1dot_max_R-u1dot_min_R)]
                    for temp_udot2 = [(u2dot_max_R-u2dot_Qv)/(u2dot_max_R-u2dot_min_R) (u2dot_Qv-u2dot_min_R)/(u2dot_max_R-u2dot_min_R)]
                        X_region2_Q{j} = X_region2_Q{j}+temp_udot1*temp_udot2*X_region2_R{k};
                        k = k+1;
                    end
                end
                j = j+1;
                k=1;
            end
        end
    elseif i == 3 % Region 3
        u1dot_min_R = u1dot_min;
        u1dot_max_R = epsilon;
        u2dot_min_R = u2dot_min;
        u2dot_max_R= epsilon;
        for u1dot_Qv = [u1dot_min -epsilon]
            for u2dot_Qv = [u2dot_min -epsilon]
                for temp_udot1 = [(u1dot_max_R-u1dot_Qv)/(u1dot_max_R-u1dot_min_R) (u1dot_Qv-u1dot_min_R)/(u1dot_max_R-u1dot_min_R)]
                    for temp_udot2 = [(u2dot_max_R-u2dot_Qv)/(u2dot_max_R-u2dot_min_R) (u2dot_Qv-u2dot_min_R)/(u2dot_max_R-u2dot_min_R)]
                        X_region3_Q{j} = X_region3_Q{j}+temp_udot1*temp_udot2*X_region3_R{k};
                        k = k+1;
                    end
                end
                j = j+1;
                k=1;
            end
        end
    else % Region 4
        u1dot_min_R = -epsilon;
        u1dot_max_R = u1dot_max;
        u2dot_min_R = u2dot_min;
        u2dot_max_R = epsilon;
        for u1dot_Qv = [epsilon u1dot_max]
            for u2dot_Qv = [u2dot_min -epsilon]
                for temp_udot1 = [(u1dot_max_R-u1dot_Qv)/(u1dot_max_R-u1dot_min_R) (u1dot_Qv-u1dot_min_R)/(u1dot_max_R-u1dot_min_R)]
                    for temp_udot2 = [(u2dot_max_R-u2dot_Qv)/(u2dot_max_R-u2dot_min_R) (u2dot_Qv-u2dot_min_R)/(u2dot_max_R-u2dot_min_R)]
                        X_region4_Q{j} = X_region4_Q{j}+temp_udot1*temp_udot2*X_region4_R{k};
                        k = k+1;
                    end
                end
                j = j+1;
                k=1;
            end
        end
    end
end
X_region_Q = {X_region1_Q, X_region2_Q, X_region3_Q, X_region4_Q};

% Define scalar gain variable delta
delta = sdpvar(1);

% Define identity and small positive definite matrices for LMIs
I_P = eye(n_x+n_u);
P_eps = 1e-3*eye(n_x+n_u);

% Add LMIs
List_LMIs = [delta >= 1e-3, P1 >= P_eps, P2 >= P_eps, P3 >= P_eps, P4 >= P_eps, ...
    alpha*I_P <= P1, P1 <= beta*I_P, ...
    alpha*I_P <= P2, P2 <= beta*I_P, ...
    alpha*I_P <= P3, P3 <= beta*I_P, ...
    alpha*I_P <= P4, P4 <= beta*I_P, ...
    exp(mu)*P2 - P1 >= 0, exp(mu)*P3 - P1 >= 0, exp(mu)*P4 - P1 >= 0, ...
    exp(mu)*P1 - P2 >= 0, exp(mu)*P3 - P2 >= 0, exp(mu)*P4 - P2 >= 0, ...
    exp(mu)*P1 - P3 >= 0, exp(mu)*P2 - P3 >= 0, exp(mu)*P4 - P3 >= 0, ...
    exp(mu)*P1 - P4 >= 0, exp(mu)*P2 - P4 >= 0, exp(mu)*P3 - P4 >= 0];

% Add LMIs for each region and vertex
for i = 1:2^n_u  % Number of regions
    for k = 1:2^n_u % Number of vertices
        A_Q = A_region_Q{i}(:,:,k);  
        A_R = A_region_R{i}(:,:,k);
        Pi = P{i};
        X_Q = X_region_Q{i}{k};  
        X_R = X_region_R{i}{k};
        
        % Add conditions
        List_LMIs = [List_LMIs, [A_Q'*Pi - C_'*X_Q' + Pi*A_Q - X_Q*C_ + lambda_minus*Pi -Pi ; -Pi -delta*I_P] <= 0];
        List_LMIs = [List_LMIs, [A_R'*Pi - C_'*X_R' + Pi*A_R - X_R*C_ - lambda_plus*Pi -Pi ; -Pi -delta*I_P] <= 0];
    end
end

% Set solver options
sdpoptions = sdpsettings('solver', 'sedumi', 'verbose', 1);

% Solve the LMIs
optimize(List_LMIs, [], sdpoptions)

% Extract numerical values
P1 = double(P1);
P2 = double(P2);
P3 = double(P3);
P4 = double(P4);

delta = value(delta);

%  Compute observer gains L = P^{-1} X at each vertex of each R-region
L_region1_R = zeros(n_x+n_u, n_y, 2^n_u);
L_region2_R = zeros(n_x+n_u, n_y, 2^n_u);
L_region3_R = zeros(n_x+n_u, n_y, 2^n_u);
L_region4_R = zeros(n_x+n_u, n_y, 2^n_u);
for k = 1:2^n_u
    L_region1_R(:,:,k) = P1\X_region1_R{k};
    L_region2_R(:,:,k) = P2\X_region2_R{k};
    L_region3_R(:,:,k) = P3\X_region3_R{k};
    L_region4_R(:,:,k) = P4\X_region4_R{k};
end
L_R = {L_region1_R, L_region2_R, L_region3_R, L_region4_R};

% Verify solution
P = {P1, P2, P3, P4};
condition = 0;

% Check conditions: exp(mu)*P{j} - P{i} >= 0
for i = 1:2^n_u
    for j = 1:2^n_u
        if  all(eig(exp(mu)*P{j} - P{i}) >= 0)
            % Condition satisfied
        else
            condition = condition + 1;
        end
    end
end

List_Eigenvalues_positive = [];
List_Eigenvalues_negative = [];
index_p = 1; index_n = 1;
for i = 1:2^n_u % Number of regions
    Pi = P{i};
    
    % Check conditions: alpha*I <= P_i <= beta*I
    List_Eigenvalues_positive(index_p) = min(real(eig(Pi - alpha*I_P))); index_p = index_p+1;
    List_Eigenvalues_positive(index_p) = min(real(eig(beta*I_P - Pi))); index_p = index_p+1;
    for k = 1:2^n_u % Number of vertices
        A_Q = A_region_Q{i}(:,:,k);  
        A_R = A_region_R{i}(:,:,k);
        X_Q = double(X_region_Q{i}{k});  
        X_R = double(X_region_R{i}{k});
        
        % Check LMIs for Q- and R-regions
        List_Eigenvalues_negative(index_n) = max(real(eig([A_Q'*Pi - C_'*X_Q' + Pi*A_Q - X_Q*C_ + lambda_minus*Pi -Pi ; -Pi -delta*I_P]))); index_n = index_n+1;
        List_Eigenvalues_negative(index_n) = max(real(eig([A_R'*Pi - C_'*X_R' + Pi*A_R - X_R*C_ - lambda_plus*Pi -Pi ; -Pi -delta*I_P]))); index_n = index_n+1;
    end
end

% Print result
if all(List_Eigenvalues_negative < 0) && all(List_Eigenvalues_positive >= 0) && condition == 0
    fprintf('Conditions verified\n');
else
    fprintf('Conditions NOT verified\n');
end


% Switching law conditions
activation_time_ratio = (lambda_plus+lambda_star)/(lambda_minus-lambda_star);
tau_a_star = mu/(lambda_star-lambda);

% Plot switching law condition summary
figure;
axis off;
hold on;
text(0.1, 1, '\bf Switching Law Conditions', 'FontSize', 20, 'Interpreter', 'latex');

% Parameters
text(0.1, 0.9, '\bf Parameters:', 'FontSize', 18, 'Interpreter', 'latex');
text(0.1, 0.85, sprintf('$\\lambda^+=$ %.3f', lambda_plus), 'FontSize', 16, 'Interpreter', 'latex');
text(0.1, 0.8, sprintf('$\\lambda^-=$ %.3f', lambda_minus), 'FontSize', 16, 'Interpreter', 'latex');
text(0.1, 0.75, sprintf('$\\lambda=$ %.3f', lambda), 'FontSize', 16, 'Interpreter', 'latex');
text(0.1, 0.7, sprintf('$\\lambda^*=$ %.3f', lambda_star), 'FontSize', 16, 'Interpreter', 'latex');

% Condition 1
text(0.1, 0.6, '$\textbf{Condition 1:}$', 'FontSize', 18, 'Interpreter', 'latex');
text(0.1, 0.525, ...
    sprintf('$\\frac{T_{\\mathcal{Q}}(t^j,t^{j+1})}{T_{\\mathcal{N}}(t^j,t^{j+1})}\\geq %.3f$', activation_time_ratio), ...
    'Interpreter', 'latex', 'FontSize', 16);

% Condition 2
text(0.1, 0.4, '$\textbf{Condition 2:}$', 'FontSize', 18, 'Interpreter', 'latex');
text(0.1, 0.35, sprintf('$\\tau_a^*=$ %.3f', tau_a_star), 'FontSize', 16, 'Interpreter', 'latex');

hold off;

%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all

% Select scenario to simulate
scenario = 1;

% Initial conditions
x0 = ones(n_x,1);
xhat0 = zeros(n_x,1); 
dhat0 = [d1_mid ; d2_mid];  
X0 = [x0 ; xhat0 ; dhat0];

% ODE solver options
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 1e-2);

% Run simulations with different choices of d^*_i
[T,Z1] = ode45(@(t,x)observer(t,x,scenario,1),0:0.1:50,X0,options);
d1_star = d1_mid; d2_star = d2_mid;
[T,Z2] = ode45(@(t,x)observer(t,x,scenario,0),0:0.1:50,X0,options);
d1_star = d1_min; d2_star = d2_min;
[T,Z3] = ode45(@(t,x)observer(t,x,scenario,0),0:0.1:50,X0,options);

% Define delays and input signals based on the selected scenario
if(scenario == 1)
    D1 = (T>=0&T<=15)*0.3+(T>15&T<=30)*0.7+(T>30)*0.5;
    D2 = (T>=0&T<=15)*0.3+(T>15&T<=30)*0.7+(T>30)*0.5;
    U1 = (T-D1);
    U2 = (T-D2);
    U1dot = (T>=0)*1;
    U2dot = (T>=0)*1;
elseif(scenario == 2)
    D1 = (T>=0&T<=15)*0.3+(T>15&T<=30)*0.7+(T>30)*0.5;
    D2 = (T>=0&T<=15)*0.3+(T>15&T<=30)*0.7+(T>30)*0.5;
    D1dot = 0;
    D2dot = 0;
    M1 = 1; w1 = 0.1;
    M2 = 1; w2 = 0.2;
elseif(scenario == 3)
    D1 = (T>=0&T<=15)*0.3+(T>15&T<=30)*0.7+(T>30)*0.5;
    D2 = (T>=0&T<=15)*0.3+(T>15&T<=30)*0.7+(T>30)*0.5;
    D1dot = 0;
    D2dot = 0;
    M1 = 10; w1 = 0.8;
    M2 = 10; w2 = 1;
elseif(scenario == 4)
    D1 = 0.5+0.3*sin(0.4*T);
    D2 = 0.5+0.3*cos(0.4*T);
    D1dot = 0.12*cos(0.4*T);
    D2dot = -0.12*sin(0.4*T);
    M1 = 1; w1 = 0.1;
    M2 = 1; w2 = 0.2;
elseif(scenario == 5)
    D1 = 0.5+0.3*sin(0.4*T);
    D2 = 0.5+0.3*cos(0.4*T);
    D1dot = 0.12*cos(0.4*T);
    D2dot = -0.12*sin(0.4*T);
    M1 = 10; w1 = 0.8;
    M2 = 10; w2 = 1;
elseif(scenario == 6)
    D1 = 0.5+0.3*sin(0.4*T);
    D2 = 0.5+0.3*cos(0.4*T);
    D1dot = 0.12*cos(0.4*T);
    D2dot = -0.12*sin(0.4*T);
    M1 = 0.1; w1 = pi/2;
    M2 = 0.1; w2 = pi/2;
elseif(scenario == 7)
    D1 = 0.5+0.3*sin(5*T);
    D2 = 0.5+0.3*cos(5*T);
    D1dot = 1.5*cos(5*T);
    D2dot = -1.5*sin(5*T);
    M1 = 1; w1 = 0.2;
    M2 = 1; w2 = 0.5;
end

% Saturated delay estimates
D1_hat_sat = max(min(Z1(:, 2*n_x+1), d1_max), d1_min);
D2_hat_sat = max(min(Z1(:, 2*n_x+2), d2_max), d2_min);

if any(scenario == [2,3,4,5,7])
    % True delayed inputs and their derivatives
    U1 = M1*sin(w1*(T-D1));
    U2 = M2*sin(w2*(T-D2));
    U1dot = M1*w1*cos(w1*(T-D1)).*(1-D1dot);
    U2dot = M2*w2*cos(w2*(T-D2)).*(1-D2dot);

    % Inputs delayed with saturated estimates of delay
    U1_est = M1*sin(w1*(T-D1_hat_sat));
    U2_est = M2*sin(w2*(T-D2_hat_sat));
    U1dot_est = M1*w1*cos(w1*(T-D1_hat_sat));
    U2dot_est = M2*w2*cos(w2*(T-D2_hat_sat));

    % Inputs delayed with d^* = d_mid
    U1_mid = M1*sin(w1*(T-d1_mid));
    U2_mid = M2*sin(w2*(T-d2_mid));
    U1dot_mid = M1*w1*cos(w1*(T-d1_mid));
    U2dot_mid = M2*w2*cos(w2*(T-d2_mid));

    % Inputs delayed with d^* = d_min
    U1_min = M1*sin(w1*(T-d1_min));
    U2_min = M2*sin(w2*(T-d2_min));
    U1dot_min = M1*w1*cos(w1*(T-d1_min));
    U2dot_min = M2*w2*cos(w2*(T-d2_min));
elseif any(scenario == [6])
    % True delayed inputs and their derivatives
    U1 = M1*sin(w1*(T-D1));
    U2 = M2*cos(w2*(T-D2));
    U1dot = M1*w1*cos(w1*(T-D1)).*(1-D1dot);
    U2dot = -M2*w2*sin(w2*(T-D2)).*(1-D2dot);

    % Inputs delayed with saturated estimates of delay
    U1_est = M1*sin(w1*(T-D1_hat_sat));
    U2_est = M2*cos(w2*(T-D2_hat_sat));
    U1dot_est = M1*w1*cos(w1*(T-D1_hat_sat));
    U2dot_est = -M2*w2*sin(w2*(T-D2_hat_sat));

    % Inputs delayed with d^* = d_mid
    U1_mid = M1*sin(w1*(T-d1_mid));
    U2_mid = M2*cos(w2*(T-d2_mid));
    U1dot_mid = M1*w1*cos(w1*(T-d1_mid));
    U2dot_mid = -M2*w2*sin(w2*(T-d2_mid));

    % Inputs delayed with d^* = d_min
    U1_min = M1*sin(w1*(T-d1_min));
    U2_min = M2*cos(w2*(T-d2_min));
    U1dot_min = M1*w1*cos(w1*(T-d1_min));
    U2dot_min = -M2*w2*sin(w2*(T-d2_min));
end

% Set default figure and font settings
set(0, 'DefaultLineLineWidth', 2, 'DefaultAxesFontSize', 28, 'DefaultTextFontSize', 28);

% Plot true states and their estimates
figure;
for i = 1:n_x
    subplot(n_x, 1, i);
    plot(T, Z1(:,i), T, Z1(:,n_x+i), T, Z2(:,n_x+i),'--', T, Z3(:,n_x+i),'-.');
    grid on
    if i == n_x
        xlabel('time (s)', 'Interpreter', 'latex');
    end
    ylabel(sprintf('$x_%d(t)$', i), 'Interpreter', 'latex');
    if i == 1
        legend('real', '$d_i^*(t) = \mathrm{sat}(\hat{d}_i(t^-))$', '$d_i^*(t) = (\underline{d}_i+\overline{d}_i)/2$',...
        '$d_i^*(t) = \underline{d}_i$', 'Interpreter', 'latex');
    end
end

% Plot true delays and their estimates
figure;
for i = 1:n_u
    subplot(n_u, 1, i);
    p = plot(T, eval(sprintf('D%d', i)), T, Z1(:,2*n_x+i), T, Z2(:,2*n_x+i),'--', T, Z3(:,2*n_x+i),'-.');
    grid on
    if i == n_u
        xlabel('time (s)', 'Interpreter', 'latex');
    end
    ylabel(sprintf('$d_%d(t)$', i), 'Interpreter', 'latex');
    if i == 1
        legend('real', '$d_i^*(t) = \mathrm{sat}(\hat{d}_i(t^-))$', '$d_i^*(t) = (\underline{d}_i+\overline{d}_i)/2$',...
        '$d_i^*(t) = \underline{d}_i$', 'Interpreter', 'latex');
    end
    if any(scenario == [1,2,3,4,5])
        ylim([0 1]);
    elseif scenario == 6
        ylim([-0.5 1.8]);
    elseif scenario == 7
        ylim([0 1.1]);
    end
end

if scenario ~= 1
    % Plot derivative of input signals
    default_colors = get(groot, 'defaultAxesColorOrder');
    custom_colors = [default_colors(2,:); default_colors(3,:); default_colors(4:end,:)];
    
    figure;
    set(gcf, 'DefaultAxesColorOrder', custom_colors);
    hold on;
    for i = 1:n_u
        subplot(n_u, 1, i);
        plot(T, eval(sprintf('U%ddot_est', i)), T, eval(sprintf('U%ddot_mid', i)), '--', T, eval(sprintf('U%ddot_min', i)), '-.');
        grid on
        if i == n_u
            xlabel('time (s)', 'Interpreter', 'latex');
        end
        ylabel(sprintf('$\\dot{u}_{%d}(t - d_{%d}^{*}(t))$', i, i), 'Interpreter', 'latex');
        if i == 1
            legend('$d_i^*(t) = \mathrm{sat}(\hat{d}_i(t^-))$', '$d_i^*(t) = (\underline{d}_i+\overline{d}_i)/2$',...
                '$d_i^*(t) = \underline{d}_i$', 'Interpreter', 'latex');
        end
    end
    hold off;
end

% Compute IAE and ISE of delays for t >= 1.5
start_idx = find(T >= 1.5, 1);
T_start = T(start_idx:end);

% Delays
IAE_d1_est = trapz(T_start, abs(D1(start_idx:end) - Z1(start_idx:end, 2*n_x+1)));
ISE_d1_est = trapz(T_start, (D1(start_idx:end) - Z1(start_idx:end, 2*n_x+1)).^2);
IAE_d2_est = trapz(T_start, abs(D2(start_idx:end) - Z1(start_idx:end, 2*n_x+2)));
ISE_d2_est = trapz(T_start, (D2(start_idx:end) - Z1(start_idx:end, 2*n_x+2)).^2);

IAE_d1_mid = trapz(T_start, abs(D1(start_idx:end) - Z2(start_idx:end, 2*n_x+1)));
ISE_d1_mid = trapz(T_start, (D1(start_idx:end) - Z2(start_idx:end, 2*n_x+1)).^2);
IAE_d2_mid = trapz(T_start, abs(D2(start_idx:end) - Z2(start_idx:end, 2*n_x+2)));
ISE_d2_mid = trapz(T_start, (D2(start_idx:end) - Z2(start_idx:end, 2*n_x+2)).^2);

IAE_d1_min = trapz(T_start, abs(D1(start_idx:end) - Z3(start_idx:end, 2*n_x+1)));
ISE_d1_min = trapz(T_start, (D1(start_idx:end) - Z3(start_idx:end, 2*n_x+1)).^2);
IAE_d2_min = trapz(T_start, abs(D2(start_idx:end) - Z3(start_idx:end, 2*n_x+2)));
ISE_d2_min = trapz(T_start, (D2(start_idx:end) - Z3(start_idx:end, 2*n_x+2)).^2);

% States
IAE_x1_est = trapz(T_start, abs(Z1(start_idx:end, 1) - Z1(start_idx:end, n_x+1)));
ISE_x1_est = trapz(T_start, (Z1(start_idx:end, 1) - Z1(start_idx:end, n_x+1)).^2);
IAE_x2_est = trapz(T_start, abs(Z1(start_idx:end, 2) - Z1(start_idx:end, n_x+2)));
ISE_x2_est = trapz(T_start, (Z1(start_idx:end, 2) - Z1(start_idx:end, n_x+2)).^2);
IAE_x3_est = trapz(T_start, abs(Z1(start_idx:end, 3) - Z1(start_idx:end, n_x+3)));
ISE_x3_est = trapz(T_start, (Z1(start_idx:end, 3) - Z1(start_idx:end, n_x+3)).^2);
IAE_x4_est = trapz(T_start, abs(Z1(start_idx:end, 4) - Z1(start_idx:end, n_x+4)));
ISE_x4_est = trapz(T_start, (Z1(start_idx:end, 4) - Z1(start_idx:end, n_x+4)).^2);
IAE_x5_est = trapz(T_start, abs(Z1(start_idx:end, 5) - Z1(start_idx:end, n_x+5)));
ISE_x5_est = trapz(T_start, (Z1(start_idx:end, 5) - Z1(start_idx:end, n_x+5)).^2);
IAE_x6_est = trapz(T_start, abs(Z1(start_idx:end, 6) - Z1(start_idx:end, n_x+6)));
ISE_x6_est = trapz(T_start, (Z1(start_idx:end, 6) - Z1(start_idx:end, n_x+6)).^2);

IAE_x1_mid = trapz(T_start, abs(Z2(start_idx:end, 1) - Z2(start_idx:end, n_x+1)));
ISE_x1_mid = trapz(T_start, (Z2(start_idx:end, 1) - Z2(start_idx:end, n_x+1)).^2);
IAE_x2_mid = trapz(T_start, abs(Z2(start_idx:end, 2) - Z2(start_idx:end, n_x+2)));
ISE_x2_mid = trapz(T_start, (Z2(start_idx:end, 2) - Z2(start_idx:end, n_x+2)).^2);
IAE_x3_mid = trapz(T_start, abs(Z2(start_idx:end, 3) - Z2(start_idx:end, n_x+3)));
ISE_x3_mid = trapz(T_start, (Z2(start_idx:end, 3) - Z2(start_idx:end, n_x+3)).^2);
IAE_x4_mid = trapz(T_start, abs(Z2(start_idx:end, 4) - Z2(start_idx:end, n_x+4)));
ISE_x4_mid = trapz(T_start, (Z2(start_idx:end, 4) - Z2(start_idx:end, n_x+4)).^2);
IAE_x5_mid = trapz(T_start, abs(Z2(start_idx:end, 5) - Z2(start_idx:end, n_x+5)));
ISE_x5_mid = trapz(T_start, (Z2(start_idx:end, 5) - Z2(start_idx:end, n_x+5)).^2);
IAE_x6_mid = trapz(T_start, abs(Z2(start_idx:end, 6) - Z2(start_idx:end, n_x+6)));
ISE_x6_mid = trapz(T_start, (Z2(start_idx:end, 6) - Z2(start_idx:end, n_x+6)).^2);

IAE_x1_min = trapz(T_start, abs(Z3(start_idx:end, 1) - Z3(start_idx:end, n_x+1)));
ISE_x1_min = trapz(T_start, (Z3(start_idx:end, 1) - Z3(start_idx:end, n_x+1)).^2);
IAE_x2_min = trapz(T_start, abs(Z3(start_idx:end, 2) - Z3(start_idx:end, n_x+2)));
ISE_x2_min = trapz(T_start, (Z3(start_idx:end, 2) - Z3(start_idx:end, n_x+2)).^2);
IAE_x3_min = trapz(T_start, abs(Z3(start_idx:end, 3) - Z3(start_idx:end, n_x+3)));
ISE_x3_min = trapz(T_start, (Z3(start_idx:end, 3) - Z3(start_idx:end, n_x+3)).^2);
IAE_x4_min = trapz(T_start, abs(Z3(start_idx:end, 4) - Z3(start_idx:end, n_x+4)));
ISE_x4_min = trapz(T_start, (Z3(start_idx:end, 4) - Z3(start_idx:end, n_x+4)).^2);
IAE_x5_min = trapz(T_start, abs(Z3(start_idx:end, 5) - Z3(start_idx:end, n_x+5)));
ISE_x5_min = trapz(T_start, (Z3(start_idx:end, 5) - Z3(start_idx:end, n_x+5)).^2);
IAE_x6_min = trapz(T_start, abs(Z3(start_idx:end, 6) - Z3(start_idx:end, n_x+6)));
ISE_x6_min = trapz(T_start, (Z3(start_idx:end, 6) - Z3(start_idx:end, n_x+6)).^2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = observer(t,x,scenario,time_varying)
global n_x n_u n_y A B C d1_min d1_max d1_star d2_min d2_max d2_star u1dot_min u1dot_max ... 
    u2dot_min u2dot_max epsilon L_R
persistent considered_region

% Define delays and input signals based on the selected scenario
if(scenario == 1)
    if(t<=15)
        d1 = 0.3;
        d2 = 0.3;
    elseif(t<=30)
        d1 = 0.7;
        d2 = 0.7;
    else
        d1 = 0.5;
        d2 = 0.5;
    end
    u1 = (t-d1);
    u2 = (t-d2);
    u1_ = (t-d1_star);
    u2_ = (t-d2_star);
    u1dot_ = 1;
    u2dot_ = 1;
elseif(scenario==2)
    if(t<=15)
        d1 = 0.3;
        d2 = 0.3;
    elseif(t<=30)
        d1 = 0.7;
        d2 = 0.7;
    else
        d1 = 0.5;
        d2 = 0.5;
    end
    M1 = 1; w1 = 0.1;
    M2 = 1; w2 = 0.2;
elseif(scenario==3)
    if(t<=15)
        d1 = 0.3;
        d2 = 0.3;
    elseif(t<=30)
        d1 = 0.7;
        d2 = 0.7;
    else
        d1 = 0.5;
        d2 = 0.5;
    end
    M1 = 10; w1 = 0.8;
    M2 = 10; w2 = 1;
elseif(scenario==4)
    d1 = 0.5+0.3*sin(0.4*t);
    d2 = 0.5+0.3*cos(0.4*t);
    M1 = 1; w1 = 0.1;
    M2 = 1; w2 = 0.2;
elseif(scenario==5)
    d1 = 0.5+0.3*sin(0.4*t);
    d2 = 0.5+0.3*cos(0.4*t);
    M1 = 10; w1 = 0.8;
    M2 = 10; w2 = 1;
elseif(scenario==6)
    d1 = 0.5+0.3*sin(0.4*t);
    d2 = 0.5+0.3*cos(0.4*t);
    M1 = 0.1; w1 = pi/2;
    M2 = 0.1; w2 = pi/2;
elseif(scenario==7)
    d1 = 0.5+0.3*sin(5*t);
    d2 = 0.5+0.3*cos(5*t);
    M1 = 1; w1 = 0.2;
    M2 = 1; w2 = 0.5;
end

if time_varying==1 % Use saturated delay estimates if time_varying = 1
    d1_hat = x(2*n_x+1); d2_hat = x(2*n_x+n_u);
    if d1_hat > d1_max
        d1_hat = d1_max;
    elseif d1_hat < d1_min
        d1_hat = d1_min;
    end
    if d2_hat > d2_max
        d2_hat = d2_max;
    elseif d2_hat < d2_min
        d2_hat = d2_min;
    end
    d1_star = d1_hat; d2_star = d2_hat;
end

if any(scenario == [2,3,4,5,7]) 
    u1 = M1*sin(w1*(t-d1));
    u2 = M2*sin(w2*(t-d2));
    u1_ = M1*sin(w1*(t-d1_star));
    u2_ = M2*sin(w2*(t-d2_star));
    u1dot_ = M1*w1*cos(w1*(t-d1_star));
    u2dot_ = M2*w2*cos(w2*(t-d2_star));
elseif any(scenario == [6])
    u1 = M1*sin(w1*(t-d1));
    u2 = M2*cos(w2*(t-d2));
    u1_ = M1*sin(w1*(t-d1_star));
    u2_ = M2*cos(w2*(t-d2_star));
    u1dot_ = M1*w1*cos(w1*(t-d1_star));
    u2dot_ = -M2*w2*sin(w2*(t-d2_star));
end

% Determine active region
previous_region = considered_region;
if(t == 0)
        % Initial region
        if(u1dot_ >= 0 && u2dot_  >=0)
            considered_region = 1;
        elseif(u1dot_  < 0 && u2dot_ >= 0)
            considered_region = 2;
        elseif(u1dot_  < 0 && u2dot_ < 0)
            considered_region = 3;
        else
            considered_region = 4;
        end
    else
        if(considered_region == 1)
            if(u1dot_ < -epsilon && u2dot_ > -epsilon)
                considered_region = 2;
            elseif(u1dot_ < -epsilon && u2dot_ < -epsilon)
                considered_region = 3;
            elseif(u1dot_ > -epsilon && u2dot_ < -epsilon)
                considered_region = 4;
            else
                considered_region = 1;
            end
        elseif(considered_region == 2)
            if(u1dot_ > epsilon && u2dot_ > -epsilon)
                considered_region = 1;
            elseif(u1dot_ < epsilon && u2dot_ < -epsilon)
                considered_region = 3;
            elseif(u1dot_ > epsilon && u2dot_ < -epsilon)
                considered_region = 4;
            else
                considered_region = 2;
            end
        elseif(considered_region == 3)
            if(u1dot_ > epsilon && u2dot_ > epsilon)
                considered_region = 1;
            elseif(u1dot_ < epsilon && u2dot_ > epsilon)
                considered_region = 2;
            elseif(u1dot_ > epsilon && u2dot_ < epsilon)
                considered_region = 4;
            else
                considered_region = 3;
            end
        elseif(considered_region == 4)
            if(u1dot_ > -epsilon && u2dot_ > epsilon)
                considered_region = 1;
            elseif(u1dot_ < -epsilon && u2dot_ > epsilon)
                considered_region = 2;
            elseif(u1dot_ < -epsilon && u2dot_ < epsilon)
                considered_region = 3;
            else
                considered_region = 4;
            end
        end
end

% Compute observer gain by bilinear interpolation based on active region
k = 1; L0 = zeros(n_x+n_u,n_y);
L_region1 = L_R{1};
L_region2 = L_R{2};
L_region3 = L_R{3};
L_region4 = L_R{4};
if(considered_region == 1)
    for temp_udot1 = [(u1dot_max-u1dot_)/(u1dot_max-(-epsilon)) (u1dot_-(-epsilon))/(u1dot_max-(-epsilon))]
        for temp_udot2 = [(u2dot_max-u2dot_)/(u2dot_max-(-epsilon)) (u2dot_-(-epsilon))/(u2dot_max-(-epsilon))]
                L0 = L0+temp_udot1*temp_udot2*L_region1(:,:,k);
                k = k+1;
        end
    end
elseif(considered_region == 2)
    for temp_udot1 = [(epsilon-u1dot_)/(epsilon-u1dot_min) (u1dot_-u1dot_min)/(epsilon-u1dot_min)]
        for temp_udot2 = [(u2dot_max-u2dot_)/(u2dot_max-(-epsilon)) (u2dot_-(-epsilon))/(u2dot_max-(-epsilon))]
                L0 = L0+temp_udot1*temp_udot2*L_region2(:,:,k);
                k = k+1;
        end
    end
elseif(considered_region == 3)
    for temp_udot1 = [(epsilon-u1dot_)/(epsilon-u1dot_min) (u1dot_-u1dot_min)/(epsilon-u1dot_min)]
        for temp_udot2 = [(epsilon-u2dot_)/(epsilon-u2dot_min) (u2dot_-u2dot_min)/(epsilon-u2dot_min)]
                L0 = L0+temp_udot1*temp_udot2*L_region3(:,:,k);
                k = k+1;
        end
    end
elseif(considered_region == 4)
    for temp_udot1 = [(u1dot_max-u1dot_)/(u1dot_max-(-epsilon)) (u1dot_-(-epsilon))/(u1dot_max-(-epsilon))]
        for temp_udot2 = [(epsilon-u2dot_)/(epsilon-u2dot_min) (u2dot_-u2dot_min)/(epsilon-u2dot_min)]
                L0 = L0+temp_udot1*temp_udot2*L_region4(:,:,k);
                k = k+1;
        end
    end
end

% Construct augmented system matrices for the observer
U = diag([u1dot_ u2dot_]); 
D_star = diag([d1_star d2_star]);

A_ = [A  -B*U ;
     zeros(n_u,n_x) zeros(n_u,n_u)];
B_ = [B B*D_star ;
      zeros(n_u,n_u) zeros(n_u,n_u)];
C_ = [C zeros(n_y,n_u)];

% State vector structure:
%   x(1:n_x)             = states of the system
%   x(n_x+1:2*n_x)       = estimates of the states
%   x(2*n_x+1:2*n_x+n_u) = estimates of the delays

dx = zeros(2*n_x+n_u,1);

% Dynamics of the true system
dx(1:n_x) = A*x(1:n_x) + B*[u1 ; u2];

% Dynamics of the observer
dx(n_x+1:2*n_x+n_u) = A_*x(n_x+1:2*n_x+n_u)+B_*[u1_ ; u2_ ; u1dot_ ; u2dot_]...
    -L0*(C_*x(n_x+1:2*n_x+n_u)-C*x(1:n_x));

end
