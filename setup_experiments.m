% Clear workspace
clear all
close all
clc

global n_x n_u n_y A B C d1_min d1_max u1dot_min u1dot_max ... 
    epsilon L_R

% Inputs signal
M = 5; % Magnitude
w = 1; % Frequency

% System parameters
Rm = 8.4; % Resistance
kt = 0.042; % Current-torque (N-m/A)
ke = 0.042; % Back-emf constant (V-s/rad)
kd = 1*10^(-5); % Propeller Drag (N-m-s/rad)
Jr = 4e-6; % Rotor inertia (kg-m^2)
mh = 0.003; % Hub mass (kg)
rh = 9/1000/2; % Hub radius (m)
Jh = 3.04e-9; % Hub inertia (kg-m^2)
Jp = 7.2*10^(-6); % Propeller moment of inertia (kg-m^2)
Jeq = Jr + Jh + Jp; % Equivalent moment of inertia (kg-m^2)

% Define system dimensions
n_x = 1;
n_u = 1;
n_y = 1;

% Define system matrices
A = -kt*ke/(Jeq*Rm) - kd/Jeq;
%B = kt/(Jeq*Rm);
C = 1;

% Experimental input matrix
B = 23*sqrt(4+A^2);

% Define delay and input derivative bounds
d1_min = 0.2; d1_max = 0.8; d1_mid = (d1_min + d1_max)/2;
u1dot_min = -5; u1dot_max = 5;

% Defining constants
epsilon = 0.05;
lambda_plus = 0.2;
lambda_minus = 0.15;
mu = 0.2;

alpha = 1;
beta = 1000;

rho = 250;

% Switching law constants
lambda = 0.25*lambda_minus;
lambda_star = 0.75*lambda_minus;

% Augmented system (symbolic)
syms u1dot u2dot real 

U = diag(u1dot); 

A_ = [A  -B*U ;
     zeros(n_u,n_x) zeros(n_u,n_u)];
C_ = [C zeros(n_y,n_u)];

% Vertex system matrices for each region
A_region1_Q = zeros(n_x+n_u, n_x+n_u, 2^n_u);
A_region2_Q = zeros(n_x+n_u, n_x+n_u, 2^n_u);
A_region1_R = zeros(n_x+n_u, n_x+n_u, 2^n_u);
A_region2_R = zeros(n_x+n_u, n_x+n_u, 2^n_u);

% Region 1
k = 1;
for temp_u1 = [epsilon u1dot_max]
    A_numeric = double(subs(A_,u1dot, temp_u1));
    A_region1_Q(:,:,k) = A_numeric;
    k = k+1;
end
k = 1;
for temp_u1 = [-epsilon u1dot_max]
    A_numeric = double(subs(A_,u1dot, temp_u1));
    A_region1_R(:,:,k) = A_numeric;
    k = k+1;
end

% Region 2
k = 1;
for temp_u1 = [u1dot_min -epsilon]
    A_numeric = double(subs(A_,u1dot, temp_u1));
    A_region2_Q(:,:,k) = A_numeric;
    k = k+1;
end
k = 1;
for temp_u1 = [u1dot_min epsilon]
    A_numeric = double(subs(A_,u1dot, temp_u1));
    A_region2_R(:,:,k) = A_numeric;
    k = k+1;
end

% Group into cell arrays
A_region_Q = {A_region1_Q, A_region2_Q};
A_region_R = {A_region1_R, A_region2_R};

% Clear YALMIP's internal database 
yalmip('clear')

% Define variables at each vertex of each R-region
for k = 1:2^n_u
    X_region1_R{k} = sdpvar(n_x+n_u,n_y,'full');
    X_region2_R{k} = sdpvar(n_x+n_u,n_y,'full');
end
X_region_R = {X_region1_R, X_region2_R};

% Define variables for symmetric Lyapunov matrices for each region
P1 = sdpvar(n_x+n_u);
P2 = sdpvar(n_x+n_u);
P = {P1, P2};

% Initialize variables at each vertex of each Q-region
for k = 1:2^n_u
    X_region1_Q{k} = 0;
    X_region2_Q{k} = 0;
end

% Bilinear interpolation of observer gain vertex matrices in Q-regions from vertex gains in R-regions
for i = 1:2^n_u % Over all regions
    j=1;
    k=1;
    if i == 1 % Region 1
        u1dot_min_R = -epsilon;
        u1dot_max_R = u1dot_max;
        for u1dot_Qv = [epsilon u1dot_max] 
            for temp_udot1 = [(u1dot_max_R-u1dot_Qv)/(u1dot_max_R-u1dot_min_R) (u1dot_Qv-u1dot_min_R)/(u1dot_max_R-u1dot_min_R)]
                X_region1_Q{j} = X_region1_Q{j}+temp_udot1*X_region1_R{k};
                k = k+1;
            end
            j = j+1;
            k=1;
        end
    elseif i == 2 % Region 2
        u1dot_min_R = u1dot_min;
        u1dot_max_R = epsilon;
        for u1dot_Qv = [u1dot_min -epsilon]
            for temp_udot1 = [(u1dot_max_R-u1dot_Qv)/(u1dot_max_R-u1dot_min_R) (u1dot_Qv-u1dot_min_R)/(u1dot_max_R-u1dot_min_R)]
                X_region2_Q{j} = X_region2_Q{j}+temp_udot1*X_region2_R{k};
                k = k+1;
            end
            j = j+1;
            k=1;
        end
    end
end
X_region_Q = {X_region1_Q, X_region2_Q};

% Define scalar gain variable delta
delta = sdpvar(1);

% Define identity and small positive definite matrices for LMIs
I_P = eye(n_x+n_u);
P_eps = 1e-3*eye(n_x+n_u);

% Add LMIs
List_LMIs = [delta >= 1e-3, P1 >= P_eps, P2 >= P_eps, ...
    alpha*I_P <= P1, P1 <= beta*I_P, ...
    alpha*I_P <= P2, P2 <= beta*I_P, ...
    exp(mu)*P2 - P1 >= 0, ...
    exp(mu)*P1 - P2 >= 0];

% Add LMIs for each region and vertex
for i = 1:2^n_u  % Number of regions
    for k = 1:2^n_u % Number of vertices
        A_Q = A_region_Q{i}(:,:,k);  
        A_R = A_region_R{i}(:,:,k);
        Pi = P{i};
        X_Q = X_region_Q{i}{k};  
        X_R = X_region_R{i}{k};
        
        % Add conditions
        List_LMIs = [List_LMIs, [rho^2*eye(n_x+n_u) X_R; X_R' eye(n_y)] >= 0];
        List_LMIs = [List_LMIs, [A_Q'*Pi - C_'*X_Q' + Pi*A_Q - X_Q*C_ + lambda_minus*Pi -Pi ; -Pi -delta*I_P] <= 0];
        List_LMIs = [List_LMIs, [A_R'*Pi - C_'*X_R' + Pi*A_R - X_R*C_ - lambda_plus*Pi -Pi ; -Pi -delta*I_P] <= 0];
    end
end

% Set solver options
sdpoptions = sdpsettings('showprogress',1,'solver','sdpt3');

% Solve the LMIs
optimize(List_LMIs, [], sdpoptions)

% Extract numerical values
P1 = double(P1);
P2 = double(P2);

delta = value(delta);

%  Compute observer gains L = P^{-1} X at each vertex of each R-region
L_region1_R = zeros(n_x+n_u, n_y, 2^n_u);
L_region2_R = zeros(n_x+n_u, n_y, 2^n_u);
for k = 1:2^n_u
    L_region1_R(:,:,k) = P1\X_region1_R{k};
    L_region2_R(:,:,k) = P2\X_region2_R{k};
end
L_R = {L_region1_R, L_region2_R};

% Verify solution
P = {P1, P2};
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
        List_Eigenvalues_positive(index_p) =  min(real(eig([rho^2*eye(n_x+n_u) X_R; X_R' eye(n_y)]))); index_p = index_p+1;
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