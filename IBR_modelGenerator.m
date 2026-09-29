function [sys_j, sys_I, sys_c] = IBR_modelGenerator(G)
% G is a cell including the subsystem models

%% Load the subsystems
for ic = 1:length(G)
    sys_j(ic).s = G{ic};
end

%% Preparation of subsystems: splitting of stable and unstable parts
for ic = 1:length(G)
    [sys_j(ic).Gs, sys_j(ic).Gns, info] = stabsep(sys_j(ic).s, 'Offset', 1e-5);
    sys_j(ic).TL = info.TL;
    sys_j(ic).TR = info.TR;
end

%% Interconnection dynamics and coupling
% Parallel connection of subsystem models
sys_b.s = blkdiag(sys_j(1:end-1).s);

AA = sys_j(end).s.A;
BB = sys_j(end).s.B;
CC = sys_j(end).s.C;

n_studyArea = length(AA);

% Coupling matrices: including the study area
dim_ucyc = 1; % Number of external inputs and outputs
K = zeros(length(sys_b.s) + dim_ucyc);

K_A = AA;
K_B = zeros(n_studyArea, size(K, 2));
K_C = zeros(size(K, 1), n_studyArea);

K_B(:, 15) = BB(:, 1); K_B(:, 16) = BB(:, 2); K_B(:, 17) = BB(:, 3);
K_B(:, 18) = BB(:, 4); K_B(:, 19) = BB(:, 7); K_B(:, 1) = BB(:, 5);
K_B(:, 2) = BB(:, 6);

K_C(1, :) = CC(6, :); K_C(2, :) = CC(7, :); K_C(15, :) = CC(2, :);
K_C(16, :) = CC(3, :); K_C(17, :) = CC(4, :); K_C(18, :) = CC(5, :);
K_C(19, :) = CC(1, :);
 
K(3,9) = 1; K(4,10) = 1; K(5,11) = 1; K(6,12) = 1;                    
K(7,13) = 1; K(8,14) = 1; K(9,3) = 1; K(10,4) = 1; K(11,5) = 1; K(12,6) = 1;
K(13,7) = 1; K(14,8) = 1;

K_11 = K(1:length(sys_b.s), 1:length(sys_b.s));
K_12 = K(1:length(sys_b.s), length(sys_b.s)+1:end);
K_21 = K(length(sys_b.s)+1:end, 1:length(sys_b.s));
K_22 = K(length(sys_b.s)+1:end, length(sys_b.s)+1:end);

K_D = [K_11, K_12; K_21, K_22];

% Dynamical interface system 
sys_I.s = ss(K_A, K_B, K_C, K_D);

% Interconnected system: lft (upper lft by default) of interface and subsystems
sys_c.s = lft(sys_b.s, sys_I.s);

sys_j = sys_j(1:end-1);