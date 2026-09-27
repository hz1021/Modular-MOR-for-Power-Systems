function [Dlv, Wout, Vout, info, Wout_inv, Vout_inv] = DK_iteration_adapt(N0_f, m_jc, p_jc, e_jc, a_jc, maxits)
% Bottum-Up AND Top-Down algorithm combined by repetitively solving of V&W
% and D, where some error bounds can be given and the remaining bounds are
% calculated. Boulds are determined by diagonal (but not scaled identity)
% matrices V and W, such that Input-Output pairs can be weighted
% differently.
%
% N0_f [matrix] is the frequency response at this particular frequency point
%
% m_jc [row matrix] indicates the dimensions of the blocks in V
%
% p_jc [row matrix] indicates the dimensions of the blocks in W
%
% e_jc [cell array] indicates the error bounds:
%   - if e_jc{k} is a negative scalar, the function will calculate the
%   corresponding bound
%   - if e_jc{k} is not a negative scalar, it is treated as a bound
%   - if e_jc{k} and e_jc{l} are identical, the corresponding errors are
%   assumed to be identical and are treated as such in the composition of 
%   V, W and D
%
% ::: OPTIONAL::: 
%
% a_jc [row matrix] indicates optional additional weighting of the bounds 
%   in the cost function of the optimization, allowing manipulation of the 
%   optimizer.
%
% maxits [positive integer] default = 50, maximum number of iterations 
%   over V&W and D


%% Processing
% Every time you call sdpvar etc, an internal database grows larger
yalmip('clear')

% Find which bounds are given, and which are unknown
I_unkwn = cellfun(@(x) isscalar(x) && x < 0, e_jc);
I_known = ~I_unkwn;
I_scalar = cellfun(@(x) isscalar(x) && x < 0 && imag(x) ~= 0, e_jc);

% Find repetitive blocks
e_ujc = e_jc(1);
ia = 1;
ic = 1;
for ii = 2:length(e_jc)
    for k = 1:length(e_ujc)
        if isequal(e_jc{ii}, e_ujc{k})
            ic = [ic, k];
            break
        end
    end
    if ~isequal(e_jc{ii}, e_ujc{k})
        ia = [ia, ii];
        ic = [ic, ii];
        e_ujc = [e_ujc, e_jc(ii)];
    end    
end

% Indices of nominal system N0_f, split over each uncertainty
CSu = cumsum([0, m_jc]);
CSy = cumsum([0, p_jc]);

% Maximum parameter value in solver
pmax = 1e4;

% Stopping tolerance for convergence of the error bounds
tol = 1e-4;

% Indices corresponding to each subsystem/uncertainty
for k = 1:length(p_jc)
    indu{k} = CSu(k)+1:CSu(k+1);
    indy{k} = CSy(k)+1:CSy(k+1);
end

% Optional parameters
if nargin > 4 && ~isempty(a_jc)
    aV = arrayfun(@(x,y) x*eye(y), a_jc, m_jc, 'UniformOutput', false);
    aW = arrayfun(@(x,y) x*eye(y), a_jc, p_jc, 'UniformOutput', false);
    aV = blkdiag(aV{:});
    aW = blkdiag(aW{:});
else
    aV = eye(sum(m_jc));
    aW = eye(sum(p_jc));
end

% Maximum number of iterations
if ~exist('maxits', 'var')
    maxits = 50;
end

% sdp variable for D-optimisation
gamma = sdpvar(1, 1);

% Set options for YALMIP mosek solver
options = sdpsettings('verbose', 0, 'solver', 'mosek');


%% Construction of D
Dls(sum(m_jc), sum(m_jc)) = sdpvar(1, 1); % Initialise sdpvar matrix ???!!!
Drs(sum(p_jc), sum(p_jc)) = sdpvar(1, 1); % Initialise sdpvar matrix

% Find the scaling of D to fix as 1: Fix the one corresponding to the last
% given, non-repeated bound
i_D1 = find(ia == find(sum(ic == ic') == 1), 1, 'last');

% Creation of sdp variables: scalar if single uncertainty, 2x2 block if
% repetitive
for ii = 1:length(e_ujc)
    % Find the uncertainties below belong to which unique bound
    ind = find(ic == ia(ii));

    if ii == i_D1
        % Set this scaling to identity
        R{ii} = 1;
    else
        % Create the block of sdpvars
        R{ii} = sdpvar(length(ind), length(ind), 'hermitian', 'complex'); % ???!!!
    end

    % Current indices of D
    iu = [indu{ind}];
    iy = [indy{ind}];

    % Add sdpvars to D using kronecker product
    Dls(iu, iu) = kron(R{ii}, eye(m_jc(ind(1))));
    Drs(iy, iy) = kron(R{ii}, eye(p_jc(ind(1))));

    % Set the initial value of D
    Dlv(iu, iu) = kron(-0.1 * ones(size(R{ii})) + eye(size(R{ii})), eye(m_jc(ind(1))));
    Drv(iy, iy) = kron(-0.1 * ones(size(R{ii})) + eye(size(R{ii})), eye(p_jc(ind(1))));

end


%% Construction of W and V
% Set E_j = W_j \Delta_j V_j, assuming square systems by W_j = V_j =
% sqrt(e_j) * I
% Create sdpvars for the subsystem errors
Wi2s(sum(p_jc), sum(p_jc)) = sdpvar(1, 1); % Initialise sdpvar matrix
Vi2s(sum(m_jc), sum(m_jc)) = sdpvar(1, 1); % Initialise sdpvar matrix
for ii = 1:length(e_ujc)
    % Find the uncertainties below belong to which unique bound
    ind = find(ic == ia(ii));

    % Current indices of V and W
    iu = [indu{ind}];
    iy = [indy{ind}];

    % Add sdpvars to W = W_j^2 using kronecker product
    if I_known(ia(ii))
        if ind == length(e_jc)
            Wc = inv(diag(sqrt(max(e_ujc{ii}, [], 1)))); % ???!!!
            Vc = inv(diag(sqrt(max(e_ujc{ii}, [], 2))));
            VeWn = norm(Vc * e_ujc{ii} * Wc);
            Wji2 = inv(Wc^2 / VeWn);
            Vji2 = inv(Vc^2 / VeWn);
        else
            Wj = diag(sqrt(max(e_ujc{ii}, [], 2)));
            Vj = diag(sqrt(max(e_ujc{ii}, [], 1)));
            VeWn = norm(Wj \ e_ujc{ii} / Vj);
            Wji2 = inv(Wj^2 * VeWn);
            Vji2 = inv(Vj^2 * VeWn);
        end

        % Add the know bound
        % Wi2s(iy, iy) = kron(eye(length(ind) * p_jc(end)), Wji2);
        % Vi2s(iu, iu) = kron(eye(length(ind) * m_jc(end)), Vji2);
        Wi2s(iy, iy) = kron(eye(length(ind)), Wji2);
        Vi2s(iu, iu) = kron(eye(length(ind)), Vji2);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Add the know bound
        % interW = eye(p_jc(end)); Windex = 3; interW(Windex, Windex) = 100;
        % interV = eye(m_jc(end)); Vindex = 2; interV(Vindex, Vindex) = 100;
        % Wi2s(iy, iy) = interW^2;
        % Vi2s(iu, iu) = interV^2;
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif I_scalar(ia(ii))
        % Create the unknown sdpvar
        ws{ii, 1} = sdpvar(1, 1);
        vs{ii, 1} = sdpvar(1, 1);
        % Add the unknown sdpvar
        Wi2s(iy, iy) = kron(eye(length(ind)), ws{ii} * eye(p_jc(ind(1))));
        Vi2s(iu, iu) = kron(eye(length(ind)), vs{ii} * eye(m_jc(ind(1))));
    else
        % Create the unknown sdpvar
        ws{ii, 1} = sdpvar(sum(p_jc(ind(1))), 1);
        vs{ii, 1} = sdpvar(sum(m_jc(ind(1))), 1);
        % Add the unknown sdpvar
        Wi2s(iy, iy) = kron(eye(length(ind)), diag(ws{ii}));
        Vi2s(iu, iu) = kron(eye(length(ind)), diag(vs{ii}));
    end
end
Wi2_val = eye(size(Wi2s));
Vi2_val = eye(size(Vi2s));


%% Optimisation of W and V using fixed D
sdps = [vertcat(ws{:}); vertcat(vs{:})];
sdp_val = 1e-10 * ones(size(sdps));
WV_val = 1e-10;
W_fail = diag(nan(sum(p_jc), 1));
W_fail([indy{I_known}], [indy{I_known}]) = Wi2s([indy{I_known}], [indy{I_known}]);
V_fail = diag(nan(sum(m_jc), 1));
V_fail([indu{I_known}], [indu{I_known}]) = Vi2s([indu{I_known}], [indu{I_known}]);
info.error = 0;
info.varRelError = 0;
info.errloc = "-";

warnState = warning('error', 'YALMIP:SuspectNonSymmetry');
for ii = 1:maxits
    info.iter = ii;

    % Define constraints
    try
        W_constr = [[Wi2s/Drv, N0_f'; N0_f, Vi2s*Dlv] >= 0; sdps <= pmax^2; sdps >= 0];
    catch ME
        if contains(ME.message, 'Suspect non-symmetry')
            % If, due to numerical inaccuracies, the LMI is interpreted as
            % something asymmetric, we use b to scale it and hope to avoid
            % element-wise constraints
            warning(warnState);
            b = sqrt(blkdiag(Wi2_val, Vi2_val));
            W_constr = [b\[Wi2s/Drv, N0_f'; N0_f, Vi2s*Dlv]/b >= 0; sdps <= pmax^2; sdps >= 0];
            warning('error', 'YALMIP:SuspectNonSymmetry');
        else
            throw(ME);
        end
    end

    % Objective function with optional weighting
    WV_obj = trace(aW * Wi2s) + trace(aV * Vi2s);

    % Solve the problem
    WV_sol = optimize(W_constr, WV_obj, options);

    errcode = WV_sol.problem;

    % Analyse error flags
    if errcode ~= 0
        info.error = errcode;
        info.errloc = "W optimisation";
        if ii > 1
            % Break and use values of the previous iteration
            fprintf('An error occured at iteration %d solving for', ii);
            cprintf('*text', ' W/V ');
            fprintf('(errCode %d). Remaining rel. error:', errcode);
            cprintf('*text', ' %.1e.\n', info.relErr);
            break;
        else
            % Return zero values and error
            info.relErr = nan;
            Wi2_val = W_fail; Vi2_val = V_fail;
            cprintf('[1, 0.4, 0]', 'An error occured at iteration %d, solving for', ii);
            cprintf('*[1, 0.4, 0]', ' W/V ');
            cprintf('[1, 0.4, 0]', '(errCode %d).', errcode);
            cprintf('*[1, 0.4, 0]', ' No bounf found.\n');
            break;
        end
    end

    % Stop the loop if the value has converged
    info.relErr = norm((value(WV_obj) - WV_val(end)) ./ WV_val(end));
    info.varRelError = norm((value(sdps) - sdp_val(:, end)) ./ sdp_val(:, end));
    if info.relErr < tol && info.varRelError < sqrt(tol)
        sdp_val(:, ii) = value(sdps);
        fprintf('Bound converged after %2d iterations. Relative error is', ii);
        cprintf('*text', ' %.2e.\n', info.relErr);
        break;
    end

    % Extract the weighting values
    WV_val(ii) = value(WV_obj);
    sdp_val(:, ii) = value(sdps);
    Wi2_val = value(Wi2s);
    Vi2_val = value(Vi2s);

    cla; semilogy(1:length(WV_val), WV_val); xlim([0, maxits]); drawnow;


    %% Optimisation of D using fixed W and V
    % Define constraints
    D_constr = [Dls - sqrt(Vi2_val) \ N0_f * (Wi2_val^-1 * Drs) * N0_f' / sqrt(Vi2_val) >= gamma; gamma >= 0; gamma <= pmax];
    for k = 1:length(R)
        D_constr = [D_constr; R{k} >= 1 / sqrt(pmax) * eye(size(R{k})); R{k} <= sqrt(pmax) * eye(size(R{k}))];
    end

    D_obj = -gamma;

    % Solve the problem
    D_sol = optimize(D_constr, D_obj, options);

    errcode = D_sol.problem;

    % Analyse error flags
    if errcode ~= 0
        info.error = errcode;
        info.errloc = "D optimisation";

        % Break and use values of previous iteration
        fprintf('An error occured at iteration %d solving for', ii - 1);
        cprintf('*text', ' D ');
        fprintf('(errCode %d). Remaining rel. error:', errcode);
        cprintf('*text', ' %.1e.\n', info.relErr);
        break;
    else
        Drv = value(Drs);
        Dlv = value(Dls);
    end
end

warning(warnState);

if ii == maxits && (info.relErr >= tol) || info.varRelError >= sqrt(tol)
    fprintf('Bounds did not converge after %d iterations. Remaining rel. error:', maxits);
    cprintf('*text', ' %.1e.\n', info.relErr);
end

% Copy the original bounds, and replace the unknown ones by the found
% bounds by inverting the weights
for ii = 1:length(p_jc)
    Wout{ii} = Wi2_val(indy{ii}, indy{ii}) ^ (-.5);
    Vout{ii} = Vi2_val(indu{ii}, indu{ii}) ^ (-.5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Wout_inv{ii} = Wi2_val(indy{ii}, indy{ii}) ^ (.5);
    Vout_inv{ii} = Vi2_val(indu{ii}, indu{ii}) ^ (.5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end