%% Structure-preserving model reduction on power systems (using balanced truncation)
%
% In this script, we consider an interconnected power system model of 2
% subsystems, i.e., G1_GFL.mat and G2_powerGrid.mat. The interconnected
% model are reduced in a modular fashion. All subsystems are reduced using
% balanced reuncation and the reduced subsystem models are assembled
% together to construct the reduced interconnected model.
% 
% The reduction steps themselves are controlled using the robust 
% performance framework: An accuracy specification is prescribed for the
% interconnected system and, using robust performance, bounds are
% calculated on each subsystems. Subsystems are then reduced to satisfy
% their associated bounds.
%
% The script is organised as follows:
% - User settings 
% - Translation of the coupled accuracy specification to the low-level
%   specifications
% - Reduction of all the subsystems and coupling
% - Evaluation of the reduced models and output of results
%
% This script is modified from Luuk Poort (TUe) as part of the modular
% model reduction PhD project
% Last updated: 10-09-2025

clc; clear; close all; home;

%% USER SETTINGS %%
%% Some hyper parameters
% Frequency range of interest
w = logspace(-1,4,1000); % in rad/s, original (1, 5, 5)
nw = length(w);

% Elimination method of balanced truncation: "Truncate" or "MatchDC"
elimMethod = "MatchDC";

%% Construct full order model of the interconnected system
% Load subsystem models
%%%%%%%%%%%%%%%%%%%%% Including the Study Area or Not %%%%%%%%%%%%%%%%%%%%%
load("linSys1.mat"); load("linSys2.mat"); % External area
load("linSys3.mat"); % Study area 
% Connect subsystem model as a diagonal transfer matrix
G = {linSys1; linSys2; linSys3};
%%%%%%%%%%%%%%%%%%%%% Including the Study Area or Not %%%%%%%%%%%%%%%%%%%%%

% Generate the subsystem models sys_j, interface model sys_I, and the
% interconnected model sys_c
[sys_j, sys_I, sys_c] = IBR_modelGenerator(G); % !!! Need to be changed for different numbers of subsystems !!!
% sys_c.s = linsysC;

% Number of the order and unstable modes of each subsystem
n_j = arrayfun(@(x) length(x.s.A), sys_j);
nns = arrayfun(@(x) length(x.Gns.A), sys_j);

% Number of subsystems Nc and the order of the interconnected system Ns
Nc = length(sys_j);
Ns = sum(n_j);

% Number of inputs m and outputs p per system
m_j = arrayfun(@(x) length(x.s.u), sys_j);
p_j = arrayfun(@(x) length(x.s.y), sys_j);
m_c = length(sys_c.s.u);
p_c = length(sys_c.s.y);

% Frequency response G(jw) at selected frequency points w
H_0 = freqresp(sys_c.s, w); % interconnected system G_c(jw)
H_B = freqresp(blkdiag(sys_j.s), w); % subsystems in diagonal form G_b(jw)
H_K = freqresp(blkdiag(sys_I.s), w); % subsystems in diagonal form K(jw)
% H_B = freqresp(blkdiag(sys_j(1).Gs, sys_j(2).s, sys_j(3).s), w); % subsystems in diagonal form G_b(jw)
for ic = 1:Nc
    H_j{ic}  = freqresp(sys_j(ic).s, w); % individual subsystem G_j(jw)
    H_sj{ic} = freqresp(sys_j(ic).Gs, w); % stable part of individual subsystems G_j(jw)
end
% H_2_update  = freqresp(linSys2_update, w); % individual subsystem G_j(jw)

%% Given error bound specification on the interconnected system model
% Error bound specified for the interconnected system, i.e.,
% frequency-weighted error ||V * E_c * W|| <= eb_c, where E_c = \hat G_c -
% G_c
eb_c = 0.25 * sys_c.s; % Original 0.25 * sys_c.s + 5e-2
eb_c_resp = freqresp(eb_c, w);      % ny-by-nu-by-nw array
eb_cr = pagenorm(eb_c_resp, 2) + 1e-2; % 1e-2
% eb_cr = pageNorm(freqresp(eb_c, w));
% eb_cr = pageNorm(freqresp(eb_c, w)) + 1e-2; % H_inf norm of eb_c (the error bound system) at selected frequency points w
% eb_cr = max(pageNorm(freqresp(eb_c, w)), 5e-4); % 1e-4

% Error bound specification pn the interconnected system, converted to the
% form ||V_c * E_c * W_c|| <= 1
W_c = ss(1); V_c = inv(eb_c);
W_cr  = repmat(eye(size(V_c.OutputName,1), size(V_c.InputName,1)), 1, 1, nw); % W_c(jw), input weighted matrix at w
V_cr  = 1 ./ eb_cr; % Original: freqresp(V_c, w); % V_c(jw), output weighted matrix at w

% Split the input and output indices into 1) related to the interconnected
% system and 2) related to the subsystems
o1 = 1:(size(sys_I.s, 1) - p_c); % input indices related to subsystems in diagonal form G_b, u_b
o2 = (size(sys_I.s, 1) - p_c + 1):size(sys_I.s, 1); % output indices related to interconnected system G_c, y_c 
i1 = 1:(size(sys_I.s, 2) - m_c); % output indices related to subsystems in diagonal form G_b, y_b
i2 = (size(sys_I.s, 2) - m_c + 1):size(sys_I.s, 2); % input indices related to interconnected system G_c, u_c 

% Define Identity matrices for later use
Io = eye(sum(p_j));
Ii = eye(sum(m_j));

%% TRANSLATION OF THE COUPLED ACCURACY SPECIFICATION TO THE LOW-LEVEL SPECIFICATIONS %%
%% Combine the weighting functions, i.e., generate W(s) and V(s)
Wr = zeros(m_c+sum(p_j), m_c+sum(p_j), nw);
Vr = zeros(p_c+sum(m_j), p_c+sum(m_j), nw);
for ii = 1:nw
    Wr(:, :, ii) = blkdiag(eye(sum(p_j)), W_cr(:, :, ii)); % W(jw), dim = (p_j + m_c) x (p_j + m_c)
    Vr(:, :, ii) = blkdiag(eye(sum(m_j)), V_cr(:, :, ii)); % V(jw), dim = (m_j + p_c) x (m_j + p_c)
end

% ----------- INFO (retained annotation from Luuk's code) --------------- %
% This is the tricky part.
% The large parfor loop calculates subsystem bounds (defined by Wb and
% Vb) that guarantee that the approximate assembly model satisfies its
% bounds (defined by W_c and V_c). These bounds are optimized in a square
% root sense (see also Lars Janssen's IFAC paper). 
% The optimization is performed iteratively, using yalmip. However, this
% optimization crashes sometimes or simply does not converge well. To
% improve the convergence, I added an additional loop ("for k = 1:.." ...),
% to use the 5 different scaling factors contained in mst. This helps  
% because the "optimal" scaling is different per frequency point. It is a 
% very primitive approach however and should be improved...
% ----------------------------- END INFO -------------------------------- %

%% Initialization for the optimization
% Informative struct to keep track of the convergence
convInfo = repmat(struct('error', [], 'errloc', [], 'iter', ...
    [], 'relErr', [], 'varRelError', []), nw, 1);

% Additional scaling factors
ms_k = [1, 3, 10, 30, 100]; % [1, 3, 10, 30, 100]

% Cell arrays to be filled with ALL (Nc subsystems + 1 interconnected
% system) bounds calculated
Wb_a = cell(Nc + 1, length(ms_k), nw); 
Vb_a = cell(Nc + 1, length(ms_k), nw);

%% Calculate bounds for each frequency point using parallelization
p4ParForLoop = gcp('nocreate');
if isempty(p4ParForLoop)
    parpool;
end  

tic
parfor ii = 1:nw
% ii = 3;
    % Cell arrays to be filled with bounds for this frequency point
    Wb_k = cell(Nc + 1, length(ms_k));
    Vb_k = cell(Nc + 1, length(ms_k));
    Wb_inv_k = cell(Nc + 1, length(ms_k));
    Vb_inv_k = cell(Nc + 1, length(ms_k));

    % Informative struct for this frequency point
    convInfo_k = repmat(struct('error', [], 'errloc', [], 'iter', [], ...
        'relErr', [], 'varRelError', []), length(ms_k), 1);

    % Loop over the various scaling factors ms_k (parfor does not allow the
    % '5' to be replace by length(ms_k))
    for k = 1 % k = 1:5
        % Extract and scale the weighting function V(jw)
        Vri = Vr(:,:,ii);
        Vri(o2, i2) = Vri(o2, i2)/ms_k(k); % !!! Need to be changed for MIMO case !!!

        % Compose the nominal system using current scaling
        Si = H_K(o1,i1,ii) / ms_k(k); % sys_I.s.D(o1,i1) / ms_k(k)
        SBr = H_B(:,:,ii) * ms_k(k); % G_B(jw) .* k
        Z = (Io - SBr * Si); % (I - G_B(jw) * S_{BB})
        Nr = Vri * ...
            full([Si/Z, (Ii-Si*SBr)\H_K(o1,i2,ii);  H_K(o2,i1,ii)/Z, zeros(size(o2,2), size(i2,2))]) * ...
            Wr(:,:,ii); % !!! Why V(jw) * N(jw) * W(jw) !!!

        % Calculate the bounds, i.e., V_j(s) and W_j(s) for each subsystems
        try
            [~, Wb_k(:,k),Vb_k(:,k),convInfo_k(k), Wb_inv_k(:,k), Vb_inv_k(:,k)] = ...
                DK_iteration_adapt(Nr,[m_j,p_c],[p_j,m_c],{-1, -2, 1}'); % 1e-15*eye(p_j(3),m_j(3))
        catch ME
            % Error handling
            if contains(ME.message,"non-symmetry")
                % A typical error: non-symmetry of the constraint (due to
                % numerical errors)
                Wb_k(:,k) = num2cell(zeros(3,1)); % !!! Need to be changed for different numbers of subsystems !!!
                Vb_k(:,k) = num2cell(zeros(3,1));
                Wb_inv_k(:,k) = num2cell(zeros(3,1));
                Vb_inv_k(:,k) = num2cell(zeros(3,1));
                convInfo_k(k).error = 99;
                convInfo_k(k).relErr = inf;
            else
                throw(ME);
            end
        end

        % Check for error occurence
        if convInfo_k(k).error ~= 0
            convInfo_k(k).relErr = inf;
        end
        
        % Save data of this frequency point
        Wb_a(:,k,ii) = Wb_k(:,k);
        Vb_a(:,k,ii) = Vb_k(:,k);
            
        % Stop trying if properly converged
        if convInfo_k(k).error == 0 && convInfo_k(k).relErr<1e-4
            break;
        end
    end

    % Keep the best-converged result and scale back
    [~,im] = min([convInfo_k.relErr]); % relErr
    for l = 1:2 % !!! Need to be changed for different numbers of subsystems !!!
        Wb{l,ii} = Wb_k{l,im}/sqrt(ms_k(im));
        Vb{l,ii} = Vb_k{l,im}/sqrt(ms_k(im));
        Wb_inv{l,ii} = Wb_inv_k{l,im}*sqrt(ms_k(im));
        Vb_inv{l,ii} = Vb_inv_k{l,im}*sqrt(ms_k(im));
    end
    convInfo(ii) = convInfo_k(im);
    ms2(ii) = ms_k(im);
end
toc

% Save all data for the Robust SubSystem (RSS) bounds in structs
cInfo.rss = convInfo;
for ii = 1:nw
    Wb{3,ii} = eye(p_c, m_c);   Vb{3,ii} = eye(p_c, m_c); % !!! Why 3? Why all equal 1? !!!
end
WTDa.rss = Wb_a;
VTDa.rss = Vb_a;
WTD.rss = Wb;
VTD.rss = Vb;
WTD.inv = Wb_inv;
VTD.inv = Vb_inv;
vld.rss = ([convInfo.error] == 0 & [convInfo.relErr]<1e-4) | [convInfo.iter]==50;
ms2s.rss = ms2;


%% REDUCTION OF ALL THE SUBSYSTEMS AND COUPLING %%

% Here, we define the G^j_u(s) and G^j_y(s) weighting matrices. These can 
% be full transfer functions (to use frequency-weighted balanced
% truncation, or just matrices, to weigh certain input-output pairs more 
% strongly. Their selection can be based on the IO-based bounds, defined by
% V and W, or you simply do not using weighting at all.
% 
% Frequency weighting (by estimating a model based on V and W) sometimes
% results in very large/small eigenvalues, which prevents Gramian
% calculation. Therefore, we simply use no weighting in this example, to
% ensure robustness.

% No weighting: identity matrices
for ic = 1:Nc
    Gy{ic} = ss(eye(p_j(ic)),'Name','Gy');
    Gu{ic} = ss(eye(m_j(ic)),'Name','Gu');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighting matrices
for ic = 1:Nc
    tic;
    W_o = [];
    W_i = [];

    for dd = 1:length(sys_j(ic).s.D)
        W_o_data = zeros(length(find(vld.rss)), 1);
        W_i_data = zeros(length(find(vld.rss)), 1);
        for ii = find(vld.rss)
            W_o_data_w = WTD.inv{ic, ii};
            W_i_data_w = VTD.inv{ic, ii};

            W_o_data(ii) = W_o_data_w(dd, dd);
            W_i_data(ii) = W_i_data_w(dd, dd);
        end

        % % 3. Smooth the data
        % window_size = 15; % Number of points to look at simultaneously
        % 
        % W_o_data = movmax(W_o_data, window_size);
        % W_i_data = movmax(W_i_data, window_size);
        % 
        % 
        % minWeight = 1e0;
        % maxWeight = 1e1;
        % z1 = log(W_o_data);
        % z2 = log(W_i_data);
        % 
        % % Remove isolated spikes
        % z1 = smoothdata(z1, 'movmedian', 5);
        % z2 = smoothdata(z2, 'movmedian', 5);
        % 
        % % Smooth the trend
        % z1 = smoothdata(z1, 'sgolay', 25);
        % z2 = smoothdata(z2, 'sgolay', 25);
        % 
        % W_o_data = exp(z1);
        % W_i_data = exp(z1);
        % W_o_data = min(max( W_o_data, minWeight), maxWeight);
        % W_i_data = min(max( W_i_data, minWeight), maxWeight);
        % 
        % W_o_data = smoothdata(W_o_data, 'sgolay', 11, 'Degree', 3); %%% 'gaussian'
        % W_i_data = smoothdata(W_i_data, 'sgolay', 11, 'Degree', 3); %%% 'gaussian'


        fitOrder = 3;

        W_o_frd = frd(W_o_data, w);
        W_i_frd = frd(W_i_data, w);

        C2o.UpperBound = [];
        C2o.LowerBound = W_o_frd;

        C2i.UpperBound = [];
        C2i.LowerBound = W_i_frd;

        W_o_siso = fitmagfrd(W_o_frd, fitOrder, 0, [], C2o); %%% fitmagfrd
        W_i_siso = fitmagfrd(W_i_frd, fitOrder, 0, [], C2i); %%% fitmagfrd


        if isempty(W_o)
            W_o = W_o_siso;
            W_i = W_i_siso;
        else
            % append() adds the new SISO system onto the diagonal
            W_o = append(W_o, W_o_siso);
            W_i = append(W_i, W_i_siso);
        end
    end
    
    Gy{ic} = W_o;
    Gu{ic} = W_i;
    toc;
end

% %%
% % Plotting to verify
% for ic = 1 %%% 1:2
%     for dd = 1:8
%         for ii = find(vld.rss)
%             W_o_data_w = WTD.inv{ic, ii};
%             W_i_data_w = VTD.inv{ic, ii};
% 
%             W_o_data(ii) = W_o_data_w(dd, dd);
%             W_i_data(ii) = W_i_data_w(dd, dd);
%         end
%         W_o_siso = Gy{ic};
%         W_o_siso = W_o_siso(dd, dd);
% 
%         W_i_siso = Gu{ic};
%         W_i_siso = W_i_siso(dd, dd);
% 
%         figure;
%         loglog(w, W_i_data, 'b-', 'LineWidth', 2); hold on;
%         [mag_fit, ~, ~] = bode(W_i_siso, w);
%         loglog(w, squeeze(mag_fit), 'r--', 'LineWidth', 2);
%         grid on;
%         legend('Computed Data |V_j|', 'Min-Phase Estimate |\hat{V}_j|');
%         title('Magnitude Fitting (fitmagfrd)');
%         xlabel('Frequency (rad/s)'); ylabel('Magnitude');
% 
%         figure;
%         loglog(w, W_o_data, 'b-', 'LineWidth', 2); hold on;
%         [mag_fit, ~, ~] = bode(W_o_siso, w);
%         loglog(w, squeeze(mag_fit), 'r--', 'LineWidth', 2);
%         grid on;
%         legend('Computed Data |W_j|', 'Min-Phase Estimate |\hat{W}_j|');
%         title('Magnitude Fitting (fitmagfrd)');
%         xlabel('Frequency (rad/s)'); ylabel('Magnitude');
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subsystem reduction using balanced truncation
for ic = 1:Nc
    R = reducespec(sys_j(ic).s, "balanced");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gy{ic} = ss(eye(p_j(ic)),'Name','Gy');
    % Gu{ic} = ss(eye(m_j(ic)),'Name','Gu');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R.Options.InputWeight  = Gu{ic};
    R.Options.OutputWeight = Gy{ic};
    
    tic;
    for ra = 1:n_j(ic)
        sysr = getrom(R, Order = ra, Method="truncate"); %%% Method="truncate"
        % Calculate the reduced model's frequency response
        H_r = freqresp(sysr, w);

        % Calculate the weighted error per frequency point, i.e.,
        % ||W^{-1}_j(jw) * E_j(jw) * V^{-1}_j(jw)||
        eEb = zeros(1, nw);
        for ii = find(vld.rss)
            % eEb(ii) = norm(WTD.rss{ic, ii} \ (H_r(:, :, ii) - H_j{ic}(:, :, ii)) /  VTD.rss{ic, ii});
            eEb(ii) = norm(WTD.inv{ic, ii} * (H_r(:, :, ii) - H_j{ic}(:, :, ii)) *  VTD.inv{ic, ii});
        end

        % Stop if all the error agrees with the bounds at all frequencies
        if all(eEb < 1)
            break;
        end
    end
    toc;
    % For the satisfactory reduced order stable part, add the unstable part
    % NOTE: add this after check the bound, because the unstable part makes
    % checking the bounds very inaccurate (due to the -2 slope created).
    % !!! Why is the case? !!!

    % Combine the reduced stable part and unstable part by summing them
    sys_j(ic).rss = sysr;
    sys_j(ic).stability = isstable(sysr);

    % Frequency response of the reduced subsystem
    Hr_j(ic).rss = freqresp(sys_j(ic).rss, w);

    % Add meta data
    sys_j(ic).rss.Name = strcat(sys_j(ic).s.Name,"_rss");
    sys_j(ic).rss.u = sys_j(ic).s.u;
    sys_j(ic).rss.y = sys_j(ic).s.y;
end
sys_j_BT = sys_j;
save sys_j_BT sys_j_BT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Subsystem reduction using IRKA
% for ic = 1:Nc
%     % Recommended: IRKA assumes stable full model for H2 reduction
%     fullPoles = eig(full(sys_j(ic).s.A));
%     if any(real(fullPoles) >= 0)
%         warning('A has unstable or marginally stable poles. Standard H2-IRKA assumes stable systems.');
%     end
% 
%     for ra = 1:n_j(ic)
% 
%         opts = struct();
%         opts.maxit = 100;
%         opts.tol = 1e-8;
%         opts.verbose = true;
%         opts.realBasis = true;
% 
%         [sysr, ~] = mimoIRKA(sys_j(ic).s.A, sys_j(ic).s.B, sys_j(ic).s.C, sys_j(ic).s.D, ra, opts);
% 
%         % Calculate the reduced model's frequency response
%         H_r = freqresp(sysr, w);
% 
%         % Calculate the weighted error per frequency point, i.e.,
%         % ||W^{-1}_j(jw) * E_j(jw) * V^{-1}_j(jw)||
%         eEb = zeros(1, nw);
%         for ii = find(vld.rss)
%             % eEb(ii) = norm(WTD.rss{ic, ii} \ (H_r(:, :, ii) - H_j{ic}(:, :, ii)) /  VTD.rss{ic, ii});
%             eEb(ii) = norm(WTD.inv{ic, ii} * (H_r(:, :, ii) - H_j{ic}(:, :, ii)) *  VTD.inv{ic, ii});
%         end
% 
%         % Stop if all the error agrees with the bounds at all frequencies
%         if all(eEb <= 1)
%             break;
%         end
%     end
%     % For the satisfactory reduced order stable part, add the unstable part
%     % NOTE: add this after check the bound, because the unstable part makes
%     % checking the bounds very inaccurate (due to the -2 slope created).
%     % !!! Why is the case? !!!
% 
%     % Combine the reduced stable part and unstable part by summing them
%     sys_j(ic).rss = sysr;
% 
%     % Frequency response of the reduced subsystem
%     Hr_j(ic).rss = freqresp(sys_j(ic).rss, w);
% 
%     % Add meta data
%     sys_j(ic).rss.Name = strcat(sys_j(ic).s.Name,"_rss");
%     sys_j(ic).rss.u = sys_j(ic).s.u;
%     sys_j(ic).rss.y = sys_j(ic).s.y;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Subsystem reduction using Loewner + AAA
for ic = 1:Nc
    % demo_mimoLoewnerWeightedAAA;
    demo_mimoLoewnerWeightedAAA_stableStrict;

    % Combine the reduced stable part and unstable part by summing them
    sys_j(ic).rss = model.sys;
    sys_j(ic).stability = isstable(model.sys);

    % Frequency response of the reduced subsystem
    Hr_j(ic).rss = freqresp(sys_j(ic).rss, w);

    % Add meta data
    sys_j(ic).rss.Name = strcat(sys_j(ic).s.Name,"_rss");
    sys_j(ic).rss.u = sys_j(ic).s.u;
    sys_j(ic).rss.y = sys_j(ic).s.y;
end
sys_j_AA = sys_j;
save sys_j_AA sys_j_AA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construction of the reduced interconnected system
sys_j(1) = sys_j_BT(1);
sys_j(2) = sys_j_AA(2);

sys_c.rss = lft(blkdiag(sys_j.rss), sys_I.s);
sys_c.rss.Name = 'RSS';

sys_c.stability =  isstable(sys_c.rss);

fprintf('\nStability of interconnection: %d\n', sys_c.stability);


%% EVALUATION OF THE REDUCED MODELS AND OUTPUT OF RESULTS %%
% Extract the names of the performed reduction methods, which allows to 
% skip some of the three reduction approaches. However, in this modified
% code there is only one approach retained
redNms = string(fieldnames(sys_c).');
redNms = setdiff(redNms,["m","s", "stability"]);

% Calculate the frequency responses of the reduced interconnected system
% and the associated error system
for nm = redNms
    Hr.(nm) = freqresp(sys_c.(nm), w); % \hat G_c(jw)
    Herr.(nm) = H_0 - Hr.(nm); % E_c(jw)
end


% %% Bode plot of FOM and mROM
% % FOM
% figure();
% [mag,~,wout] = bode(sys_c.s);
% magdb = 20*log10(mag);
% semilogx(wout, squeeze(magdb), 'b-', 'LineWidth', 1.6);
% grid on;
% hold on;
% % mROM
% [mag,~,wout] = bode(sys_c.rss);
% magdb = 20*log10(mag);
% semilogx(wout, squeeze(magdb), 'r--', 'LineWidth', 1.6);
% xlim([w(1), w(end)]);
% xlabel('Frequency (rad/s)', 'interpreter', 'latex', 'fontSize', 13.5);
% ylabel('Magnitude (dB)', 'interpreter', 'latex', 'fontSize', 13.5);
% 
% legend("FOM", "mROM", 'interpreter', 'latex', 'fontSize', 13.5);
% % Min and max limit of the plot
% [mag,~,~] = bode(sys_c.s, w);
% magdb = 20*log10(mag);
% yl = ylim;
% min_magdb = yl(1);
% max_magdb = yl(end);


% %% Plot of the interconnected system ROM, FOM, and specs
% 
% % Add additional bounds to the frequency grid for patch plotting
% fe = [1e-10, w, 1e10];
% 
% % Prespecified bounds for patch plotting:
% % gbm represents the system +- the bounds
% % ebm represents only the bounds itself
% e_j = zeros(p_c, 1); e_j(oCn) = 1;
% eb_cr_ji = zeros(1, nw);
% for jj = 1:1:nw
%     eb_cr_ji(jj) = norm(freqresp(V_c, w(jj)) * e_j);
% end
% 
% gbm  = max([eps, pageNorm(H_0(oCn, iCn, :)) - eb_cr_ji, eps, eps, fliplr(pageNorm(H_0(oCn, iCn, :)) + eb_cr_ji), eps], eps);
% ebm  = [eps, eb_cr_ji, eps * ones(1, numel(fe) + 1)];
% 
% % Interconnected system G_c, \hat G_c
% figure();
% ax.Gc = axes('Layer', 'top');
% hold on;
% patch([fe fliplr(fe)], 20 * log10(gbm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
% semilogx(w, 20 * log10(pageNorm(H_0(oCn, iCn, :))), 'b-', 'linewidth', 1.6);
% semilogx(w, 20 * log10(pageNorm(Hr.rss(oCn, iCn, :))), 'r--', 'linewidth', 1.6);
% 
% leg_1 = 'Acc. spec.';
% leg_2 = sprintf('FOM from %d to %d', iCn, oCn); 
% leg_3 = sprintf('ROM from %d to %d', iCn, oCn);
% lgdc = legend(leg_1, leg_2, leg_3);
% xlabel(gca, 'Frequency (rad/s)', 'interpreter', 'latex');
% ylabel('Magnitude (dB)', 'interpreter', 'latex');
% xlim(gca, [w(1), w(end)]);
% % ylim(gca, [min_magdb, max_magdb]); % !!! needs to be changed according to different interconnected systems !!!
% % xticklabels({});
% grid on;
% set(gca, 'XScale', 'log');
% 
% % Interconnected system error E_c
% figure();
% ax.Ec = axes('Layer', 'top');
% hold on;
% patch([fe fliplr(fe)], 20 * log10(ebm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
% semilogx(w, 20 * log10(pageNorm(H_0(oCn, iCn, :) - Hr.rss(oCn, iCn, :))), 'r--', 'linewidth', 1.6);
% 
% lgdce = legend('Acc. spec.', sprintf('Err. sys. from %d to %d', iCn, oCn));
% xlabel(gca, 'Frequency (rad/s)', 'interpreter', 'latex');
% ylabel(ax.Ec, 'Err. magn. (dB)', 'interpreter', 'latex');
% xlim(gca, [w(1), w(end)]);
% % ylim(gca, [db(eps), 0]); % !!! needs to be changed according to different interconnected systems !!!
% % set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
% grid on;

% %% Plot of the subsystem ROM, FOM, and specs
% 
% % Select the subsystem to plot
% ic = 2;
% 
% % Bounds on (\hat G_j - G_j)
% for ii = 1:nw
%     % eb_VW(ic, ii) = norm(WTD.rss{ic, ii} * ones(length(WTD.rss{ic, ii}), length(VTD.rss{ic, ii})) * VTD.rss{ic, ii}) ...
%         % / norm(ones(length(WTD.rss{ic, ii}), length(VTD.rss{ic, ii}))); % !!! Why is the case? Could make the error more consevertive. !!!
%     eb_VW(ic, ii) = norm(WTD.rss{ic, ii}) * norm(VTD.rss{ic, ii});
% end
% 
% % Prespecified bounds for patch plotting:
% % gbm represents the system +- the bounds
% % ebm represents only the bounds itself
% 
% oCn = 2;
% iCn = 1;
% 
% e_jr = zeros(p_j(ic), 1); e_jr(oCn) = 1; e_ir = zeros(m_j(ic), 1); e_ir(iCn) = 1;
% eb_VW_ji = zeros(1, nw);
% for jj = 1:1:nw
%     eb_VW_ji(jj) = norm(WTD.rss{ic, jj} * e_jr) * norm(VTD.rss{ic, jj} * e_ir);
% end
% 
% gbm  = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_ji, eps, eps, fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_ji), eps], eps); % Original: pageNorm(H_j{ic}) - eb_VW(ic, :)
% ebm  = [eps, eb_VW_ji, eps * ones(1, numel(fe) + 1)]; % Original: eb_VW(ic, :)
% 
% % Subsystems G_j, \hat G_j
% figure();
% ax.Gj = axes('Layer','top');
% hold on;
% patch([fe fliplr(fe)], 20 * log10(gbm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
% semilogx(w, 20 * log10(pageNorm(H_j{ic}(oCn, iCn, :))), 'b-', 'linewidth', 1.6);
% semilogx(w, 20 * log10(pageNorm(Hr_j(ic).rss(oCn, iCn, :))), 'r--', 'linewidth', 1.6);
% 
% leg_sub1 = 'Acc. spec.';
% leg_sub2 = sprintf('sub_FOM from %d to %d', iCn, oCn); 
% leg_sub3 = sprintf('sub_ROM from %d to %d', iCn, oCn);
% lgdj = legend(leg_sub1, leg_sub2, leg_sub3);
% % lgdj = legend('Acc. spec.', 'FOM', 'RSS');
% xlabel(gca, 'Frequency (rad/s)', 'interpreter', 'latex');
% ylabel(ax.Gj, 'Magnitude (dB)', 'interpreter', 'latex');
% xlim(gca, [w(1), w(end)]);
% % ylim(gca, [db(eps), 50]); % !!! needs to be changed according to different interconnected systems !!!
% grid on;
% set(gca, 'XScale', 'log');
% 
% % Subsystems error E_j
% figure();
% ax.Ej = axes('Layer','top');
% hold on;
% patch([fe fliplr(fe)], 20 * log10(ebm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
% semilogx(w, 20 * log10(pageNorm(H_j{ic}(oCn, iCn, :) - Hr_j(ic).rss(oCn, iCn, :))), 'r--', 'linewidth', 1.6);
% 
% lgdje = legend('Acc. spec.', sprintf('Err. sys. from %d to %d', iCn, oCn));
% xlabel(gca, 'Frequency (rad/s)', 'interpreter', 'latex');
% ylabel(ax.Ej, 'Err. magn. (dB)', 'interpreter', 'latex');
% xlim(gca, [w(1), w(end)]);
% % ylim(gca, [db(eps), 0]); % !!! needs to be changed according to different interconnected systems !!!
% grid on;
% set(gca, 'XScale', 'log');


% %% Save ROM of subsystems (\hat G_j) for time-domain simulation
% area_2_rom = sys_j(1).rss;
% area_3_rom = sys_j(2).rss;
% save area_2_rom area_2_rom;
% save area_3_rom area_3_rom;

%% Plot of the interconnected system ROM, FOM, and specs

% Add additional bounds to the frequency grid for patch plotting
fe = [1e-10, w, 1e10];

% Prespecified bounds for patch plotting:
% gbm represents the system +- the bounds
% ebm represents only the bounds itself
%%
eb_cr = squeeze(eb_cr)'; %%
%%
gbm  = max([eps, pageNorm(H_0) - eb_cr, eps, eps, fliplr(pageNorm(H_0) + eb_cr), eps], eps);
ebm  = [eps, eb_cr, eps * ones(1, numel(fe) + 1)];

% Interconnected system G_c, \hat G_c
figure();
ax.Gc = subplot(2, 1, 1);
set(ax.Gc, 'Layer', 'top');
hold(ax.Gc, 'on');
patch(ax.Gc, [fe fliplr(fe)], db(gbm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
semilogx(ax.Gc, w, db(pageNorm(H_0)), 'b-', 'linewidth', 1.6);
semilogx(ax.Gc, w, db(pageNorm(Hr.rss)), 'r--', 'linewidth', 1.6);

lgdc = legend(ax.Gc, 'Acc. spec.', 'FOM', 'ROM');
xlabel(ax.Gc, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.Gc, 'Magnitude (dB)', 'interpreter', 'latex');
xlim(ax.Gc, [w(1), w(end)]);
% ylim(gca, [min_magdb, max_magdb]); % !!! needs to be changed according to different interconnected systems !!!
% xticklabels({});
grid(ax.Gc, 'on');
set(ax.Gc, 'XScale', 'log');

% % Interconnected system error E_c
% figure();
% ax.Ec = axes('Layer', 'top');
% hold on;
% patch([fe fliplr(fe)], db(ebm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
% semilogx(w, db(pageNorm(H_0 - Hr.rss)), 'r--', 'linewidth', 1.6);
% 
% lgdce = legend('Acc. spec.', 'Err. sys.');
% xlabel(gca, 'Frequency (rad/s)', 'interpreter', 'latex');
% ylabel(ax.Ec, 'Err. magn. (dB)', 'interpreter', 'latex');
% xlim(gca, [w(1), w(end)]);
% % ylim(gca, [db(eps), 0]); % !!! needs to be changed according to different interconnected systems !!!
% % set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
% grid on;

% Normalised interconnected system error E_c
eEc = zeros(1, nw);
for ii = find(vld.rss)
    eEc(ii) = norm(V_cr(ii) * (H_0(:, :, ii) - Hr.rss(:, :, ii)));
end

% figure();
ax.nEc = subplot(2, 1, 2);
set(ax.nEc, 'Layer', 'top');
hold(ax.nEc, 'on');
patch(ax.nEc, [fe fliplr(fe)], [ones(size(fe)), fliplr(eps*ones(size(fe)))], [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
semilogx(ax.nEc, w, eEc, 'r--', 'linewidth', 1.6);

lgdje = legend(ax.nEc, 'Acc. spec.', 'Err. sys.');
xlabel(ax.nEc, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.nEc, 'Err. magn. (dB)', 'interpreter', 'latex');
xlim(ax.nEc, [1e-1, 1e4]);
ylim(ax.nEc, [1e-7*0.3, 1e0*1.3]); % !!! needs to be changed according to different interconnected systems !!!
grid(ax.nEc, 'on');
set(ax.nEc, 'XScale', 'log');
set(ax.nEc, 'YScale', 'log');
%%
fig = figure();
clf(fig);

tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---------- Top tile ----------
ax.Gc = nexttile;

set(ax.Gc, 'Layer', 'top');
hold(ax.Gc, 'on');

patch(ax.Gc, ...
    [fe fliplr(fe)], ...
    db(gbm), ...
    [.8 .8 1], ...
    'LineWidth', 1.6, ...
    'LineStyle', ':', ...
    'EdgeColor', [0 0 0]);

semilogx(ax.Gc, w, db(pageNorm(H_0)), 'b-', 'linewidth', 2.0);
semilogx(ax.Gc, w, db(pageNorm(Hr.rss)), 'r--', 'linewidth', 2.0);

legend(ax.Gc, 'Red. err. spec.', 'FOM', 'ROM', 'interpreter', 'latex');
% xlabel(ax.Gc, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.Gc, 'Magnitude (dB)', 'interpreter', 'latex');
xlim(ax.Gc, [1e-1, 1e4]);
ylim(ax.Gc, [-90, 0]);
grid(ax.Gc, 'on');
set(ax.Gc, 'XScale', 'log');
set(ax.Gc, 'XTickLabel', []);

% ---------- Bottom tile ----------
ax.nEc = nexttile;
set(ax.nEc, 'Layer', 'top');
hold(ax.nEc, 'on');

patch(ax.nEc, ...
    [fe fliplr(fe)], ...
    [ones(size(fe)), fliplr(eps * ones(size(fe)))], ...
    [.8 .8 1], ...
    'LineWidth', 1.2, ...
    'LineStyle', ':', ...
    'EdgeColor', [0 0 0]);

loglog(ax.nEc, [1e-1, 1e4], [1, 1], 'k:', 'LineWidth', 2.0);

semilogx(ax.nEc, w, eEc, 'r--', 'linewidth', 2.0);

legend(ax.nEc, 'Red. err. spec.', 'Weigh. red. err.', 'interpreter', 'latex');
xlabel(ax.nEc, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.nEc, 'Frequency-weighted reduction error', 'interpreter', 'latex');
xlim(ax.nEc, [1e-1, 1e4]);
ylim(ax.nEc, [1e-7 * 0.2, 1e0 * 1.3]);
grid(ax.nEc, 'on');
set(ax.nEc, 'XScale', 'log');
set(ax.nEc, 'YScale', 'log');

% Optional figure size for two-column-width publication figure
set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

%% Plot of the subsystem ROM, FOM, and specs

% Select the subsystem to plot
ic = 1;

% Bounds on (\hat G_j - G_j)
for ii = 1:nw
    eb_VW(ic, ii) = norm(WTD.rss{ic, ii} * ones(length(WTD.rss{ic, ii}), length(VTD.rss{ic, ii})) * VTD.rss{ic, ii}) ...
        / norm(ones(length(WTD.rss{ic, ii}), length(VTD.rss{ic, ii}))); % !!! Why is the case? Could make the error more consevertive. !!!
    % eb_VW(ic, ii) = norm(WTD.rss{ic, ii}) * norm(VTD.rss{ic, ii});
end

% Prespecified bounds for patch plotting:
% gbm represents the system +- the bounds
% ebm represents only the bounds itself
gbm  = max([eps, pageNorm(H_j{ic}) - eb_VW(ic, :), eps, eps, fliplr(pageNorm(H_j{ic}) + eb_VW(ic, :)), eps], eps); % Original: pageNorm(H_j{ic}) - eb_VW(ic, :)
ebm  = [eps, eb_VW(ic, :), eps * ones(1, numel(fe) + 1)]; % Original: eb_VW(ic, :)

% Subsystems G_j, \hat G_j
figure();
ax.Gj = axes('Layer','top');
hold on;
patch([fe fliplr(fe)], db(gbm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
semilogx(w, db(pageNorm(H_j{ic})), 'b-', 'linewidth', 1.6);
semilogx(w, db(pageNorm(Hr_j(ic).rss)), 'r--', 'linewidth', 1.6);

lgdj = legend('Acc. spec.', 'FOM', 'ROM');
xlabel(gca, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.Gj, 'Magnitude (dB)', 'interpreter', 'latex');
xlim(gca, [w(1), w(end)]);
% ylim(gca, [db(eps), 50]); % !!! needs to be changed according to different interconnected systems !!!
grid on;
set(gca, 'XScale', 'log');


% Subsystems error E_j
figure();
ax.Ej = axes('Layer','top');
hold on;
patch([fe fliplr(fe)], db(ebm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
semilogx(w, db(pageNorm(H_j{ic} - Hr_j(ic).rss)), 'r--', 'linewidth', 1.6);

lgdje = legend('Acc. spec.', 'Err. sys.');
xlabel(gca, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.Ej, 'Err. magn. (dB)', 'interpreter', 'latex');
xlim(gca, [w(1), w(end)]);
% ylim(gca, [db(eps), 0]); % !!! needs to be changed according to different interconnected systems !!!
grid on;
set(gca, 'XScale', 'log');


% Normalised subsystems error E_j
eEj = zeros(1, nw);
for ii = find(vld.rss)
    eEj(ii) = norm(WTD.inv{ic, ii} * (Hr_j(ic).rss(:, :, ii) - H_j{ic}(:, :, ii)) *  VTD.inv{ic, ii});
end

figure();
ax.nEj = axes('Layer','top');
hold on;
patch([fe fliplr(fe)], [ones(size(fe)), fliplr(eps*ones(size(fe)))], [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);
semilogx(w, eEj, 'r--', 'linewidth', 1.6);

lgdje = legend('Acc. spec.', 'Err. sys.');
xlabel(gca, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.nEj, 'Err. magn. (dB)', 'interpreter', 'latex');
xlim(gca, [w(1), w(end)]);
ylim(gca, [1e-2*0.9, 1e0+0.1]); % !!! needs to be changed according to different interconnected systems !!!
grid on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');


%% Plot of the subsystem ROM, FOM, and specs
% Compact two-panel version, same style as interconnected system plot

% Select the subsystem to plot
ic = 2;

% Force frequency grid to row vector
w = w(:).';

% Add additional bounds to the frequency grid for patch plotting
fe = [1e-10, w, 1e10];

% Subsystem FOM and ROM magnitudes
GjMag  = pageNorm(H_j{ic});
HrjMag = pageNorm(Hr_j(ic).rss);

GjMag  = GjMag(:).';
HrjMag = HrjMag(:).';

% Bounds on (\hat G_j - G_j)
eb_VW = zeros(1, nw);

for ii = 1:nw
    eb_VW(ii) = norm(WTD.rss{ic, ii} * ...
        100*ones(length(WTD.rss{ic, ii}), length(VTD.rss{ic, ii})) * ...
        VTD.rss{ic, ii}) ...
        / norm(ones(length(WTD.rss{ic, ii}), length(VTD.rss{ic, ii})));

    % Alternative:
    eb_VW(ii) = norm(WTD.rss{ic, ii}) * norm(VTD.rss{ic, ii});
end

eb_VW = eb_VW(:).';

% Prespecified bounds for patch plotting
% gbm represents subsystem magnitude +- bounds
gbm = max([ ...
    eps, ...
    GjMag - eb_VW, ...
    eps, ...
    eps, ...
    fliplr(GjMag + eb_VW), ...
    eps], eps);

% Normalised subsystem error E_j
% Use NaN instead of 0 because log-scale y-axis cannot show zero
eEj = NaN(1, nw);

for ii = find(vld.rss)
    eEj(ii) = norm(WTD.inv{ic, ii} * ...
        (Hr_j(ic).rss(:, :, ii) - H_j{ic}(:, :, ii)) * ...
        VTD.inv{ic, ii});
end


% Compact two-panel subsystem plot

figure();
clf;

t = tiledlayout(2, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% ============================================================
% Top tile: subsystem FOM, ROM, and reduction error specification
% ============================================================
ax.Gj = nexttile(t);
set(ax.Gj, 'Layer', 'top');
hold(ax.Gj, 'on');

patch(ax.Gj, ...
    [fe fliplr(fe)], ...
    db(gbm), ...
    [.8 .8 1], ...
    'LineWidth', 1.2, ...
    'LineStyle', ':', ...
    'EdgeColor', [0 0 0]);

semilogx(ax.Gj, w, db(GjMag),  'b-',  'linewidth', 1.6);
semilogx(ax.Gj, w, db(HrjMag), 'r--', 'linewidth', 1.6);

legend(ax.Gj, ...
    'Red. err. spec.', ...
    'Weigh. FOM', ...
    'Weigh. ROM', ...
    'interpreter', 'latex');

ylabel(ax.Gj, 'Magnitude (dB)', 'interpreter', 'latex');

xlim(ax.Gj, [w(1), w(end)]);
ylim(ax.Gj, [-5, 35]);   % adjust manually if needed

grid(ax.Gj, 'on');
set(ax.Gj, 'XScale', 'log');

% Hide top x tick labels for compactness
xticklabels(ax.Gj, []);


% ============================================================
% Bottom tile: frequency-weighted normalised reduction error
% ============================================================
ax.nEj = nexttile(t);
set(ax.nEj, 'Layer', 'top');
hold(ax.nEj, 'on');

patch(ax.nEj, ...
    [fe fliplr(fe)], ...
    [ones(size(fe)), fliplr(eps * ones(size(fe)))], ...
    [.8 .8 1], ...
    'LineWidth', 1.2, ...
    'LineStyle', ':', ...
    'EdgeColor', [0 0 0]);

semilogx(ax.nEj, w, eEj, 'r--', 'linewidth', 1.6);

legend(ax.nEj, ...
    'Red. err. spec.', ...
    'Weigh. red. err.', ...
    'interpreter', 'latex');

xlabel(ax.nEj, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.nEj, 'Frequency-weighted reduction error', ...
    'interpreter', 'latex');

xlim(ax.nEj, [w(1), w(end)]);
ylim(ax.nEj, [1e-3, 1e0 * 1.1]);

grid(ax.nEj, 'on');
set(ax.nEj, 'XScale', 'log');
set(ax.nEj, 'YScale', 'log');

% Link frequency axes
linkaxes([ax.Gj, ax.nEj], 'x');

%% Poles
% % Example systems
% sys_j(1) = sys_j_BT(1);
% sys_j(2) = sys_j_AA(2);
% 
% sys_c.Exernal = lft(blkdiag(sys_j.s), sys_I.s.D);
% sys_c.rssExernal = lft(blkdiag(sys_j.rss), sys_I.s.D);
% sys_c.rssExernal.Name = 'RSSExternal';


sys1 = sys_c.s;
sys2 = sys_c.rss;

% Extract poles
p1 = pole(sys1);
p2 = pole(sys2);

% Find pole closest to imaginary axis for each system
[~, idx1] = min(abs(real(p1)));
[~, idx2] = min(abs(real(p2)));

p1_close = p1(idx1);
p2_close = p2(idx2);

fig = figure();

% Main pole plot
hold on;

plot(real(p1), imag(p1), 'bo', ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2);

plot(real(p2), imag(p2), 'rx', ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2);


% Axes
xline(0, 'k--', 'LineWidth', 1);
yline(0, 'k--', 'LineWidth', 1);

grid on;
box on;
axis equal;

xlabel('Real axis', 'interpreter', 'latex');
ylabel('Imaginary axis', 'interpreter', 'latex');
% title('Pole Plot of Two Systems');

legend('Poles of full power grid model', ...
       'Poles of reduced power grid model', ...
        'interpreter', 'latex'); % 'Location', 'best',

hold off;
% ylim([-2000, 2000]);
% xlim([-2000, 200]);

% Inset position: [left bottom width height]
axes('Position', [0.2 0.58 0.30 0.30]);
hold on;

plot(real(p1), imag(p1), 'bo', ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2);

plot(real(p2), imag(p2), 'rx', ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2);


xline(0, 'k--', 'LineWidth', 1);
yline(0, 'k--', 'LineWidth', 1);

grid on;
box on;

% Zoom around both closest poles
xlim([-200, 20]);
ylim([-500, 500]);

% title('Inset Zoom');
hold off;

% Optional figure size for two-column-width publication figure
set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);
