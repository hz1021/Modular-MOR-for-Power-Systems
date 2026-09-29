%% Modular Model Order Reduction for Power Systems
%
% In this script, we consider an interconnected power system model of 3
% subsystems, i.e., two external areas linSys1.mat and linSys2.mat, and a 
% study area linSys3.mat. The interconnected model are reduced in a modular
% fashion. All subsystems are reduced using either balanced reuncation or 
% greedy Loewner, and the reduced subsystem models are assembled together 
% to construct the reduced interconnected model.
% 
% The reduction steps themselves are controlled using the robust 
% performance framework: An accuracy specification is prescribed for the
% interconnected system and, using robust performance, accuracy 
% specifications are calculated on each subsystems. Subsystems are then 
% reduced to satisfy their associated accuracy specifications.
%
% The script is organised as follows:
% (1) Construct the interconnected power grid, i.e., external areas and
% their interconnection with the study area
% (2) Define the global accuracy specification 
% (3) Translation of the coupled accuracy specification to local accuracy 
% specifications 
% (4) Reduce all the subsystems and reconnect
% (5) Evaluation of the reduced models and output of results
%
% This script is modified from Luuk Poort (TUe) as part of the modular
% model reduction PhD project
% Last updated: 28-09-2026

clc; clear; close all; home;

%% (1) Construct the Interconnected Power Grid %%
%%% Some hyper parameters
% Frequency range of interest
w = logspace(-1,4,1000); % in rad/s
nw = length(w);

%%% Construct full order model of the interconnected system
% Load subsystem models
load("linSys1.mat"); load("linSys2.mat"); % External area
load("linSys3.mat"); % Study area 
% Connect subsystem model as a diagonal transfer matrix
G = {linSys1; linSys2; linSys3};

% Generate the subsystem models sys_j, interface model sys_I, and the
% interconnected model sys_c
[sys_j, sys_I, sys_c] = IBR_modelGenerator(G);

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

for ic = 1:Nc
    H_j{ic}  = freqresp(sys_j(ic).s, w); % individual subsystem G_j(jw)
    H_sj{ic} = freqresp(sys_j(ic).Gs, w); % stable part of individual subsystems G_j(jw)
end

%% (2) Define the Global Accuracy Specification %%
% Error bound specified for the interconnected system, i.e.,
% frequency-weighted error ||V * E_c * W|| <= eb_c, where E_c = \hat G_c -
% G_c
eb_c = 0.25 * sys_c.s;
eb_c_resp = freqresp(eb_c, w); % ny-by-nu-by-nw array
eb_cr = pagenorm(eb_c_resp, 2) + 1e-2;

% Error bound specification pn the interconnected system, converted to the
% form ||V_c * E_c * W_c|| <= 1
W_c = ss(1); V_c = inv(eb_c);
W_cr  = repmat(eye(size(V_c.OutputName,1), size(V_c.InputName,1)), 1, 1, nw); % W_c(jw), input weighted matrix at w
V_cr  = 1 ./ eb_cr; % Original: freqresp(V_c, w); % V_c(jw), output weighted matrix at w

% Split the input and output indices into 1. related to the interconnected
% system and 2. related to the subsystems
o1 = 1:(size(sys_I.s, 1) - p_c); % input indices related to subsystems in diagonal form G_b, u_b
o2 = (size(sys_I.s, 1) - p_c + 1):size(sys_I.s, 1); % output indices related to interconnected system G_c, y_c 
i1 = 1:(size(sys_I.s, 2) - m_c); % output indices related to subsystems in diagonal form G_b, y_b
i2 = (size(sys_I.s, 2) - m_c + 1):size(sys_I.s, 2); % input indices related to interconnected system G_c, u_c 

% Define Identity matrices for later use
Io = eye(sum(p_j));
Ii = eye(sum(m_j));

%% (3) Translation of Accuracy Specification from Global to Local %%
%%% Combine the weighting functions, i.e., generate W(s) and V(s)
Wr = zeros(m_c+sum(p_j), m_c+sum(p_j), nw);
Vr = zeros(p_c+sum(m_j), p_c+sum(m_j), nw);
for ii = 1:nw
    Wr(:, :, ii) = blkdiag(eye(sum(p_j)), W_cr(:, :, ii)); % W(jw), dim = (p_j + m_c) x (p_j + m_c)
    Vr(:, :, ii) = blkdiag(eye(sum(m_j)), V_cr(:, :, ii)); % V(jw), dim = (m_j + p_c) x (m_j + p_c)
end

%%% Initialization for the optimization
% Informative struct to keep track of the convergence
convInfo = repmat(struct('error', [], 'errloc', [], 'iter', ...
    [], 'relErr', [], 'varRelError', []), nw, 1);

% Additional scaling factors
ms_k = 1;

% Cell arrays to be filled with ALL (Nc subsystems + 1 interconnected
% system) bounds calculated
Wb_a = cell(Nc + 1, length(ms_k), nw); 
Vb_a = cell(Nc + 1, length(ms_k), nw);

%%% Calculate bounds for each frequency point using parallelization
p4ParForLoop = gcp('nocreate');
if isempty(p4ParForLoop)
    parpool;
end  

tic
parfor ii = 1:nw
    % Cell arrays to be filled with bounds for this frequency point
    Wb_k = cell(Nc + 1, length(ms_k));
    Vb_k = cell(Nc + 1, length(ms_k));
    Wb_inv_k = cell(Nc + 1, length(ms_k));
    Vb_inv_k = cell(Nc + 1, length(ms_k));

    % Informative struct for this frequency point
    convInfo_k = repmat(struct('error', [], 'errloc', [], 'iter', [], ...
        'relErr', [], 'varRelError', []), length(ms_k), 1);

    for k = 1
        % Extract and scale the weighting function V(jw)
        Vri = Vr(:,:,ii);
        Vri(o2, i2) = Vri(o2, i2)/ms_k(k);

        % Compose the nominal system using current scaling
        Si = H_K(o1,i1,ii) / ms_k(k); % sys_I.s.D(o1,i1) / ms_k(k)
        SBr = H_B(:,:,ii) * ms_k(k); % G_B(jw) .* k
        Z = (Io - SBr * Si); % (I - G_B(jw) * S_{BB})
        Nr = Vri * ...
            full([Si/Z, (Ii-Si*SBr)\H_K(o1,i2,ii);  H_K(o2,i1,ii)/Z, zeros(size(o2,2), size(i2,2))]) * ...
            Wr(:,:,ii);

        % Calculate the bounds, i.e., V_j(s) and W_j(s) for each subsystems
        try
            [~, Wb_k(:,k),Vb_k(:,k),convInfo_k(k), Wb_inv_k(:,k), Vb_inv_k(:,k)] = ...
                DK_iteration_adapt(Nr,[m_j,p_c],[p_j,m_c],{-1, -2, 1}');
        catch ME
            % Error handling
            if contains(ME.message,"non-symmetry")
                % A typical error: non-symmetry of the constraint (due to
                % numerical errors)
                Wb_k(:,k) = num2cell(zeros(3,1));
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
    end

    % Keep the best-converged result and scale back
    [~,im] = min([convInfo_k.relErr]);
    for l = 1:2
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
    Wb{3,ii} = eye(p_c, m_c);   Vb{3,ii} = eye(p_c, m_c);
end
WTDa.rss = Wb_a;
VTDa.rss = Vb_a;
WTD.rss = Wb;
VTD.rss = Vb;
WTD.inv = Wb_inv;
VTD.inv = Vb_inv;
vld.rss = ([convInfo.error] == 0 & [convInfo.relErr]<1e-4) | [convInfo.iter]==50;
ms2s.rss = ms2;


%% (4) Reduce All the Subsystems and Reconnect %%

% Here, we define the G^j_u(s) and G^j_y(s) weighting matrices. These can 
% be full transfer functions (to use frequency-weighted balanced
% truncation, or just matrices, to weigh certain input-output pairs more 
% strongly. Their selection can be based on the IO-based bounds, defined by
% V and W, or you simply do not using weighting at all.

% No weighting: identity matrices
for ic = 1:Nc
    Gy{ic} = ss(eye(p_j(ic)),'Name','Gy');
    Gu{ic} = ss(eye(m_j(ic)),'Name','Gu');
end


% With weighting: fit the wrighting matrices
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

        fitOrder = 3;

        W_o_frd = frd(W_o_data, w);
        W_i_frd = frd(W_i_data, w);

        C2o.UpperBound = [];
        C2o.LowerBound = W_o_frd;

        C2i.UpperBound = [];
        C2i.LowerBound = W_i_frd;

        W_o_siso = fitmagfrd(W_o_frd, fitOrder, 0, [], C2o);
        W_i_siso = fitmagfrd(W_i_frd, fitOrder, 0, [], C2i);


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


%%% Subsystem reduction using frequency weighted balanced truncation
for ic = 1:Nc
    R = reducespec(sys_j(ic).s, "balanced");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Unannotate this block if using balanced truncation
    % Gy{ic} = ss(eye(p_j(ic)),'Name','Gy');
    % Gu{ic} = ss(eye(m_j(ic)),'Name','Gu');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R.Options.InputWeight  = Gu{ic};
    R.Options.OutputWeight = Gy{ic};
    
    tic;
    for ra = 1:n_j(ic)
        sysr = getrom(R, Order = ra, Method="truncate");
        % Calculate the reduced model's frequency response
        H_r = freqresp(sysr, w);

        % Calculate the weighted error per frequency point, i.e.,
        % ||W^{-1}_j(jw) * E_j(jw) * V^{-1}_j(jw)||
        eEb = zeros(1, nw);
        for ii = find(vld.rss)
            eEb(ii) = norm(WTD.inv{ic, ii} * (H_r(:, :, ii) - H_j{ic}(:, :, ii)) *  VTD.inv{ic, ii});
        end

        % Stop if all the error agree with the bounds at all frequencies
        if all(eEb < 1)
            break;
        end
    end
    toc;

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

%% Subsystem reduction using greedy Loewner
for ic = 1:Nc

    greedyLoewner;

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

%%% Construction of the reduced interconnected system
sys_j(1) = sys_j_BT(1);
sys_j(2) = sys_j_AA(2);

sys_c.rss = lft(blkdiag(sys_j.rss), sys_I.s);
sys_c.rss.Name = 'RSS';

sys_c.stability =  isstable(sys_c.rss);

fprintf('\nStability of interconnection: %d\n', sys_c.stability);


%% (5) Evaluation of the Reduced Models and Output of Results (Fig. 6)%%
% Extract the names of the performed reduction methods
redNms = string(fieldnames(sys_c).');
redNms = setdiff(redNms,["m","s", "stability"]);

% Calculate the frequency responses of the reduced interconnected system
% and the associated error system
for nm = redNms
    Hr.(nm) = freqresp(sys_c.(nm), w); % \hat G_c(jw)
    Herr.(nm) = H_0 - Hr.(nm); % E_c(jw)
end

%%% Plot of the interconnected system ROM, FOM, and specs

% Add additional bounds to the frequency grid for patch plotting
fe = [1e-10, w, 1e10];

% Prespecified bounds for patch plotting:
% gbm represents the system +- the bounds
% ebm represents only the bounds itself

eb_cr = squeeze(eb_cr)';

gbm  = max([eps, pageNorm(H_0) - eb_cr, eps, eps, fliplr(pageNorm(H_0) + eb_cr), eps], eps);

% Normalised interconnected system error E_c
eEc = zeros(1, nw);
for ii = find(vld.rss)
    eEc(ii) = norm(V_cr(ii) * (H_0(:, :, ii) - Hr.rss(:, :, ii)));
end

fig = figure();
clf(fig);

tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top tile
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
ylabel(ax.Gc, 'Magnitude (dB)', 'interpreter', 'latex');
xlim(ax.Gc, [1e-1, 1e4]);
ylim(ax.Gc, [-90, 0]);
grid(ax.Gc, 'on');
set(ax.Gc, 'XScale', 'log');
set(ax.Gc, 'XTickLabel', []);

% Bottom tile
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


%%% Plot of the subsystem ROM, FOM, and specs
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

% Top tile: subsystem FOM, ROM, and reduction error specification
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


% Bottom tile: frequency-weighted normalised reduction error
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


%%% Poles
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

legend('Poles of full power grid model', ...
       'Poles of reduced power grid model', ...
        'interpreter', 'latex');

hold off;

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
