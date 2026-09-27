%%
t = 0:5e-5:1;

u = zeros(size(t));      % initial value 0
u(t > 0.1) = 0.03;          % rise to 0.4 after t = 1
u(t > 0.2) = 0.01;          % drop to 0.2 after t = 4

sys1 = lft(blkdiag(sys_j.rss), sys_I.s);
sys2 = lft(blkdiag(sys_j(1:2).s), sys_I.s);

tic; y1 = lsim(sys1, u, t); toc;
tic; y2 = lsim(sys2, u, t); toc;

% plot(t, u, t, y1, t, y2, 'LineWidth', 1.5);
% grid on;
% xlabel('Time [s]');
% ylabel('u(t)');
% title('Piecewise step input');

%%
% 1. Define time and input vectors
t = 0:5e-5:1;
u = zeros(size(t));      
u(t > 0.1) = 0.03;       
u(t > 0.2) = 0.01;       

% 2. Build the interconnected systems
% (Assuming sys_j and sys_I are already loaded in your workspace)
sys1 = lft(blkdiag(sys_j.rss), sys_I.s);
sys2 = lft(blkdiag(sys_j(1:2).s), sys_I.s);

% 3. Setup benchmarking parameters
num_runs = 10;
time_sys1 = zeros(num_runs, 1);
time_sys2 = zeros(num_runs, 1);

fprintf('Performing warm-up run to initialize JIT compiler...\n');
% Run once without timing so MATLAB can cache the functions
y1_warmup = lsim(sys1, u, t);
y2_warmup = lsim(sys2, u, t);

fprintf('Running %d iterations...\n', num_runs);

% 4. Execute the loop
for i = 1:num_runs
    % Time the Reduced Order Model (sys1)
    tic; 
    y1 = lsim(sys1, u, t); 
    time_sys1(i) = toc;
    
    % Time the Full Order Model (sys2)
    tic; 
    y2 = lsim(sys2, u, t); 
    time_sys2(i) = toc;
end

% 5. Calculate and display the averages
avg_time_sys1 = mean(time_sys1);
avg_time_sys2 = mean(time_sys2);

fprintf('\n=== Benchmarking Results (%d runs) ===\n', num_runs);
fprintf('Average time sys1 (ROM): %.4f seconds\n', avg_time_sys1);
fprintf('Average time sys2 (FOM): %.4f seconds\n', avg_time_sys2);

% Optional: Calculate the speedup factor
speedup = avg_time_sys2 / avg_time_sys1;
fprintf('Speedup factor:          %.2fx faster\n', speedup);

% %%
% 
% % Full-order and reduced-order systems
% % You should already have these:
% Gc = sys_c.s;
% Gr = sys_c.rss;
% 
% % E = minreal(Gc - Gr);
% E = Gc - Gr;
% 
% % Frequency response of Gc
% Gc_resp = squeeze(H_0);
% 
% % Define Wc^{-1}(jw)
% Wcinv = 0.25*abs(Gc_resp(:)) + 0.01;
% % Wcinv = abs(0.25*Gc_resp(:) + 0.01);
% 
% % Compute epsilon_r
% epsilon_r = max(Wcinv);
% 
% fprintf('epsilon_r = %.6e\n', epsilon_r);
% 
% % Frequency response of error system
% E_resp = squeeze(Herr.rss);
% sigE = abs(E_resp(:));
% 
% % Check sampled weighted bound
% weighted_ratio = sigE ./ Wcinv;
% gamma_weighted_sampled = max(weighted_ratio);
% 
% fprintf('Sampled max |E|/|Wc^{-1}| = %.6e\n', gamma_weighted_sampled);
% 
% % Plot frequency-domain error envelope
% figure;
% loglog(w, sigE, 'LineWidth', 1.5);
% hold on;
% loglog(w, Wcinv, '--', 'LineWidth', 1.5);
% grid on;
% xlabel('Frequency [rad/s]');
% ylabel('Magnitude');
% legend('|E_c(j\omega)|', 'W_c^{-1}(j\omega)', 'Location', 'best');
% title('Frequency-domain error bound');
% 
% % Time-domain input
% t = (0:0.00001:1)';
% 
% u = zeros(length(t), 1);
% u(t > 0.1 & t <= 0.2) = 0.03;
% % u(t > 0.2 & t <= 0.3) = 0.01;
% u(t > 0.2) = 0.01;
% 
% % Simulate full and reduced systems
% y  = lsim(Gc, u, t);
% yr = lsim(Gr, u, t);
% 
% e = y - yr;
% 
% % L2 norm calculation
% u_L2 = sqrt(trapz(t, u.^2));
% e_L2 = sqrt(trapz(t, e.^2));
% 
% ratio = e_L2 / u_L2;
% normalised_ratio = e_L2 / (epsilon_r*u_L2);
% 
% fprintf('||u||_2              = %.6e\n', u_L2);
% fprintf('||y - yr||_2         = %.6e\n', e_L2);
% fprintf('||y - yr||_2/||u||_2 = %.6e\n', ratio);
% fprintf('ratio / epsilon_r    = %.6e\n', normalised_ratio);
% 
% % Plot input
% figure;
% plot(t, u, 'LineWidth', 1.5);
% grid on;
% xlabel('Time [s]');
% ylabel('u(t)');
% title('Finite-energy input signal');
% 
% % Plot output error
% figure;
% plot(t, e, 'LineWidth', 1.5);
% grid on;
% xlabel('Time [s]');
% ylabel('e(t) = y(t) - y_r(t)');
% title('Time-domain output error');
% 
% % Plot cumulative L2 ratio
% cum_u_energy = cumtrapz(t, u.^2);
% cum_e_energy = cumtrapz(t, e.^2);
% 
% R = sqrt(cum_e_energy) ./ (epsilon_r*sqrt(cum_u_energy));
% R(cum_u_energy < 1e-12) = NaN;
% 
% fig = figure;
% semilogy(t, R, 'LineWidth', 1.6);
% hold on;
% yline(1, '--', 'LineWidth', 1.2);
% grid on;
% xlabel('Time ($s$)', 'Interpreter', 'latex', 'FontSize', 13.5);
% ylabel('Normalised ratio', 'Interpreter', 'latex', 'FontSize', 13.5);
% % Optional figure size for two-column-width publication figure
% set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);
% % title('Verification of ||y-y_r||_2 <= \epsilon_r ||u||_2');
% % legend('R(t)', 'Bound = 1', 'Interpreter', 'latex');

%%
ic = 2;

Gpos = freqresp(linSys2, w); % Gpos = freqresp(sys_j_AA(ic).s, w);
GfitPos = freqresp(sys_j_AA(ic).rss, w);

weightedErr = zeros(numel(w),1);
unweightedErr = zeros(numel(w),1);

% Weights
WoErrBound = zeros(size(sys_j(ic).s.D, 1), size(sys_j(ic).s.D, 1), numel(w));
WiErrBound = zeros(size(sys_j(ic).s.D, 2), size(sys_j(ic).s.D, 2), numel(w));

for k = 1:numel(w)
    WoErrBound(:,:,k) = WTD.inv{ic, k};
    WiErrBound(:,:,k) = VTD.inv{ic, k};
end

for k = 1:numel(w)
    Rk = GfitPos(:,:,k) - Gpos(:,:,k);
    weightedErr(k) = norm(WoErrBound(:,:,k)*Rk*WiErrBound(:,:,k),2);
    unweightedErr(k) = norm(Rk,2);
end

% fprintf('\nFinal descriptor order: %d\n',model.order);
fprintf('Max positive-frequency weighted error: %.4e\n',max(weightedErr));
fprintf('Grid specification satisfied: %d\n',max(weightedErr) <= 1);
% fprintf('Imaginary leakage after real transform: %.4e\n',info.realificationImaginaryLeakage);
% fprintf('Max imag part of returned matrices: %.4e\n',max([ ...
%     norm(imag(model.E),'fro'), norm(imag(model.A),'fro'), ...
%     norm(imag(model.B),'fro'), norm(imag(model.C),'fro'), norm(imag(model.D),'fro') ]));

figure;
semilogx(w,weightedErr,'LineWidth',1.2); grid on;
yline(1,'--');
xlabel('Frequency [rad/s]');
ylabel('||Wo(jw)(Gr(jw)-G(jw))Wi(jw)||_2');
title('Point-wise weighted error');

% model.system = dss(model.A, model.B, model.C, model.D, model.E);
% isstable(model.system)

figure;
semilogx(w,unweightedErr,'LineWidth',1.2); grid on;
xlabel('Frequency [rad/s]');
ylabel('||Gr(jw)-G(jw)||_2');
title('Unweighted error');


%%
Gpos = freqresp(sys_c.s, w);
GfitPos = freqresp(sys_c.rss, w);

weightedErr = zeros(numel(w),1);
unweightedErr = zeros(numel(w),1);

% Weights
WoErrBound = V_cr;
WiErrBound = W_cr;

for k = 1:numel(w)
    Rk = GfitPos(:,:,k) - Gpos(:,:,k);
    weightedErr(k) = norm(WoErrBound(:,:,k)*Rk*WiErrBound(:,:,k),2);
    unweightedErr(k) = norm(Rk,2);
end

% fprintf('\nFinal descriptor order: %d\n',model.order);
fprintf('Max positive-frequency weighted error: %.4e\n',max(weightedErr));
fprintf('Grid specification satisfied: %d\n',max(weightedErr) <= 1);
% fprintf('Imaginary leakage after real transform: %.4e\n',info.realificationImaginaryLeakage);
% fprintf('Max imag part of returned matrices: %.4e\n',max([ ...
%     norm(imag(model.E),'fro'), norm(imag(model.A),'fro'), ...
%     norm(imag(model.B),'fro'), norm(imag(model.C),'fro'), norm(imag(model.D),'fro') ]));

figure;
semilogx(w,weightedErr,'LineWidth',1.2); grid on;
yline(1,'--');
xlabel('Frequency [rad/s]');
ylabel('||Wo(jw)(Gr(jw)-G(jw))Wi(jw)||_2');
title('Point-wise weighted error');

% model.system = dss(model.A, model.B, model.C, model.D, model.E);
% isstable(model.system)

figure;
semilogx(w,unweightedErr,'LineWidth',1.2); grid on;
xlabel('Frequency [rad/s]');
ylabel('||Gr(jw)-G(jw)||_2');
title('Unweighted error');


%% Poles
% Example systems
ic = 1;

sys1 = sys_j_BT(ic).s;
sys2 = sys_j_BT(ic).rss;

% Extract poles
p1 = pole(sys1);
p2 = pole(sys2);

%%%
% Frequency band: ±50 Hz
fBand = 50;                 % Hz
wBand = 2*pi*fBand;         % rad/s

% Select poles whose imaginary part is within ±50 Hz
idxBand1 = abs(imag(p1)) <= wBand;
idxBand2 = abs(imag(p2)) <= wBand;

p1_band = p1(idxBand1);
p2_band = p2(idxBand2);

% Dominant pole = rightmost pole in this frequency band
[~, idxDom1] = max(real(p1_band));
[~, idxDom2] = max(real(p2_band));

p1_dom = p1_band(idxDom1);
p2_dom = p2_band(idxDom2);

% Export useful quantities
domPoleTable = table( ...
    ["FOM"; "ROM"], ...
    [real(p1_dom); real(p2_dom)], ...
    [imag(p1_dom); imag(p2_dom)], ...
    [imag(p1_dom)/(2*pi); imag(p2_dom)/(2*pi)], ...
    [abs(p1_dom); abs(p2_dom)], ...
    [-real(p1_dom)/abs(p1_dom); -real(p2_dom)/abs(p2_dom)], ...
    'VariableNames', {'Model','RealPart','ImagPart_rad_s','ImagPart_Hz','AbsPole','DampingRatio'} ...
);

disp(domPoleTable);
%%%

% Find pole closest to imaginary axis for each system
[~, idx1] = min(abs(real(p1)));
[~, idx2] = min(abs(real(p2)));

p1_close = p1(idx1);
p2_close = p2(idx2);

figure;

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
% ylim([-3000, 3000]);

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
xlim([-1200, 200]);
ylim([-500, 500]);

% title('Inset Zoom');
hold off;

%%
%% Compact comparison figure:
% Two upper magnitude plots and one lower normalized-error plot

ic = 2;

oCn = 2;
iCn = 1;

e_jr = zeros(p_j(ic), 1);
e_jr(oCn) = 1;

e_ir = zeros(m_j(ic), 1);
e_ir(iCn) = 1;

eb_VW_ji = zeros(1, nw);

for jj = 1:nw
    eb_VW_ji(jj) = norm(WTD.inv{ic, jj} \ e_jr) * ...
                    norm(VTD.inv{ic, jj} \ e_ir);
end

gbm  = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_ji, eps, eps, fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_ji), eps], eps); % Original: pageNorm(H_j{ic}) - eb_VW(ic, :)
ebm  = [eps, eb_VW_ji, eps * ones(1, numel(fe) + 1)]; % Original: eb_VW(ic, :)


% -------------------------------------------------------------------------
% Create compact tiled figure
% -------------------------------------------------------------------------
fig = figure();
clf(fig);

tl = tiledlayout(fig, 2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');


% ========================= Upper-left tile ===============================
ax.Gj = nexttile(tl, 1);
hold(ax.Gj, 'on');
set(ax.Gj, 'Layer', 'top');

patch(ax.Gj, [fe fliplr(fe)], 20 * log10(gbm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);

semilogx(ax.Gj, w, 20 * log10(pageNorm(H_j{ic}(oCn, iCn, :))), ...
    'b-', 'LineWidth', 2.0);

semilogx(ax.Gj, w, 20 * log10(pageNorm(Hr_j(ic).rss(oCn, iCn, :))), ...
    'r--', 'LineWidth', 2.0);

ylabel(ax.Gj, 'Magnitude (dB)', 'Interpreter', 'latex');
% xlabel(ax.Gj, 'Frequency (rad/s)', 'interpreter', 'latex');

legend(ax.Gj, ...
    'Acc. spec.', 'FOM', 'ROM', ...
    'Interpreter', 'latex', ...
    'Location', 'best');

grid(ax.Gj, 'on');
set(ax.Gj, 'XScale', 'log');
xlim(ax.Gj, [w(1), w(end)]);
ylim(ax.Gj, [-70, 0]);

set(ax.Gj, 'XTickLabel', []);



% ========================= Upper-right tile ==============================
oCn = 4;
iCn = 5;

e_jr = zeros(p_j(ic), 1);
e_jr(oCn) = 1;

e_ir = zeros(m_j(ic), 1);
e_ir(iCn) = 1;

eb_VW_ji = zeros(1, nw);

for jj = 1:nw
    eb_VW_ji(jj) = norm(WTD.inv{ic, jj} \ e_jr) * ...
                    norm(VTD.inv{ic, jj} \ e_ir);
end

gbm  = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_ji, eps, eps, fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_ji), eps], eps); % Original: pageNorm(H_j{ic}) - eb_VW(ic, :)
ebm  = [eps, eb_VW_ji, eps * ones(1, numel(fe) + 1)]; % Original: eb_VW(ic, :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax.Gj = nexttile(tl, 2);
hold(ax.Gj, 'on');
set(ax.Gj, 'Layer', 'top');

patch(ax.Gj, [fe fliplr(fe)], 20 * log10(gbm), [.8 .8 1], 'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);

semilogx(ax.Gj, w, 20 * log10(pageNorm(H_j{ic}(oCn, iCn, :))), ...
    'b-', 'LineWidth', 2.0);

semilogx(ax.Gj, w, 20 * log10(pageNorm(Hr_j(ic).rss(oCn, iCn, :))), ...
    'r--', 'LineWidth', 2.0);

% ylabel(ax.Gj, 'Magnitude (dB)', 'Interpreter', 'latex');
% xlabel(ax.Gj, 'Frequency (rad/s)', 'interpreter', 'latex');

legend(ax.Gj, ...
    'Acc. spec.', 'FOM', 'ROM', ...
    'Interpreter', 'latex', ...
    'Location', 'best');

grid(ax.Gj, 'on');
set(ax.Gj, 'XScale', 'log');
xlim(ax.Gj, [w(1), w(end)]);
ylim(ax.Gj, [-70, 0]);

set(ax.Gj, 'XTickLabel', []);


% ========================= Bottom spanning tile ==========================
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax.nEj = nexttile(tl, [1, 2]);
set(ax.nEj, 'Layer', 'top');
hold(ax.nEj, 'on');

patch(ax.nEj, ...
    [fe fliplr(fe)], ...
    [ones(size(fe)), fliplr(eps * ones(size(fe)))], ...
    [.8 .8 1], ...
    'LineWidth', 1.2, ...
    'LineStyle', ':', ...
    'EdgeColor', [0 0 0]);

semilogx(ax.nEj, w, eEj, 'r--', 'linewidth', 2.0);

legend(ax.nEj, ...
    'Red. err. spec.', ...
    'Weigh. red. err.', ...
    'interpreter', 'latex');

xlabel(ax.nEj, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.nEj, 'Weighted error', ...
    'interpreter', 'latex');

xlim(ax.nEj, [w(1), w(end)]);
ylim(ax.nEj, [1e-2*0.6, 1e0 * 1.1]);

grid(ax.nEj, 'on');
set(ax.nEj, 'XScale', 'log');
set(ax.nEj, 'YScale', 'log');

% Link frequency axes
% linkaxes([ax.Gj, ax.nEj], 'x');


% -------------------------------------------------------------------------
% Shared formatting
% -------------------------------------------------------------------------
% set([ax.Gj, ax.Gj, ax.Err], ...
%     'FontName', 'Times New Roman', ...
%     'FontSize', 10, ...
%     'TickLabelInterpreter', 'latex', ...
%     'Box', 'on');

% linkaxes([ax.Gj, ax.Gj, ax.Err], 'x');

% Optional figure size for two-column-width publication figure
set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

% Optional vector export
% exportgraphics(fig, 'compact_MOR_comparison.pdf', 'ContentType', 'vector');

%%
% -------------------------------------------------------------------------
% Preparation: Format frequency grids
% -------------------------------------------------------------------------
w = w(:).';
fe = [1e-10, w, 1e10]; % Extended bounds for patch plotting

Hr_j(1).rss = freqresp(sys_j_BT(1).rss, w);
Hr_j(2).rss = freqresp(sys_j_AA(2).rss, w);

% -------------------------------------------------------------------------
% Create compact 3x1 tiled figure
% -------------------------------------------------------------------------
fig = figure();
clf(fig);
tl = tiledlayout(fig, 3, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% ========================= TOP TILE: ic = 1 (4 to 4) =====================
ic = 1;
oCn = 7;
iCn = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WoSys1 = freqresp(Gy{1}, w);
WiSys1 = freqresp(Gu{1}, w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate bounds for Tile 1
e_jr = zeros(p_j(ic), 1); e_jr(oCn) = 1;
e_ir = zeros(m_j(ic), 1); e_ir(iCn) = 1;
%%%
eb_VW_ji = zeros(1, nw);
for jj = 1:nw
    eb_VW_ji(jj) = norm(WTD.inv{ic, jj} \ e_jr) * norm(VTD.inv{ic, jj} \ e_ir);
    % eb_VW_ji(jj) = norm(WoSys1(:, :, jj) \ e_jr) * norm(WiSys1(:, :, jj) \ e_ir);
end
%%%
%%%
% Allocate two different bounds
eb_VW_original = zeros(1, nw);   % from computed pointwise WTD/VTD
eb_VW_fitted   = zeros(1, nw);   % from fitted weighting systems WoSys1/WiSys1

for jj = 1:nw
    % Bound from computed frequency-dependent weights
    eb_VW_original(jj) = norm(WTD.inv{ic,jj} \ e_jr) * ...
                         norm(VTD.inv{ic,jj} \ e_ir);

    % Bound from fitted weighting systems
    eb_VW_fitted(jj) = norm(WoSys1(:,:,jj) \ e_jr) * ...
                       norm(WiSys1(:,:,jj) \ e_ir);
end
%%%

gbm = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_ji, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_ji), eps], eps);

gbm_original = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_original, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_original), eps], eps);

gbm_fitted = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_fitted, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_fitted), eps], eps);

ax.top = nexttile(tl, 1);
hold(ax.top, 'on');
set(ax.top, 'Layer', 'top');

% Magnitude bound patch
% patch(ax.top, [fe fliplr(fe)], 20 * log10(gbm), [.8 .8 1], ...
%     'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);

p_original = patch(ax.top, [fe fliplr(fe)], 20 * log10(gbm_original), [.8 .8 1], ...
    'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);

p_fitted = patch(ax.top, [fe fliplr(fe)], 20 * log10(gbm_fitted), [1 .85 .75], ...
    'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0.75 0.25 0]);


% Magnitude plots
semilogx(ax.top, w, 20 * log10(pageNorm(H_j{ic}(oCn, iCn, :))), 'b-', 'LineWidth', 2.0);
semilogx(ax.top, w, 20 * log10(pageNorm(Hr_j(ic).rss(oCn, iCn, :))), 'r--', 'LineWidth', 2.0);

ylabel(ax.top, 'Magnitude (dB)', 'Interpreter', 'latex');
legend(ax.top, sprintf('Acc. spec. %d', ic), sprintf('Fit. Acc. spec. %d', ic), 'FOM', 'ROM', ...
    'Interpreter', 'latex', 'Location', 'best');
grid(ax.top, 'on');
set(ax.top, 'XScale', 'log');
xlim(ax.top, [1e-1, 1e4]);
ylim(ax.top, [-20, 30]);
set(ax.top, 'XTickLabel', []); % Hide X-axis labels for clean stacking

% ========================= MIDDLE TILE: ic = 2 (1 to 2) ==================
ic = 2;
oCn = 2;
iCn = 1;

% Calculate bounds for Tile 2
e_jr = zeros(p_j(ic), 1); e_jr(oCn) = 1;
e_ir = zeros(m_j(ic), 1); e_ir(iCn) = 1;
eb_VW_ji = zeros(1, nw);
for jj = 1:nw
    eb_VW_ji(jj) = norm(WTD.inv{ic, jj} \ e_jr) * norm(VTD.inv{ic, jj} \ e_ir);
end
gbm = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_ji, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_ji), eps], eps);

ax.mid = nexttile(tl, 2);
hold(ax.mid, 'on');
set(ax.mid, 'Layer', 'top');

% Magnitude bound patch
patch(ax.mid, [fe fliplr(fe)], 20 * log10(gbm), [.8 .8 1], ...
    'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);

% Magnitude plots
semilogx(ax.mid, w, 20 * log10(pageNorm(H_j{ic}(oCn, iCn, :))), 'b-', 'LineWidth', 2.0);
semilogx(ax.mid, w, 20 * log10(pageNorm(Hr_j(ic).rss(oCn, iCn, :))), 'r--', 'LineWidth', 2.0);

ylabel(ax.mid, 'Magnitude (dB)', 'Interpreter', 'latex');
legend(ax.mid, sprintf('Acc. spec. %d', ic), 'FOM', 'ROM', ...
    'Interpreter', 'latex', 'Location', 'best');
grid(ax.mid, 'on');
set(ax.mid, 'XScale', 'log');
xlim(ax.mid, [1e-1, 1e4]);
ylim(ax.mid, [-50, 10]);
set(ax.mid, 'XTickLabel', []); % Hide X-axis labels for clean stacking
% xlabel(ax.mid, 'Frequency (rad/s)', 'interpreter', 'latex');

% ========================= BOTTOM TILE: Weighted Error ===================
% Calculate normalized weighted error E_j for both ic = 1 and ic = 2
eEj1 = NaN(1, nw);
eEj2 = NaN(1, nw);

for ii = find(vld.rss)
    % Error for ic = 1
    eEj1(ii) = norm(WTD.inv{1, ii} * ...
        (Hr_j(1).rss(:, :, ii) - H_j{1}(:, :, ii)) * VTD.inv{1, ii});
    % Error for ic = 2
    eEj2(ii) = norm(WTD.inv{2, ii} * ...
        (Hr_j(2).rss(:, :, ii) - H_j{2}(:, :, ii)) * VTD.inv{2, ii});
end

ax.bot = nexttile(tl, 3);
set(ax.bot, 'Layer', 'top');
hold(ax.bot, 'on');

% Plot the constant error specification line at y = 1 (No patch)
loglog(ax.bot, [w(1), w(end)], [1, 1], 'k:', 'LineWidth', 1.5);

% Plot the errors for both subsystems
loglog(ax.bot, w, eEj1, 'r--', 'linewidth', 2.0);
loglog(ax.bot, w, eEj2, 'b--', 'linewidth', 2.0);

% Labels and Formatting
legend(ax.bot, 'Acc. spec.', 'Weigh. err. 1', 'Weigh. err. 2)', ...
    'interpreter', 'latex', 'Location', 'best');
xlabel(ax.bot, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.bot, 'Weighted error', 'interpreter', 'latex');

xlim(ax.bot, [1e-1, 1e4]);
ylim(ax.bot, [1e-2, 1e0 * 1.3]);
grid(ax.bot, 'on');
set(ax.bot, 'XScale', 'log');
set(ax.bot, 'YScale', 'log');

% -------------------------------------------------------------------------
% Figure Sizing for IEEE Publication
% -------------------------------------------------------------------------
% Adjusted height from 4.6 to 6.5 to properly fit 3 vertically stacked tiles
set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 6.5]);

%%

% 3. Setup benchmarking parameters
num_runs = 10;
time_sys1 = zeros(num_runs, 1);

fprintf('Performing warm-up run to initialize JIT compiler...\n');
% Run once without timing so MATLAB can cache the functions
sim("simulink_files\subsystem_118_bus");

fprintf('Running %d iterations...\n', num_runs);

% 4. Execute the loop
for i = 1:num_runs
    % Time the Reduced Order Model (sys1)
    tic; 
    sim("simulink_files\subsystem_118_bus");
    time_sys1(i) = toc;
end

% 5. Calculate and display the averages
avg_time_sys1 = mean(time_sys1);

fprintf('\n=== Benchmarking Results (%d runs) ===\n', num_runs);
fprintf('Average time sys1: %.4f seconds\n', avg_time_sys1);

%%
%% CHECK UPDATED INTERCONNECTED ERROR: || V_c (Ghat_c^+ - G_c^+) W_c || <= 1

tol = 1e-8;

% ------------------------------------------------------------
% 1. Build updated full-order interconnected model G_c^+
% ------------------------------------------------------------
sys_j_plus = sys_j;

% Replace subsystem 2 by the updated full-order subsystem
% Use your actual variable name here:
sys_j_plus(2).s = linSys2_update;          % G_2^+

sys_c_plus.s = lft(blkdiag(sys_j_plus.s), sys_I.s);
sys_c_plus.s.Name = 'Gc_plus';

% ------------------------------------------------------------
% 2. Build updated reduced interconnected model Ghat_c^+
% ------------------------------------------------------------
sys_j_plus_r = sys_j_plus;

% Keep old ROMs for subsystems 1 and 3
sys_j_plus_r(1).rss = sys_j_AA(1).rss;

% Replace subsystem 2 ROM by the new reduced model
% Use your actual variable name here:
sys_j_plus_r(2).rss = model.sys;   % \hat G_2^+

sys_c_plus.rss = lft(blkdiag(sys_j_plus_r.rss), sys_I.s);
sys_c_plus.rss.Name = 'Gchat_plus';

% ------------------------------------------------------------
% 3. Frequency responses
% ------------------------------------------------------------
Hc_plus    = freqresp(sys_c_plus.s,   w);   % G_c^+(jw)
Hchat_plus = freqresp(sys_c_plus.rss, w);   % \hat G_c^+(jw)

Ec_plus = Hchat_plus - Hc_plus;             % E_c^+(jw)

% ------------------------------------------------------------
% 4. Choose the global error specification
% ------------------------------------------------------------
% Option A: reuse the old global specification
% Then keep V_cr and W_cr from your old code.
V_cr_plus = V_cr;
W_cr_plus = W_cr;    % usually identity in your code


% % Option B: recompute the bound around the updated full model G_c^+
% % This matches your existing choice eb_c = 0.25 * sys_c.s + 1e-2.
% eb_c_plus_resp = freqresp(0.25 * sys_c_plus.s, w);
% eb_cr_plus = pagenorm(eb_c_plus_resp, 2) + 1e-2;
% 
% % Normalized scalar output weighting:
% V_cr_plus = 1 ./ squeeze(eb_cr_plus);
% W_cr_plus = W_cr;    % usually identity in your code

% ------------------------------------------------------------
% 5. Check normalized error <= 1
% ------------------------------------------------------------
eEc_plus = NaN(1, nw);
rawEc_plus = NaN(1, nw);

for ii = 1:nw
    Eii = Ec_plus(:, :, ii);

    rawEc_plus(ii) = norm(Eii, 2);

    % Your current code uses scalar V_cr(ii).
    % This corresponds to V_c(jw) = scalar * I.
    % Vci = V_cr_plus(ii) * eye(size(Eii, 1));
    Vci = V_cr_plus(:, :, ii);

    % Your W_c is identity, but keep this general.
    if ndims(W_cr_plus) == 3
        Wci = W_cr_plus(:, :, ii);
    else
        Wci = eye(size(Eii, 2));
    end

    eEc_plus(ii) = norm(Vci * Eii * Wci, 2);
end

% ------------------------------------------------------------
% 6. Pass/fail result
% ------------------------------------------------------------
[peakErr, idxPeak] = max(eEc_plus);

pass_updated_global = all(eEc_plus <= 1 + tol);

fprintf('\nUpdated interconnected error check:\n');
fprintf('  pass = %d\n', pass_updated_global);
fprintf('  max ||Vc*(Ghat_c^+ - G_c^+)*Wc|| = %.6e\n', peakErr);
fprintf('  peak frequency = %.6e rad/s\n', w(idxPeak));

figure;
loglog(w, eEc_plus, 'r--', 'LineWidth', 2.0); hold on;
loglog([w(1), w(end)], [1, 1], 'k:', 'LineWidth', 2.0);
grid on;

xlabel('Frequency (rad/s)', 'interpreter', 'latex');
ylabel('$\bar\sigma(V_c(\hat G_c^+ - G_c^+)W_c)$', 'interpreter', 'latex');
legend('Updated weighted error', 'Specification', 'interpreter', 'latex');

title('Check of updated interconnected reduction error');


%%
% -------------------------------------------------------------------------
% Preparation: Format frequency grids
% -------------------------------------------------------------------------
w = w(:).';
fe = [1e-10, w, 1e10]; % Extended bounds for patch plotting

Hr_j(1).rss = freqresp(sys_j_BT(1).rss, w);
Hr_j(2).rss = freqresp(sys_j_AA(2).rss, w);

% -------------------------------------------------------------------------
% Create compact 3x1 tiled figure
% -------------------------------------------------------------------------
fig = figure();
clf(fig);
tl = tiledlayout(fig, 3, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% ========================= TOP TILE: ic = 1 (4 to 4) =====================
ic = 1;
oCn = 7;
iCn = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WoSys1 = freqresp(Gy{1}, w);
WiSys1 = freqresp(Gu{1}, w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate bounds for Tile 1
e_jr = zeros(p_j(ic), 1); e_jr(oCn) = 1;
e_ir = zeros(m_j(ic), 1); e_ir(iCn) = 1;
%%%
eb_VW_ji = zeros(1, nw);
for jj = 1:nw
    eb_VW_ji(jj) = norm(WTD.inv{ic, jj} \ e_jr) * norm(VTD.inv{ic, jj} \ e_ir);
    % eb_VW_ji(jj) = norm(WoSys1(:, :, jj) \ e_jr) * norm(WiSys1(:, :, jj) \ e_ir);
end
%%%
%%%
% Allocate two different bounds
eb_VW_original = zeros(1, nw);   % from computed pointwise WTD/VTD
eb_VW_fitted   = zeros(1, nw);   % from fitted weighting systems WoSys1/WiSys1

for jj = 1:nw
    % Bound from computed frequency-dependent weights
    eb_VW_original(jj) = norm(WTD.inv{ic,jj} \ e_jr) * ...
                         norm(VTD.inv{ic,jj} \ e_ir);

    % Bound from fitted weighting systems
    eb_VW_fitted(jj) = norm(WoSys1(:,:,jj) \ e_jr) * ...
                       norm(WiSys1(:,:,jj) \ e_ir);
end
%%%

gbm = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_ji, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_ji), eps], eps);

gbm_original = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_original, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_original), eps], eps);

gbm_fitted = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_fitted, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_fitted), eps], eps);

ax.top = nexttile(tl, 1);
hold(ax.top, 'on');
set(ax.top, 'Layer', 'top');

% Magnitude bound patch
% patch(ax.top, [fe fliplr(fe)], 20 * log10(gbm), [.8 .8 1], ...
%     'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);

p_original = patch(ax.top, [fe fliplr(fe)], 20 * log10(gbm_original), [.8 .8 1], ...
    'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);

p_fitted = patch(ax.top, [fe fliplr(fe)], 20 * log10(gbm_fitted), [1 .85 .75], ...
    'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0.75 0.25 0]);


% Magnitude plots
semilogx(ax.top, w, 20 * log10(pageNorm(H_j{ic}(oCn, iCn, :))), 'b-', 'LineWidth', 2.0);
semilogx(ax.top, w, 20 * log10(pageNorm(Hr_j(ic).rss(oCn, iCn, :))), 'r--', 'LineWidth', 2.0);

ylabel(ax.top, 'Magnitude (dB)', 'Interpreter', 'latex');
legend(ax.top, sprintf('Acc. spec. %d', ic), sprintf('Fit. Acc. spec. %d', ic), 'FOM', 'ROM', ...
    'Interpreter', 'latex', 'Location', 'best');
grid(ax.top, 'on');
set(ax.top, 'XScale', 'log');
xlim(ax.top, [1e-1, 1e4]);
ylim(ax.top, [-20, 30]);
set(ax.top, 'XTickLabel', []); % Hide X-axis labels for clean stacking

% ========================= MIDDLE TILE: ic = 2 (1 to 2) ==================
ic = 2;
oCn = 2;
iCn = 1;

% Calculate bounds for Tile 2
e_jr = zeros(p_j(ic), 1); e_jr(oCn) = 1;
e_ir = zeros(m_j(ic), 1); e_ir(iCn) = 1;
eb_VW_ji = zeros(1, nw);
for jj = 1:nw
    eb_VW_ji(jj) = norm(WTD.inv{ic, jj} \ e_jr) * norm(VTD.inv{ic, jj} \ e_ir);
end
gbm = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_ji, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_ji), eps], eps);

ax.mid = nexttile(tl, 2);
hold(ax.mid, 'on');
set(ax.mid, 'Layer', 'top');

% Magnitude bound patch
patch(ax.mid, [fe fliplr(fe)], 20 * log10(gbm), [.8 .8 1], ...
    'LineWidth', 1.2, 'LineStyle', ':', 'EdgeColor', [0 0 0]);

% Magnitude plots
semilogx(ax.mid, w, 20 * log10(pageNorm(H_j{ic}(oCn, iCn, :))), 'b-', 'LineWidth', 2.0);
semilogx(ax.mid, w, 20 * log10(pageNorm(Hr_j(ic).rss(oCn, iCn, :))), 'r--', 'LineWidth', 2.0);

ylabel(ax.mid, 'Magnitude (dB)', 'Interpreter', 'latex');
legend(ax.mid, sprintf('Acc. spec. %d', ic), 'FOM', 'ROM', ...
    'Interpreter', 'latex', 'Location', 'best');
grid(ax.mid, 'on');
set(ax.mid, 'XScale', 'log');
xlim(ax.mid, [1e-1, 1e4]);
ylim(ax.mid, [-50, 10]);
set(ax.mid, 'XTickLabel', []); % Hide X-axis labels for clean stacking
xlabel(ax.mid, 'Frequency (rad/s)', 'interpreter', 'latex');

% ========================= BOTTOM TILE: Weighted Error ===================
% Calculate normalized weighted error E_j for both ic = 1 and ic = 2
eEj1 = NaN(1, nw);
eEj2 = NaN(1, nw);

for ii = find(vld.rss)
    % Error for ic = 1
    eEj1(ii) = norm(WTD.inv{1, ii} * ...
        (Hr_j(1).rss(:, :, ii) - H_j{1}(:, :, ii)) * VTD.inv{1, ii});
    % Error for ic = 2
    eEj2(ii) = norm(WTD.inv{2, ii} * ...
        (Hr_j(2).rss(:, :, ii) - H_j{2}(:, :, ii)) * VTD.inv{2, ii});
end

ax.bot = nexttile(tl, 3);
set(ax.bot, 'Layer', 'top');
hold(ax.bot, 'on');

% Plot the constant error specification line at y = 1 (No patch)
loglog(ax.bot, [w(1), w(end)], [1, 1], 'k:', 'LineWidth', 1.5);

% Plot the errors for both subsystems
loglog(ax.bot, w, eEj1, 'r--', 'linewidth', 2.0);
loglog(ax.bot, w, eEj2, 'b--', 'linewidth', 2.0);

% Labels and Formatting
legend(ax.bot, 'Acc. spec.', 'Weigh. err. 1', 'Weigh. err. 2)', ...
    'interpreter', 'latex', 'Location', 'best');
xlabel(ax.bot, 'Frequency (rad/s)', 'interpreter', 'latex');
ylabel(ax.bot, 'Weighted error', 'interpreter', 'latex');

xlim(ax.bot, [1e-1, 1e4]);
ylim(ax.bot, [1e-2, 1e0 * 1.3]);
grid(ax.bot, 'on');
set(ax.bot, 'XScale', 'log');
set(ax.bot, 'YScale', 'log');

% -------------------------------------------------------------------------
% Freeze subplot geometry, remove bottom tile, preserve/recreate legends
% -------------------------------------------------------------------------

drawnow;

% Save exact positions of the original top and middle axes
pos_top = ax.top.Position;
pos_mid = ax.mid.Position;

% Delete the bottom subplot
delete(ax.bot);

% Re-parent the two remaining axes directly to the figure
ax.top.Parent = fig;
ax.mid.Parent = fig;

% Delete the tiled layout
delete(tl);

% Restore units and exact original positions
ax.top.Units = 'normalized';
ax.mid.Units = 'normalized';

ax.top.Position = pos_top;
ax.mid.Position = pos_mid;

% -------------------------------------------------------------------------
% Restore axis formatting
% -------------------------------------------------------------------------

% Top subplot: keep x tick labels hidden
set(ax.top, 'XTickLabel', []);

% Middle subplot: show x tick labels
set(ax.mid, 'XTickLabelMode', 'auto');

xlabel(ax.mid, ...
    'Frequency (rad/s)', ...
    'Interpreter', 'latex');

% -------------------------------------------------------------------------
% Recreate legends
% -------------------------------------------------------------------------

legend(ax.top, ...
    'Acc. spec. 1', ...
    'Fit. Acc. spec. 1', ...
    'FOM 1', ...
    'ROM 1', ...
    'Interpreter', 'latex', ...
    'Location', 'best');

legend(ax.mid, ...
    'Acc. spec. 2', ...
    'FOM 2', ...
    'ROM 2', ...
    'Interpreter', 'latex', ...
    'Location', 'best');

% -------------------------------------------------------------------------
% Keep original figure size
% -------------------------------------------------------------------------

set(fig, ...
    'Units', 'inches', ...
    'Position', [1, 1, 7.16, 6.5]);