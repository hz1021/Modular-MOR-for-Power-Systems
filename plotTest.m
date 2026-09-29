%% Table II: Computational comparison: Linear case
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


%% Fig. 8: Poles
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


%% Fig. 7: External-area subsystem frequency response: FOM vs. ROM
% Preparation: Format frequency grids
w = w(:).';
fe = [1e-10, w, 1e10]; % Extended bounds for patch plotting

Hr_j(1).rss = freqresp(sys_j_BT(1).rss, w);
Hr_j(2).rss = freqresp(sys_j_AA(2).rss, w);

% Create compact 3x1 tiled figure
fig = figure();
clf(fig);
tl = tiledlayout(fig, 3, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% TOP TILE: ic = 1 (4 to 4)
ic = 1;
oCn = 7;
iCn = 7;

WoSys1 = freqresp(Gy{1}, w);
WiSys1 = freqresp(Gu{1}, w);

% Calculate bounds for Tile 1
e_jr = zeros(p_j(ic), 1); e_jr(oCn) = 1;
e_ir = zeros(m_j(ic), 1); e_ir(iCn) = 1;

eb_VW_ji = zeros(1, nw);
for jj = 1:nw
    eb_VW_ji(jj) = norm(WTD.inv{ic, jj} \ e_jr) * norm(VTD.inv{ic, jj} \ e_ir);
end

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


gbm = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_ji, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_ji), eps], eps);

gbm_original = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_original, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_original), eps], eps);

gbm_fitted = max([eps, pageNorm(H_j{ic}(oCn, iCn, :)) - eb_VW_fitted, eps, eps, ...
           fliplr(pageNorm(H_j{ic}(oCn, iCn, :)) + eb_VW_fitted), eps], eps);

ax.top = nexttile(tl, 1);
hold(ax.top, 'on');
set(ax.top, 'Layer', 'top');


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

% MIDDLE TILE: ic = 2 (1 to 2)
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


% BOTTOM TILE: Weighted Error
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


% Adjusted height from 4.6 to 6.5 to properly fit 3 vertically stacked tiles
set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 6.5]);

%% Table II: Computational comparison: Nonlinear case

main_sub_118_bus;

num_runs = 10;
time_sys1 = zeros(num_runs, 1);

fprintf('Performing warm-up run to initialize JIT compiler...\n');
% Run once without timing so MATLAB can cache the functions
sim("subsystem_118_bus");

fprintf('Running %d iterations...\n', num_runs);

for i = 1:num_runs
    % Time the Reduced Order Model (sys1)
    tic; 
    sim("subsystem_118_bus");
    time_sys1(i) = toc;
end

% 5. Calculate and display the averages
avg_time_sys1 = mean(time_sys1);

fprintf('\n=== Benchmarking Results (%d runs) ===\n', num_runs);
fprintf('Average time sys1: %.4f seconds\n', avg_time_sys1);