%%% Load equilibrium data
load nominalU.mat;
load nominalY.mat;
load nominalInter.mat

Tstop = 1;
%%% Sim the .slx file
P_input = 0.03;

%% Fig. 9: Active power and error
out = sim("subsystem_118_bus");

%%% Extract time-domain data from .slx simulation
tout = squeeze(out.tout);

P_ref = squeeze(out.input.Data);

P_ROM = squeeze(out.P_ROM.Data);
P_linFOM = squeeze(out.P_linFOM.Data);
P_nonlinFOM = squeeze(out.P_nonlinFOM.Data);
P_ROM1 = squeeze(out.P_ROM1.Data);
P_linFOM1 = squeeze(out.P_linFOM1.Data);

sub3_5_1_ROM = squeeze(out.sub3_5_1_ROM.Data);
sub3_5_1_linFOM = squeeze(out.sub3_5_1_linFOM.Data);
sub3_5_1 = squeeze(out.sub3_5_1.Data);
sub3_5_1_ROM1 = squeeze(out.sub3_5_1_ROM1.Data);
sub3_5_1_linFOM1 = squeeze(out.sub3_5_1_linFOM1.Data);

sub1_6_1_ROM = squeeze(out.sub1_6_1_ROM.Data);
sub1_6_1_linFOM = squeeze(out.sub1_6_1_linFOM.Data);
sub1_6_1 = squeeze(out.sub1_6_1.Data);
sub1_6_1_ROM1 = squeeze(out.sub1_6_1_ROM1.Data);
sub1_6_1_linFOM1 = squeeze(out.sub1_6_1_linFOM1.Data);

sub2_5_1_ROM = squeeze(out.sub2_5_1_ROM.Data);
sub2_5_1_linFOM = squeeze(out.sub2_5_1_linFOM.Data);
sub2_5_1 = squeeze(out.sub2_5_1.Data);
sub2_5_1_ROM1 = squeeze(out.sub2_5_1_ROM1.Data);
sub2_5_1_linFOM1 = squeeze(out.sub2_5_1_linFOM1.Data);

sub3_6_ROM = squeeze(out.sub3_6_ROM.Data);
sub3_6_linFOM = squeeze(out.sub3_6_linFOM.Data);
sub3_6 = squeeze(out.sub3_6_norm.Data);
sub3_6_ROM1 = squeeze(out.sub3_6_ROM1.Data);
sub3_6_linFOM1 = squeeze(out.sub3_6_linFOM1.Data);


fig = figure();
clf(fig);

t = tiledlayout(3, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% Top tile: Input Power
ax.Pin = nexttile(t);
hold(ax.Pin, 'on');
grid(ax.Pin, 'on');


V_GFL1 = squeeze(out.P_GFL1.Data);
V_GFL1_LinFOM = squeeze(out.P_GFL1_linFOM.Data);
V_GFL1_LinROM = squeeze(out.P_GFL1_linROM.Data);


plot(ax.Pin, tout, V_GFL1, 'b-', 'linewidth', 2);
plot(ax.Pin, tout, V_GFL1_LinFOM, 'k:',  'linewidth', 2);
plot(ax.Pin, tout, V_GFL1_LinROM,    'r--',  'linewidth', 2);

ylabel(ax.Pin, 'Active power (p.u.)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

legend(ax.Pin, ...
    '$\mathrm{FOM}_{\mathrm{NL-NL}}$', ...
    '$\mathrm{FOM}_{\mathrm{Lin-Lin}}$', ...
    '$\mathrm{ROM}_{\mathrm{Lin-Lin}}$', ...
    'FontSize', 14, ...
    'Interpreter', 'latex');

% Remove top x tick labels for clean stacking
xticklabels(ax.Pin, []);
hold(ax.Pin, 'off');


% Middle tile: Active Power Responses
ax.P = nexttile(t);
hold(ax.P, 'on');
grid(ax.P, 'on');

P_ROM = squeeze(out.P_ROM.Data);
P_linFOM = squeeze(out.P_linFOM.Data);
P_nonlinFOM = squeeze(out.P_nonlinFOM.Data);
P_ROM1 = squeeze(out.P_ROM1.Data);
P_linFOM1 = squeeze(out.P_linFOM1.Data);

plot(ax.P, tout, P_nonlinFOM, 'b-', 'linewidth', 2);
plot(ax.P, tout, P_linFOM1, 'k:',  'linewidth', 2);
plot(ax.P, tout, P_ROM1,    'r--', 'linewidth', 2);

ylabel(ax.P, 'Active power (p.u.)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

legend(ax.P, ...
    '$\mathrm{FOM}_\mathrm{NL-NL}$', ...
    '$\mathrm{FOM}_\mathrm{Lin-Lin}$', ...
    '$\mathrm{ROM}_\mathrm{Lin-Lin}$', ...
    'FontSize', 14, ...
    'Interpreter', 'latex', ...
    'Location', 'best');

% Remove middle x tick labels for clean stacking
xticklabels(ax.P, []);
hold(ax.P, 'off');


% Bottom tile: Power error
ax.Pe = nexttile(t);
hold(ax.Pe, 'on');
grid(ax.Pe, 'on');

% Note: Since these are multiplied by 100, they are percentages.
Perr_nonlin = abs(P_ROM1 - P_nonlinFOM) ./ abs(P_nonlinFOM) .* 100;
Perr_lin    = abs(P_ROM1 - P_linFOM1)    ./ abs(P_linFOM1)    .* 100;

plot(ax.Pe, tout, Perr_nonlin, 'b-', 'linewidth', 2);
plot(ax.Pe, tout, Perr_lin,    'k:', 'linewidth', 2);

xlabel(ax.Pe, 'Time (s)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

% Added (%) to the label for clarity since data is multiplied by 100
ylabel(ax.Pe, 'Relative error (\%)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

legend(ax.Pe, ...
    'vs. $\mathrm{FOM}_\mathrm{NL-NL}$', ...
    'vs. $\mathrm{FOM}_\mathrm{Lin-Lin}$', ...
    'FontSize', 14, ...
    'Interpreter', 'latex', ...
    'Location', 'best');

hold(ax.Pe, 'off');


% Final Formatting
% Link all three x-axes so zooming works perfectly across the column
linkaxes([ax.Pin, ax.P, ax.Pe], 'x');

% Adjust figure height to accommodate 3 rows (increased from 4.6 to 6.5)
set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 6.5]);

% Optional: Use exportgraphics for perfect PDF LaTeX rendering
% drawnow;
% exportgraphics(fig, 'Power_Response_Plot.pdf', 'ContentType', 'vector');



%% Fig. 10: Self-clearing inverter blocking
%%% Additional settings
selectedGFL = 2;

Ttrip_gfl  = 0.1;
Tclear_gfl = 0.11;

numberOfGFLs = sub(3).n_gfl;

gfl3_enable_before = ones(numberOfGFLs,1);

gfl3_enable_blocked = gfl3_enable_before;
gfl3_enable_blocked(selectedGFL) = 0;

gfl3_enable_zero = zeros(numberOfGFLs,1);

gfl3_enable_clear_increment = ...
    gfl3_enable_before - gfl3_enable_blocked;

gfl3_enable_after = gfl3_enable_blocked;

%%% Plot
out = sim("subsystem_118_bus_self_clearing_blocking");

tout = squeeze(out.tout);
VBus72 = squeeze(out.self_clear_scope8_1.Data);
VBus72_LinFOM = squeeze(out.self_clear_scope8_2.Data);
VBus72_LinROM = squeeze(out.self_clear_scope8_3.Data);

fig = figure();
clf;

t = tiledlayout(2, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax.P = nexttile(t);
hold(ax.P, 'on');
grid(ax.P, 'on');

plot(ax.P, tout, VBus72, 'b-', 'linewidth', 2);
plot(ax.P, tout, VBus72_LinFOM, 'k:',  'linewidth', 2);
plot(ax.P, tout, VBus72_LinROM,    'r--',  'linewidth', 2);

ylabel(ax.P, 'Voltage (p.u.)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

legend(ax.P, ...
    '$\mathrm{FOM}_{\mathrm{NL-NL}}$', ...
    '$\mathrm{FOM}_{\mathrm{Lin-NL}}$', ...
    '$\mathrm{ROM}_{\mathrm{Lin-NL}}$', ...
    'FontSize', 14, ...
    'Interpreter', 'latex');

% Remove top x tick labels for compactness
xticklabels(ax.P, []);

hold(ax.P, 'off');

%%%
% Add zoomed inset to the top axes
zoomInterval = [0.1, 0.15];

% Ensure axes positions have been calculated by tiledlayout
drawnow;

% Identify samples inside the zoom interval
idxZoom = tout >= zoomInterval(1) & tout <= zoomInterval(2);

if any(idxZoom)

    % Position the inset relative to the upper axes
    mainPos = ax.P.Position;

    insetPos = [ ...
        mainPos(1) + 0.57*mainPos(3), ...
        mainPos(2) + 0.43*mainPos(4), ...
        0.38*mainPos(3), ...
        0.48*mainPos(4)];

    axInset = axes(fig, 'Position', insetPos);
    hold(axInset, 'on');
    grid(axInset, 'on');
    box(axInset, 'on');

    plot(axInset, tout, VBus72, ...
        'b-', 'LineWidth', 1.8);

    plot(axInset, tout, VBus72_LinFOM, ...
        'k:', 'LineWidth', 1.8);

    plot(axInset, tout, VBus72_LinROM, ...
        'r--', 'LineWidth', 1.8);

    xlim(axInset, zoomInterval);

    % Automatically choose a useful vertical range
    yZoom = [ ...
        VBus72(idxZoom); ...
        VBus72_LinFOM(idxZoom); ...
        VBus72_LinROM(idxZoom)];

    yMin = min(yZoom, [], 'all');
    yMax = max(yZoom, [], 'all');

    yRange = yMax - yMin;

    if yRange > 0
        yMargin = 0.10*yRange;
    else
        yMargin = max(abs(yMin)*0.01, 1e-4);
    end

    ylim(axInset, [yMin-yMargin, yMax+yMargin]);

    axInset.FontSize = 12;
    axInset.TickLabelInterpreter = 'latex';

    % xlabel(axInset, '$t$ (s)', ...
    %     'Interpreter', 'latex', ...
    %     'FontSize', 10);
    % 
    % ylabel(axInset, '$V$ (p.u.)', ...
    %     'Interpreter', 'latex', ...
    %     'FontSize', 10);

    hold(axInset, 'off');

    % Mark the enlarged interval on the main upper plot
    rectangle(ax.P, ...
        'Position', [ ...
            zoomInterval(1), ...
            yMin-yMargin, ...
            diff(zoomInterval), ...
            yRange+2*yMargin], ...
        'LineStyle', '--', ...
        'LineWidth', 1.0, ...
        'EdgeColor', [0.35, 0.35, 0.35]);

end
%%%


ax.Pe = nexttile(t);
hold(ax.Pe, 'on');
grid(ax.Pe, 'on');

Perr_nonlin = abs(VBus72_LinROM - VBus72)./abs(VBus72).*100;
Perr_lin    = abs(VBus72_LinROM - VBus72_LinFOM)./abs(VBus72_LinFOM).*100;

plot(ax.Pe, tout, Perr_nonlin, 'b-', 'linewidth', 2);
plot(ax.Pe, tout, Perr_lin,    'k:', 'linewidth', 2);

xlabel(ax.Pe, 'Time (s)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

ylabel(ax.Pe, 'Relative error ($\%$)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

legend(ax.Pe, ...
    'vs. $\mathrm{FOM}_{\mathrm{NL-NL}}$', ...
    'vs. $\mathrm{FOM}_{\mathrm{Lin-NL}}$', ...
    'FontSize', 14, ...
    'Interpreter', 'latex');

hold(ax.Pe, 'off');

% Link x-axes
linkaxes([ax.P, ax.Pe], 'x');

% Optional figure size for two-column-width publication figure
set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);
