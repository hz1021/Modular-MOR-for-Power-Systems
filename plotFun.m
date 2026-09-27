%%
load nominalU.mat;
load nominalY.mat;
load nominalInter.mat

Tstop = 1;
%%
for P_input = 0.03 % [0.01, 0.05, 0.1, 0.5, 1]
    out = sim("subsystem_118_bus");
    %% Extract time-domain data from .slx simulation
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


    %% Power
    fig = figure();
    subplot(2, 1, 1);
    hold on;
    grid on;


    plot(tout, P_ref, 'g-.', 'linewidth', 1.6);
    plot(tout, P_ROM, 'r--', 'linewidth', 1.6);
    plot(tout, P_nonlinFOM, 'b-', 'linewidth', 1.6);
    plot(tout, P_linFOM, 'k:', 'linewidth', 1.6);
    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Power (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('$P_{ROM}$', '$P_{nonlinFOM}$', '$P_{linFOM}$', 'input', 'FontSize', 13.5, 'Interpreter', 'latex');
    title("Input " + P_input);

    hold off;

    %% Power error
    subplot(2, 1, 2);
    hold on;
    grid on;

    plot(tout, abs(P_ROM - P_nonlinFOM)./abs(P_nonlinFOM).*100, 'g:', 'linewidth', 1.6);
    plot(tout, abs(P_ROM - P_linFOM)./abs(P_linFOM).*100, 'k:', 'linewidth', 1.6);
    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Power Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('vs. $P_{nonlinFOM}$', 'vs. $P_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    grid on;
    title("Input " + P_input);
    % ytickformat('percentage');

    hold off;

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    %%
    %% Power and power error: compact version

    fig = figure();
    clf;
    
    t = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Power
    % ============================================================
    ax.P = nexttile(t);
    hold(ax.P, 'on');
    grid(ax.P, 'on');
    
    % plot(ax.P, tout, P_input, 'g-.', 'linewidth', 2);
    plot(ax.P, tout, P_nonlinFOM, 'b-', 'linewidth', 2);
    plot(ax.P, tout, P_linFOM, 'k:',  'linewidth', 2);
    plot(ax.P, tout, P_ROM,    'r--',  'linewidth', 2);
    
    ylabel(ax.P, 'Active power (p.u.)', ...
        'FontSize', 18, ...
        'Interpreter', 'latex');
    
    legend(ax.P, ...
        '$\mathrm{FOM}_\mathrm{NL-NL}$', ...
        '$\mathrm{FOM}_\mathrm{Lin-NL}$', ...
        '$\mathrm{ROM}_\mathrm{Lin-NL}$', ...
        'FontSize', 14, ...
        'Interpreter', 'latex', ...
        'Location', 'best');
    
    % title(ax.P, "Input " + P_input);
    
    % Remove top x tick labels for compactness
    xticklabels(ax.P, []);
    
    hold(ax.P, 'off');
    
    
    % ============================================================
    % Bottom tile: Power error
    % ============================================================
    ax.Pe = nexttile(t);
    hold(ax.Pe, 'on');
    grid(ax.Pe, 'on');
    
    Perr_nonlin = abs(P_ROM - P_nonlinFOM) ./ abs(P_nonlinFOM) .* 100;
    Perr_lin    = abs(P_ROM - P_linFOM)    ./ abs(P_linFOM)    .* 100;
    
    plot(ax.Pe, tout, Perr_nonlin, 'b-', 'linewidth', 2);
    plot(ax.Pe, tout, Perr_lin,    'k:', 'linewidth', 2);
    
    xlabel(ax.Pe, 'Time (s)', ...
        'FontSize', 18, ...
        'Interpreter', 'latex');
    
    ylabel(ax.Pe, 'Relative error', ...
        'FontSize', 18, ...
        'Interpreter', 'latex');
    
    legend(ax.Pe, ...
        'vs. $\mathrm{FOM}_\mathrm{NL-NL}$', ...
        'vs. $\mathrm{FOM}_\mathrm{Lin-NL}$', ...
        'FontSize', 14, ...
        'Interpreter', 'latex', ...
        'Location', 'best');
    
    hold(ax.Pe, 'off');
    
    % Link x-axes
    linkaxes([ax.P, ax.Pe], 'x');

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);
    

    %%
    %% Power, Response, and Error: 3-Tile Version
    fig = figure();
    clf(fig);
    
    t = tiledlayout(3, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Input Power
    % ============================================================
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
        '$\mathrm{FOM}_{\mathrm{Lin-NL}}$', ...
        '$\mathrm{ROM}_{\mathrm{Lin-NL}}$', ...
        'FontSize', 14, ...
        'Interpreter', 'latex');
    
    % Remove top x tick labels for clean stacking
    xticklabels(ax.Pin, []);
    hold(ax.Pin, 'off');
    
    % ============================================================
    % Middle tile: Active Power Responses
    % ============================================================
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
    
    % ============================================================
    % Bottom tile: Power error
    % ============================================================
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
    
    % ============================================================
    % Final Formatting
    % ============================================================
    % Link all three x-axes so zooming works perfectly across the column
    linkaxes([ax.Pin, ax.P, ax.Pe], 'x');
    
    % Adjust figure height to accommodate 3 rows (increased from 4.6 to 6.5)
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 6.5]);
    
    % Optional: Use exportgraphics for perfect PDF LaTeX rendering
    % drawnow;
    % exportgraphics(fig, 'Power_Response_Plot.pdf', 'ContentType', 'vector');

    %% sub3/5_1: Vout1
    fig = figure();
    subplot(2, 1, 1);
    hold on;
    grid on;

    plot(tout, sub3_5_1_ROM, 'r--', 'linewidth', 1.6);
    plot(tout, sub3_5_1, 'b-', 'linewidth', 1.6);
    plot(tout, sub3_5_1_linFOM, 'k:', 'linewidth', 1.6);

    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Voltage (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('$V_{ROM}$', '$V_{nonlinFOM}$', '$V_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    title("Input " + P_input);
    
    hold off;
    
    %% sub3/5_1 error
    subplot(2, 1, 2);
    hold on;

    plot(tout, abs(sub3_5_1_ROM - sub3_5_1)./abs(sub3_5_1).*100, 'g:', 'linewidth', 1.6);
    plot(tout, abs(sub3_5_1_ROM - sub3_5_1_linFOM)./abs(sub3_5_1_linFOM).*100, 'k:', 'linewidth', 1.6);
    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Voltage Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('vs. $V_{nonlinFOM}$ error', 'vs. $V_{linFOM}$ error', 'FontSize', 13.5, 'Interpreter', 'latex');
    grid on;
    title("Input " + P_input);
    % ytickformat('percentage');

    hold off;

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    %%
    fig = figure();
    clf;
    
    t = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Vout1
    % ============================================================
    ax.Vout1 = nexttile(t);
    hold(ax.Vout1, 'on');
    grid(ax.Vout1, 'on');
    
    plot(ax.Vout1, tout, sub3_5_1,    'b-', 'linewidth', 1.6);
    plot(ax.Vout1, tout, sub3_5_1_linFOM,        'k:',  'linewidth', 1.6);
    plot(ax.Vout1, tout, sub3_5_1_ROM, 'r--',  'linewidth', 1.6);
    
    ylabel(ax.Vout1, 'Voltage (p.u.)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.Vout1, ...
        'Nonlinear FOM', ...
        'Linear FOM', ...
        'ROM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    % title(ax.Vout1, "Input " + P_input);
    
    % Hide top x-axis tick labels for compactness
    xticklabels(ax.Vout1, []);
    
    hold(ax.Vout1, 'off');
    
    
    % ============================================================
    % Bottom tile: Vout1 error
    % ============================================================
    ax.Vout1Err = nexttile(t);
    hold(ax.Vout1Err, 'on');
    grid(ax.Vout1Err, 'on');
    
    Verr_nonlin = abs(sub3_5_1_ROM - sub3_5_1) ./ abs(sub3_5_1) .* 100;
    Verr_lin    = abs(sub3_5_1_ROM - sub3_5_1_linFOM) ./ abs(sub3_5_1_linFOM) .* 100;
    
    plot(ax.Vout1Err, tout, Verr_nonlin, 'g:', 'linewidth', 1.6);
    plot(ax.Vout1Err, tout, Verr_lin,    'k:', 'linewidth', 1.6);
    
    xlabel(ax.Vout1Err, 'Time ($s$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    ylabel(ax.Vout1Err, 'Comparison error ($\%$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.Vout1Err, ...
        'vs. Nonlinear FOM', ...
        'vs. Linear FOM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    hold(ax.Vout1Err, 'off');
    
    % Link x-axes
    linkaxes([ax.Vout1, ax.Vout1Err], 'x');

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    %% sub1/6_1: I
    fig = figure();
    
    subplot(2, 1, 1);
    hold on;
    grid on;

    plot(tout, sub1_6_1_ROM, 'r--', 'linewidth', 1.6);
    plot(tout, sub1_6_1, 'b-', 'linewidth', 1.6);
    plot(tout, sub1_6_1_linFOM, 'k:', 'linewidth', 1.6);

    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Current (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('$I_{ROM}$', '$I_{nonlinFOM}$', '$I_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    title("Input " + P_input);
    
    hold off;
    
    %% sub1/6_1 error
    subplot(2, 1, 2);
    hold on;

    plot(tout, abs(sub1_6_1_ROM - sub1_6_1)./abs(sub1_6_1).*100, 'g:', 'linewidth', 1.6);
    plot(tout, abs(sub1_6_1_ROM - sub1_6_1_linFOM)./abs(sub1_6_1_linFOM).*100, 'k:', 'linewidth', 1.6);
    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Current Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('vs. $I_{nonlinFOM}$ error', 'vs. $I_{linFOM}$ error', 'FontSize', 13.5, 'Interpreter', 'latex');
    grid on;
    title("Input " + P_input);
    % ytickformat('percentage');

    hold off;

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    %% sub2/5_1: Vout
    fig = figure();
    subplot(2, 1, 1);
    hold on;
    grid on;


    plot(tout, sub2_5_1_ROM, 'r--', 'linewidth', 1.6);
    plot(tout, sub2_5_1, 'b-', 'linewidth', 1.6);
    plot(tout, sub2_5_1_linFOM, 'k:', 'linewidth', 1.6);

    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Voltage (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('$V_{ROM}$', '$V_{nonlinFOM}$', '$V_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    title("Input " + P_input);
    
    hold off;
    
    %% sub2/5_1 error
    subplot(2, 1, 2);
    hold on;

    plot(tout, abs(sub2_5_1_ROM - sub2_5_1)./abs(sub2_5_1).*100, 'g:', 'linewidth', 1.6);
    plot(tout, abs(sub2_5_1_ROM - sub2_5_1_linFOM)./abs(sub2_5_1_linFOM).*100, 'k:', 'linewidth', 1.6);
    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Voltage Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('vs. $V_{nonlinFOM}$ error', 'vs. $V_{linFOM}$ error', 'FontSize', 13.5, 'Interpreter', 'latex');
    grid on;
    title("Input " + P_input);
    % ytickformat('percentage');

    hold off;

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

        %% 2_5 error: compact version

    fig = figure();
    clf;
    
    t = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Power
    % ============================================================
    ax.P = nexttile(t);
    hold(ax.P, 'on');
    grid(ax.P, 'on');

    plot(ax.P, tout, sub2_5_1_ROM, 'b-', 'linewidth', 1.6);
    plot(ax.P, tout, sub2_5_1, 'k:',  'linewidth', 1.6);
    plot(ax.P, tout, sub2_5_1_linFOM,    'r--',  'linewidth', 1.6);
    
    ylabel(ax.P, 'Voltage (p.u.)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.P, ...
        'Nonlinear FOM', ...
        'Linear FOM', ...
        'ROM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    % title(ax.P, "Input " + P_input);
    
    % Remove top x tick labels for compactness
    xticklabels(ax.P, []);
    
    hold(ax.P, 'off');
    
    
    % ============================================================
    % Bottom tile: Power error
    % ============================================================
    ax.Pe = nexttile(t);
    hold(ax.Pe, 'on');
    grid(ax.Pe, 'on');

    Perr_nonlin = abs(sub2_5_1_ROM - sub2_5_1)./abs(sub2_5_1).*100;
    Perr_lin    = abs(sub2_5_1_ROM - sub2_5_1_linFOM)./abs(sub2_5_1_linFOM).*100;
    
    plot(ax.Pe, tout, Perr_nonlin, 'g:', 'linewidth', 1.6);
    plot(ax.Pe, tout, Perr_lin,    'k:', 'linewidth', 1.6);
    
    xlabel(ax.Pe, 'Time ($s$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    ylabel(ax.Pe, 'Comparison error ($\%$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.Pe, ...
        'vs. Nonlinear FOM', ...
        'vs. Linear FOM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    hold(ax.Pe, 'off');
    
    % Link x-axes
    linkaxes([ax.P, ax.Pe], 'x');

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    %% sub3/6: I_PCC
    fig = figure();
    subplot(2, 1, 1);
    hold on;
    grid on;


    plot(tout, sub3_6_ROM, 'r--', 'linewidth', 1.6);
    plot(tout, sub3_6, 'b-', 'linewidth', 1.6);
    plot(tout, sub3_6_linFOM, 'k:', 'linewidth', 1.6);

    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Current PCC (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('$I\_PCC_{ROM}$', '$I\_PCC_{nonlinFOM}$', '$I\_PCC_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    title("Input " + P_input);
    
    hold off;
    
    %% sub3/6: I_PCC error
    subplot(2, 1, 2);
    hold on;

    plot(tout, abs(sub3_6_ROM - sub3_6)./abs(sub3_6).*100, 'g:', 'linewidth', 1.6);
    plot(tout, abs(sub3_6_ROM - sub3_6_linFOM)./abs(sub3_6_linFOM).*100, 'k:', 'linewidth', 1.6);
    xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    ylabel('Current PCC Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    legend('vs. $I\_PCC_{nonlinFOM}$ error', 'vs. $I\_PCC_{linFOM}$ error', 'FontSize', 13.5, 'Interpreter', 'latex');
    grid on;
    title("Input " + P_input);
    % ytickformat('percentage');

    hold off;

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %% Power
    % figure();
    % subplot(2, 1, 1);
    % hold on;
    % grid on;
    % 
    % plot(tout, P_ROM1, 'r--', 'linewidth', 1.6);
    % plot(tout, P_nonlinFOM, 'b-', 'linewidth', 1.6);
    % plot(tout, P_linFOM1, 'k:', 'linewidth', 1.6);
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Power (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('$P_{ROM}$', '$P_{nonlinFOM}$', '$P_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    % title("Input " + P_input);
    % 
    % hold off;
    % 
    % %% Power error
    % subplot(2, 1, 2);
    % hold on;
    % grid on;
    % 
    % plot(tout, abs(P_ROM1 - P_nonlinFOM)./abs(P_nonlinFOM).*100, 'g:', 'linewidth', 1.6);
    % plot(tout, abs(P_ROM1 - P_linFOM1)./abs(P_linFOM1).*100, 'k:', 'linewidth', 1.6);
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Power Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('vs. $P_{nonlinFOM}$', 'vs. $P_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    % grid on;
    % title("Input " + P_input);
    % % ytickformat('percentage');
    % 
    % hold off;

    %%
    %% Power and power error: compact version

    fig = figure();
    clf;
    
    t = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Power
    % ============================================================
    ax.P = nexttile(t);
    hold(ax.P, 'on');
    grid(ax.P, 'on');
    
    plot(ax.P, tout, P_nonlinFOM, 'b-', 'linewidth', 2);
    plot(ax.P, tout, P_linFOM1, 'k:',  'linewidth', 2);
    plot(ax.P, tout, P_ROM1,    'r--',  'linewidth', 2);
    
    ylabel(ax.P, 'Active power (p.u.)', ...
        'FontSize', 18, ...
        'Interpreter', 'latex');
    
    legend(ax.P, ...
        'NL FOM + NL $G_0$', ...
        'Lin FOM + Lin $G_0$', ...
        'ROM + Lin $G_0$', ...
        'FontSize', 14, ...
        'Interpreter', 'latex');
    
    % title(ax.P, "Input " + P_input);
    
    % Remove top x tick labels for compactness
    xticklabels(ax.P, []);
    
    hold(ax.P, 'off');
    
    
    % ============================================================
    % Bottom tile: Power error
    % ============================================================
    ax.Pe = nexttile(t);
    hold(ax.Pe, 'on');
    grid(ax.Pe, 'on');
    
    Perr_nonlin = abs(P_ROM1 - P_nonlinFOM) ./ abs(P_nonlinFOM) .* 100;
    Perr_lin    = abs(P_ROM1 - P_linFOM1)    ./ abs(P_linFOM1)    .* 100;
    
    plot(ax.Pe, tout, Perr_nonlin, 'g:', 'linewidth', 2);
    plot(ax.Pe, tout, Perr_lin,    'k:', 'linewidth', 2);
    
    xlabel(ax.Pe, 'Time (s)', ...
        'FontSize', 18, ...
        'Interpreter', 'latex');
    
    ylabel(ax.Pe, 'Comparison error', ...
        'FontSize', 18, ...
        'Interpreter', 'latex');
    
    legend(ax.Pe, ...
        'vs. NL FOM', ...
        'vs. Lin FOM', ...
        'FontSize', 14, ...
        'Interpreter', 'latex');
    
    hold(ax.Pe, 'off');
    
    % Link x-axes
    linkaxes([ax.P, ax.Pe], 'x');

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);
    
    % %% sub3/5_1: Vout1
    % figure();
    % subplot(2, 1, 1);
    % hold on;
    % grid on;
    % 
    % plot(tout, sub3_5_1_ROM1, 'r--', 'linewidth', 1.6);
    % plot(tout, sub3_5_1, 'b-', 'linewidth', 1.6);
    % plot(tout, sub3_5_1_linFOM1, 'k:', 'linewidth', 1.6);
    % 
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Voltage (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('$V_{ROM}$', '$V_{nonlinFOM}$', '$V_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    % title("Input " + P_input);
    % 
    % hold off;
    % 
    % %% sub3/5_1 error
    % subplot(2, 1, 2);
    % hold on;
    % 
    % plot(tout, abs(sub3_5_1_ROM1 - sub3_5_1)./abs(sub3_5_1).*100, 'g:', 'linewidth', 1.6);
    % plot(tout, abs(sub3_5_1_ROM1 - sub3_5_1_linFOM1)./abs(sub3_5_1_linFOM1).*100, 'k:', 'linewidth', 1.6);
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Voltage Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('vs. $V_{nonlinFOM}$ error', 'vs. $V_{linFOM}$ error', 'FontSize', 13.5, 'Interpreter', 'latex');
    % grid on;
    % title("Input " + P_input);
    % % ytickformat('percentage');
    % 
    % hold off;

    %%
    fig = figure();
    clf;
    
    t = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Vout1
    % ============================================================
    ax.Vout1 = nexttile(t);
    hold(ax.Vout1, 'on');
    grid(ax.Vout1, 'on');
    
    plot(ax.Vout1, tout, sub3_5_1,    'b-', 'linewidth', 1.6);
    plot(ax.Vout1, tout, sub3_5_1_linFOM1,        'k:',  'linewidth', 1.6);
    plot(ax.Vout1, tout, sub3_5_1_ROM1, 'r--',  'linewidth', 1.6);
    
    ylabel(ax.Vout1, 'Voltage (p.u.)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.Vout1, ...
        'Nonlinear FOM', ...
        'Linear FOM', ...
        'ROM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    % title(ax.Vout1, "Input " + P_input);
    
    % Hide top x-axis tick labels for compactness
    xticklabels(ax.Vout1, []);
    
    hold(ax.Vout1, 'off');
    
    
    % ============================================================
    % Bottom tile: Vout1 error
    % ============================================================
    ax.Vout1Err = nexttile(t);
    hold(ax.Vout1Err, 'on');
    grid(ax.Vout1Err, 'on');
    
    Verr_nonlin = abs(sub3_5_1_ROM1 - sub3_5_1) ./ abs(sub3_5_1) .* 100;
    Verr_lin    = abs(sub3_5_1_ROM1 - sub3_5_1_linFOM1) ./ abs(sub3_5_1_linFOM1) .* 100;
    
    plot(ax.Vout1Err, tout, Verr_nonlin, 'g:', 'linewidth', 1.6);
    plot(ax.Vout1Err, tout, Verr_lin,    'k:', 'linewidth', 1.6);
    
    xlabel(ax.Vout1Err, 'Time ($s$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    ylabel(ax.Vout1Err, 'Comparison error ($\%$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.Vout1Err, ...
        'vs. Nonlinear FOM', ...
        'vs. Linear FOM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    hold(ax.Vout1Err, 'off');
    
    % Link x-axes
    linkaxes([ax.Vout1, ax.Vout1Err], 'x');

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    % %% sub1/6_1: I
    % figure();
    % 
    % subplot(2, 1, 1);
    % hold on;
    % grid on;
    % 
    % plot(tout, sub1_6_1_ROM1, 'r--', 'linewidth', 1.6);
    % plot(tout, sub1_6_1, 'b-', 'linewidth', 1.6);
    % plot(tout, sub1_6_1_linFOM1, 'k:', 'linewidth', 1.6);
    % 
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Current (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('$I_{ROM}$', '$I_{nonlinFOM}$', '$I_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    % title("Input " + P_input);
    % 
    % hold off;
    % 
    % %% sub1/6_1 error
    % subplot(2, 1, 2);
    % hold on;
    % 
    % plot(tout, abs(sub1_6_1_ROM1 - sub1_6_1)./abs(sub1_6_1).*100, 'g:', 'linewidth', 1.6);
    % plot(tout, abs(sub1_6_1_ROM1 - sub1_6_1_linFOM1)./abs(sub1_6_1_linFOM1).*100, 'k:', 'linewidth', 1.6);
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Current Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('vs. $I_{nonlinFOM}$ error', 'vs. $I_{linFOM}$ error', 'FontSize', 13.5, 'Interpreter', 'latex');
    % grid on;
    % title("Input " + P_input);
    % % ytickformat('percentage');
    % 
    % hold off;

    %%
    %% 1_6_1 error: compact version

    fig = figure();
    clf;
    
    t = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Power
    % ============================================================
    ax.P = nexttile(t);
    hold(ax.P, 'on');
    grid(ax.P, 'on');
    
    plot(ax.P, tout, sub1_6_1_ROM1, 'b-', 'linewidth', 1.6);
    plot(ax.P, tout, sub1_6_1, 'k:',  'linewidth', 1.6);
    plot(ax.P, tout, sub1_6_1_linFOM1,    'r--',  'linewidth', 1.6);
    
    ylabel(ax.P, 'Current (p.u.)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.P, ...
        'Nonlinear FOM', ...
        'Linear FOM', ...
        'ROM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    % title(ax.P, "Input " + P_input);
    
    % Remove top x tick labels for compactness
    xticklabels(ax.P, []);
    
    hold(ax.P, 'off');
    
    
    % ============================================================
    % Bottom tile: Power error
    % ============================================================
    ax.Pe = nexttile(t);
    hold(ax.Pe, 'on');
    grid(ax.Pe, 'on');

    
    Perr_nonlin = abs(sub1_6_1_ROM1 - sub1_6_1)./abs(sub1_6_1) .* 100;
    Perr_lin    = abs(sub1_6_1_ROM1 - sub1_6_1_linFOM1)./abs(sub1_6_1_linFOM1)    .* 100;
    
    plot(ax.Pe, tout, Perr_nonlin, 'g:', 'linewidth', 1.6);
    plot(ax.Pe, tout, Perr_lin,    'k:', 'linewidth', 1.6);
    
    xlabel(ax.Pe, 'Time ($s$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    ylabel(ax.Pe, 'Comparison error ($\%$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.Pe, ...
        'vs. Nonlinear FOM', ...
        'vs. Linear FOM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    hold(ax.Pe, 'off');
    
    % Link x-axes
    linkaxes([ax.P, ax.Pe], 'x');

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    % %% sub2/5_1: Vout
    % figure();
    % subplot(2, 1, 1);
    % hold on;
    % grid on;
    % 
    % 
    % plot(tout, sub2_5_1_ROM1, 'r--', 'linewidth', 1.6);
    % plot(tout, sub2_5_1, 'b-', 'linewidth', 1.6);
    % plot(tout, sub2_5_1_linFOM1, 'k:', 'linewidth', 1.6);
    % 
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Voltage (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('$V_{ROM}$', '$V_{nonlinFOM}$', '$V_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    % title("Input " + P_input);
    % 
    % hold off;
    % 
    % %% sub2/5_1 error
    % subplot(2, 1, 2);
    % hold on;
    % 
    % plot(tout, abs(sub2_5_1_ROM1 - sub2_5_1)./abs(sub2_5_1).*100, 'g:', 'linewidth', 1.6);
    % plot(tout, abs(sub2_5_1_ROM1 - sub2_5_1_linFOM1)./abs(sub2_5_1_linFOM1).*100, 'k:', 'linewidth', 1.6);
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Voltage Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('vs. $V_{nonlinFOM}$ error', 'vs. $V_{linFOM}$ error', 'FontSize', 13.5, 'Interpreter', 'latex');
    % grid on;
    % title("Input " + P_input);
    % % ytickformat('percentage');
    % 
    % hold off;

    %% 2_5 error: compact version

    fig = figure();
    clf;
    
    t = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Power
    % ============================================================
    ax.P = nexttile(t);
    hold(ax.P, 'on');
    grid(ax.P, 'on');

    plot(ax.P, tout, sub2_5_1_ROM1, 'b-', 'linewidth', 1.6);
    plot(ax.P, tout, sub2_5_1, 'k:',  'linewidth', 1.6);
    plot(ax.P, tout, sub2_5_1_linFOM1,    'r--',  'linewidth', 1.6);
    
    ylabel(ax.P, 'Voltage (p.u.)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.P, ...
        'Nonlinear FOM', ...
        'Linear FOM', ...
        'ROM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    % title(ax.P, "Input " + P_input);
    
    % Remove top x tick labels for compactness
    xticklabels(ax.P, []);
    
    hold(ax.P, 'off');
    
    
    % ============================================================
    % Bottom tile: Power error
    % ============================================================
    ax.Pe = nexttile(t);
    hold(ax.Pe, 'on');
    grid(ax.Pe, 'on');

    Perr_nonlin = abs(sub2_5_1_ROM1 - sub2_5_1)./abs(sub2_5_1).*100;
    Perr_lin    = abs(sub2_5_1_ROM1 - sub2_5_1_linFOM1)./abs(sub2_5_1_linFOM1).*100;
    
    plot(ax.Pe, tout, Perr_nonlin, 'g:', 'linewidth', 1.6);
    plot(ax.Pe, tout, Perr_lin,    'k:', 'linewidth', 1.6);
    
    xlabel(ax.Pe, 'Time ($s$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    ylabel(ax.Pe, 'Comparison error ($\%$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.Pe, ...
        'vs. Nonlinear FOM', ...
        'vs. Linear FOM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    hold(ax.Pe, 'off');
    
    % Link x-axes
    linkaxes([ax.P, ax.Pe], 'x');

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    % %% sub3/6: I_PCC
    % figure();
    % subplot(2, 1, 1);
    % hold on;
    % grid on;
    % 
    % 
    % plot(tout, sub3_6_ROM1, 'r--', 'linewidth', 1.6);
    % plot(tout, sub3_6, 'b-', 'linewidth', 1.6);
    % plot(tout, sub3_6_linFOM1, 'k:', 'linewidth', 1.6);
    % 
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Current PCC (p.u.)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('$I\_PCC_{ROM}$', '$I\_PCC_{nonlinFOM}$', '$I\_PCC_{linFOM}$', 'FontSize', 13.5, 'Interpreter', 'latex');
    % title("Input " + P_input);
    % 
    % hold off;
    % 
    % %% sub3/6: I_PCC error
    % subplot(2, 1, 2);
    % hold on;
    % 
    % plot(tout, abs(sub3_6_ROM1 - sub3_6)./abs(sub3_6).*100, 'g:', 'linewidth', 1.6);
    % plot(tout, abs(sub3_6_ROM1 - sub3_6_linFOM1)./abs(sub3_6_linFOM1).*100, 'k:', 'linewidth', 1.6);
    % xlabel('Time ($s$)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % ylabel('Current PCC Error (%)', 'FontSize', 13.5, 'Interpreter', 'latex');
    % legend('vs. $I\_PCC_{nonlinFOM}$ error', 'vs. $I\_PCC_{linFOM}$ error', 'FontSize', 13.5, 'Interpreter', 'latex');
    % grid on;
    % title("Input " + P_input);
    % % ytickformat('percentage');
    % 
    % hold off;

    %% 3_6 error: compact version

    fig = figure();
    clf;
    
    t = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Power
    % ============================================================
    ax.P = nexttile(t);
    hold(ax.P, 'on');
    grid(ax.P, 'on');

    plot(ax.P, tout, sub3_6_ROM1, 'b-', 'linewidth', 1.6);
    plot(ax.P, tout, sub3_6, 'k:',  'linewidth', 1.6);
    plot(ax.P, tout, sub3_6_linFOM1,    'r--',  'linewidth', 1.6);
    
    ylabel(ax.P, 'Current (p.u.)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.P, ...
        'Nonlinear FOM', ...
        'Linear FOM', ...
        'ROM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    % title(ax.P, "Input " + P_input);
    
    % Remove top x tick labels for compactness
    xticklabels(ax.P, []);
    
    hold(ax.P, 'off');
    
    
    % ============================================================
    % Bottom tile: Power error
    % ============================================================
    ax.Pe = nexttile(t);
    hold(ax.Pe, 'on');
    grid(ax.Pe, 'on');

    Perr_nonlin = abs(sub3_6_ROM1 - sub3_6)./abs(sub3_6).*100;
    Perr_lin    = abs(sub3_6_ROM1 - sub3_6_linFOM1)./abs(sub3_6_linFOM1).*100;
    
    plot(ax.Pe, tout, Perr_nonlin, 'g:', 'linewidth', 1.6);
    plot(ax.Pe, tout, Perr_lin,    'k:', 'linewidth', 1.6);
    
    xlabel(ax.Pe, 'Time ($s$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    ylabel(ax.Pe, 'Comparison error ($\%$)', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    legend(ax.Pe, ...
        'vs. Nonlinear FOM', ...
        'vs. Linear FOM', ...
        'FontSize', 13.5, ...
        'Interpreter', 'latex');
    
    hold(ax.Pe, 'off');
    
    % Link x-axes
    linkaxes([ax.P, ax.Pe], 'x');

    % Optional figure size for two-column-width publication figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 4.6]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%%
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


%%
tout = squeeze(out.tout);
V_GFL1 = squeeze(out.P_GFL1.Data); % P_GFL1, sub3_5_1, sub1_6_1, sub2_5_1, sub3_6, P_nonlinFOM
V_GFL1_LinFOM = squeeze(out.P_GFL1_linFOM.Data);
V_GFL1_LinROM = squeeze(out.P_GFL1_linROM.Data);

% V_GFL1 = squeeze(out.sub3_5_1.Data); % P_GFL1
% V_GFL1_LinFOM = squeeze(out.sub3_5_1_linFOM.Data);
% V_GFL1_LinROM = squeeze(out.sub3_5_1_ROM.Data);
% 
% V_GFL1 = squeeze(out.sub3_6_norm.Data); % P_GFL1
% V_GFL1_LinFOM = squeeze(out.sub3_6_linFOM.Data);
% V_GFL1_LinROM = squeeze(out.sub3_6_ROM.Data);

V_GFL1 = squeeze(out. P_nonlinFOM.Data);
V_GFL1_LinFOM = squeeze(out. P_linFOM.Data);
V_GFL1_LinROM = squeeze(out. P_ROM.Data);


fig = figure();
clf;

t = tiledlayout(2, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax.P = nexttile(t);
hold(ax.P, 'on');
grid(ax.P, 'on');

plot(ax.P, tout, V_GFL1, 'b-', 'linewidth', 2);
plot(ax.P, tout, V_GFL1_LinFOM, 'k:',  'linewidth', 2);
plot(ax.P, tout, V_GFL1_LinROM,    'r--',  'linewidth', 2);

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

% %%%
% % Add zoomed inset to the top axes
% zoomInterval = [0.1, 0.15];
% 
% % Ensure axes positions have been calculated by tiledlayout
% drawnow;
% 
% % Identify samples inside the zoom interval
% idxZoom = tout >= zoomInterval(1) & tout <= zoomInterval(2);
% 
% if any(idxZoom)
% 
%     % Position the inset relative to the upper axes
%     mainPos = ax.P.Position;
% 
%     insetPos = [ ...
%         mainPos(1) + 0.57*mainPos(3), ...
%         mainPos(2) + 0.43*mainPos(4), ...
%         0.38*mainPos(3), ...
%         0.48*mainPos(4)];
% 
%     axInset = axes(fig, 'Position', insetPos);
%     hold(axInset, 'on');
%     grid(axInset, 'on');
%     box(axInset, 'on');
% 
%     plot(axInset, tout, V_GFL1, ...
%         'b-', 'LineWidth', 1.8);
% 
%     plot(axInset, tout, V_GFL1_LinFOM, ...
%         'k:', 'LineWidth', 1.8);
% 
%     plot(axInset, tout, V_GFL1_LinROM, ...
%         'r--', 'LineWidth', 1.8);
% 
%     xlim(axInset, zoomInterval);
% 
%     % Automatically choose a useful vertical range
%     yZoom = [ ...
%         VBus72(idxZoom); ...
%         VBus72_LinFOM(idxZoom); ...
%         VBus72_LinROM(idxZoom)];
% 
%     yMin = min(yZoom, [], 'all');
%     yMax = max(yZoom, [], 'all');
% 
%     yRange = yMax - yMin;
% 
%     if yRange > 0
%         yMargin = 0.10*yRange;
%     else
%         yMargin = max(abs(yMin)*0.01, 1e-4);
%     end
% 
%     ylim(axInset, [yMin-yMargin, yMax+yMargin]);
% 
%     axInset.FontSize = 12;
%     axInset.TickLabelInterpreter = 'latex';
% 
% 
%     hold(axInset, 'off');
% 
%     % Mark the enlarged interval on the main upper plot
%     rectangle(ax.P, ...
%         'Position', [ ...
%             zoomInterval(1), ...
%             yMin-yMargin, ...
%             diff(zoomInterval), ...
%             yRange+2*yMargin], ...
%         'LineStyle', '--', ...
%         'LineWidth', 1.0, ...
%         'EdgeColor', [0.35, 0.35, 0.35]);
% 
% end
% %%%


ax.Pe = nexttile(t);
hold(ax.Pe, 'on');
grid(ax.Pe, 'on');

Perr_nonlin = abs(V_GFL1_LinROM - V_GFL1)./abs(V_GFL1).*100;
Perr_lin    = abs(V_GFL1_LinROM - V_GFL1_LinFOM)./abs(V_GFL1_LinFOM).*100;

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



%% Power, Response, and Error: 3-Tile Version
    fig = figure();
    clf(fig);
    
    t = tiledlayout(3, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    
    % ============================================================
    % Top tile: Input Power
    % ============================================================
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
        '$\mathrm{FOM}_{\mathrm{Lin-NL}}$', ...
        '$\mathrm{ROM}_{\mathrm{Lin-NL}}$', ...
        'FontSize', 14, ...
        'Interpreter', 'latex');
    
    % Remove top x tick labels for clean stacking
    xticklabels(ax.Pin, []);
    hold(ax.Pin, 'off');
    
    % ============================================================
    % Middle tile: Active Power Responses
    % ============================================================
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
    
    % ============================================================
    % Bottom tile: Power error
    % ============================================================
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
    
    % ============================================================
    % Final Formatting
    % ============================================================
    % Link all three x-axes so zooming works perfectly across the column
    linkaxes([ax.Pin, ax.P, ax.Pe], 'x');
    
    % Adjust figure height to accommodate 3 rows (increased from 4.6 to 6.5)
    set(fig, 'Units', 'inches', 'Position', [1, 1, 7.16, 6.5]);
    
    % Optional: Use exportgraphics for perfect PDF LaTeX rendering
    % drawnow;
    % exportgraphics(fig, 'Power_Response_Plot.pdf', 'ContentType', 'vector');


    %% Power, Response, and Error: Error as Inset in Second Plot

fig = figure();
clf(fig);

% IMPORTANT:
% Keep the original 3x1 layout.
% This preserves the original size and gap of Fig. 1 and Fig. 2.
t = tiledlayout(fig, 3, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');


% ============================================================
% Top tile: Input Power
% ============================================================
ax.Pin = nexttile(t, 1);

hold(ax.Pin, 'on');
grid(ax.Pin, 'on');

V_GFL1        = squeeze(out.P_GFL1.Data);
V_GFL1_LinFOM = squeeze(out.P_GFL1_linFOM.Data);
V_GFL1_LinROM = squeeze(out.P_GFL1_linROM.Data);

plot(ax.Pin, tout, V_GFL1, ...
    'b-', 'LineWidth', 2);

plot(ax.Pin, tout, V_GFL1_LinFOM, ...
    'k:', 'LineWidth', 2);

plot(ax.Pin, tout, V_GFL1_LinROM, ...
    'r--', 'LineWidth', 2);

ylabel(ax.Pin, 'Active power (p.u.)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

legend(ax.Pin, ...
    '$\mathrm{FOM}_{\mathrm{NL-NL}}$', ...
    '$\mathrm{FOM}_{\mathrm{Lin-NL}}$', ...
    '$\mathrm{ROM}_{\mathrm{Lin-NL}}$', ...
    'FontSize', 14, ...
    'Interpreter', 'latex', ...
    'Location', 'best');

% Keep exactly as original
xticklabels(ax.Pin, []);

hold(ax.Pin, 'off');


% ============================================================
% Middle tile: Active Power Responses
% ============================================================
ax.P = nexttile(t, 2);

hold(ax.P, 'on');
grid(ax.P, 'on');

P_ROM       = squeeze(out.P_ROM.Data);
P_linFOM    = squeeze(out.P_linFOM.Data);
P_nonlinFOM = squeeze(out.P_nonlinFOM.Data);
P_ROM1      = squeeze(out.P_ROM1.Data);
P_linFOM1   = squeeze(out.P_linFOM1.Data);

h1 = plot(ax.P, tout, P_nonlinFOM, ...
    'b-', 'LineWidth', 2);

h2 = plot(ax.P, tout, P_linFOM1, ...
    'k:', 'LineWidth', 2);

h3 = plot(ax.P, tout, P_ROM1, ...
    'r--', 'LineWidth', 2);

ylabel(ax.P, 'Active power (p.u.)', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');

legend(ax.P, ...
    [h1 h2 h3], ...
    { ...
    '$\mathrm{FOM}_{\mathrm{NL-NL}}$', ...
    '$\mathrm{FOM}_{\mathrm{Lin-Lin}}$', ...
    '$\mathrm{ROM}_{\mathrm{Lin-Lin}}$' ...
    }, ...
    'FontSize', 14, ...
    'Interpreter', 'latex', ...
    'Location', 'best');

% Now this is the lowest main plot, so show the x tick labels
set(ax.P, 'XTickLabelMode', 'auto');

hold(ax.P, 'off');


% ============================================================
% Calculate relative errors
% ============================================================

Perr_nonlin = abs(P_ROM1 - P_nonlinFOM) ./ ...
              abs(P_nonlinFOM) .* 100;

Perr_lin = abs(P_ROM1 - P_linFOM1) ./ ...
           abs(P_linFOM1) .* 100;


% ============================================================
% Create ERROR INSET inside middle plot
% ============================================================

% First force MATLAB to finish positioning the tiled axes
drawnow;

% Position of the middle main axes
posP = ax.P.Position;

% ------------------------------------------------------------
% Inset size and position relative to middle axes
%
% [left bottom width height]
%
% Change these four numbers if you want to move/resize the inset.
% ------------------------------------------------------------

insetPos = [ ...
    posP(1) + 0.34 * posP(3), ...   % left
    posP(2) + 0.12 * posP(4), ...   % bottom
    0.42 * posP(3), ...             % width
    0.52 * posP(4)];                % height

ax.Pe = axes(fig, ...
    'Position', insetPos);

hold(ax.Pe, 'on');
grid(ax.Pe, 'on');
box(ax.Pe, 'on');

plot(ax.Pe, tout, Perr_nonlin, ...
    'b-', ...
    'LineWidth', 1.4);

plot(ax.Pe, tout, Perr_lin, ...
    'k:', ...
    'LineWidth', 1.4);

ylabel(ax.Pe, 'Relative error (\%)', ...
    'FontSize', 10, ...
    'Interpreter', 'latex');

xlabel(ax.Pe, 'Time (s)', ...
    'FontSize', 10, ...
    'Interpreter', 'latex');

legend(ax.Pe, ...
    'vs. $\mathrm{FOM}_{\mathrm{NL-NL}}$', ...
    'vs. $\mathrm{FOM}_{\mathrm{Lin-Lin}}$', ...
    'FontSize', 8, ...
    'Interpreter', 'latex', ...
    'Location', 'best');

set(ax.Pe, ...
    'FontSize', 9, ...
    'Layer', 'top');

set(ax.Pe, 'XTickLabel', []);

hold(ax.Pe, 'off');


% ============================================================
% Add xlabel to the middle MAIN axes WITHOUT changing tile gap
% ============================================================

% Do NOT use xlabel(ax.P,...) here if MATLAB changes your tiled spacing.
% Instead use text, which does not cause tiledlayout to resize the axes.

text(ax.P, ...
    0.5, -0.18, ...
    'Time (s)', ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16, ...
    'Interpreter', 'latex');


% ============================================================
% Final Formatting
% ============================================================

% Link only main figures.
% Do NOT include inset, otherwise zooming can become inconvenient.
linkaxes([ax.Pin, ax.P], 'x');

% Keep original 3-panel figure dimensions
set(fig, ...
    'Units', 'inches', ...
    'Position', [1, 1, 7.16, 6.5]);