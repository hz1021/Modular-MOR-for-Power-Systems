format short g;
w = 2*pi*50;
Imax = 15;

%%%%% sigmaplot settings %%%%%%%%%%%%
fs = 14;
p = sigmaoptions('cstprefs');
p.FreqUnits = 'Hz';
p.Grid = 'on';
p.Title.FontSize = 12;
p.YLabel.FontSize = fs;
p.XLabel.FontSize = fs;
wmin = 0.1*2*pi; wmax = 1e4*2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% bodeplot settings %%%%%%%%%%%%
opts = bodeoptions;
opts.PhaseVisible = 'off';
opts.FreqUnits = 'Hz';
% opts.Xlim = {[1 100]};
% opts.Ylim = {[-10 10]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Simulation settings %%%%%%%%%%%%%%%%%
Ts = 50e-6; % samping time. Should be entered in all Simulink files.
Tstop = 1;
%%%%%%%%%%%%%%%%%%%%%%%%