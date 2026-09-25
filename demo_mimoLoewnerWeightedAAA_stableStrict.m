%DEMO_MIMOLOEWNERWEIGHTEDAAA_STABLESTRICT
% Stable strict weighted adaptive MIMO Loewner/AAA demo.
%
% This demo assumes your workspace already contains:
%   sys_j, H_j, WTD, VTD, w, ic
% as in your original demo_mimoLoewnerWeightedAAA.m.

% A stable subsystem.
A = sys_j(ic).s.A; %#ok<NASGU>
B = sys_j(ic).s.B; %#ok<NASGU>
C = sys_j(ic).s.C; %#ok<NASGU>
D = sys_j(ic).s.D;

Gpos = H_j{ic};
% Gpos = H_2_update;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation-direction weights.
% Leave these empty for unweighted Loewner tangential data.
WoPos = [];
WiPos = [];

% 2. Exact weight
WoPos = zeros(size(D, 1), size(D, 1), numel(w));
WiPos = zeros(size(D, 2), size(D, 2), numel(w));
for k = 1:numel(w)
    WoPos(:,:,k) = WTD.inv{ic, k};
    WiPos(:,:,k) = VTD.inv{ic, k};
end

% If you want to use fitted interpolation-direction weights instead, uncomment:
% WoPos = zeros(size(D,1),size(D,1),numel(w));
% WiPos = zeros(size(D,2),size(D,2),numel(w));
% for k = 1:numel(w)
%     WoPos(:,:,k) = freqresp(Gy{ic},w(k));
%     WiPos(:,:,k) = freqresp(Gu{ic},w(k));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard frequency-weighted error-bound weights.
WoErrBound = zeros(size(D,1),size(D,1),numel(w));
WiErrBound = zeros(size(D,2),size(D,2),numel(w));
for k = 1:numel(w)
    WoErrBound(:,:,k) = WTD.inv{ic,k};
    WiErrBound(:,:,k) = VTD.inv{ic,k};
end

% Mirror data for real realization.
[s,G,Wo,Wi,WoErrBoundFull,WiErrBoundFull] = mirrorMimoWeightedFrequencyData( ...
    w,Gpos,WoPos,WiPos,WoErrBound,WiErrBound);

opts = struct();
opts.D = D;
opts.specTol = 1;                  % strict post-stabilization bound
opts.maxSupportPairs = 1000;
opts.maxOrder = 200;               % increase if no stable feasible ROM is found
opts.order = [];                   % auto order
opts.svtol = 1e-11;
opts.normType = 'sigma';
opts.forceReal = true;
opts.verbose = true;
opts.makeDSS = true;
opts.sideSelection = 'balanced';   % balanced or alternate

% Strict stable/feasible return.
opts.enforceStableROM = true;
opts.rawFeasibleStabilizationRatio = 1;     % attempt stabilization only after raw ROM satisfies bound
opts.postStabilizationWeightedRatio = 1;    % accept only if stabilized ROM also satisfies bound
opts.requireStableFeasibleReturn = true;

% Use post-stabilized residual for next greedy point whenever stabilization
% was attempted. This helps correct the error introduced by stabilization.
opts.greedyModelSource = 'postStabilizedWhenAttempted';

% Stability enforcement and accuracy recovery.
opts.stabilityTol = 1e-10;
opts.EcondTol = 1e-12;
opts.regularizeSingularE = true;
opts.EsvdTol = 1e-12;
opts.stableUpdateMatrix = 'C';              % refit C after pole stabilization
opts.refitAfterStabilization = true;
opts.stableUpdateFitWeightMode = 'scalarWeighted';
opts.lsRidge = 1e-12;

% Loosely reduce unweighted error after a stable feasible ROM is found.
opts.weightedGuard = 1;
opts.unweightedMixWeight = 0;
opts.nearFeasibleDirectionMode = 'weighted';
opts.feasibleDirectionMode = 'weighted';
opts.unweightedMetric = 'relative';
opts.continueAfterStableFeasible = true;
opts.extraBlocksAfterStableFeasible = 0;
opts.unweightedPatience = 10;
opts.unweightedImproveTol = 1e-3;

tic;
[model,info] = mimoLoewnerWeightedAAA_stableStrict( ...
    s,G,Wo,Wi,opts,WoErrBoundFull,WiErrBoundFull);
toc;

GfitPos = evalMimoLoewnerModel(model,1i*w);

weightedErr = zeros(numel(w),1);
unweightedErr = zeros(numel(w),1);
unweightedRelErr = zeros(numel(w),1);
for k = 1:numel(w)
    Rk = GfitPos(:,:,k) - Gpos(:,:,k);
    weightedErr(k) = norm(WoErrBound(:,:,k)*Rk*WiErrBound(:,:,k),2);
    unweightedErr(k) = norm(Rk,2);
    unweightedRelErr(k) = unweightedErr(k)/max(1,norm(Gpos(:,:,k),2));
end

fprintf('\nReturned ROM order: %d\n',model.order);
fprintf('Returned ROM stable: %d\n',model.stableOnReturn);
fprintf('Max positive-frequency weighted error: %.6e\n',max(weightedErr));
fprintf('Grid specification satisfied: %d\n',max(weightedErr) <= opts.specTol);
fprintf('Max positive-frequency unweighted abs error: %.6e\n',max(unweightedErr));
fprintf('Max positive-frequency unweighted rel error: %.6e\n',max(unweightedRelErr));
fprintf('Max real pole: %.6e\n',max(real(eig(model.A,model.E))));

if opts.makeDSS && ~isempty(model.sys)
    fprintf('isstable(model.sys): %d\n',isstable(model.sys));
end

figure;
semilogx(w,weightedErr,'LineWidth',1.2); grid on;
yline(opts.specTol,'--');
xlabel('Frequency [rad/s]');
ylabel('||Wo(jw)(Gr(jw)-G(jw))Wi(jw)||_2');
title('Strict post-stabilization weighted error');

figure;
semilogx(w,unweightedErr,'LineWidth',1.2); grid on;
xlabel('Frequency [rad/s]');
ylabel('||Gr(jw)-G(jw)||_2');
title('Unweighted absolute error');

figure;
semilogx(w,unweightedRelErr,'LineWidth',1.2); grid on;
xlabel('Frequency [rad/s]');
ylabel('||Gr-G||_2 / max(1,||G||_2)');
title('Unweighted relative error');
