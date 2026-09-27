function [model,info] = mimoLoewnerWeightedAAA_stableStrict(s,G,Wo,Wi,opts,WoErrBound,WiErrBound)
%MIMOLOEWNERWEIGHTEDAAA_STABLESTRICT
% Stable strict frequency-weighted adaptive MIMO Loewner/AAA reduction.
%
% This is a drop-in variant of mimoLoewnerWeightedAAA with stricter return
% logic:
%
%   1) The greedy Loewner/AAA search is still driven by grid residuals.
%   2) Stabilization is attempted only when the RAW Loewner ROM is already
%      frequency-weighted feasible, unless rawFeasibleStabilizationRatio > 1.
%   3) After stabilization/refitting, the weighted grid bound is checked
%      again.
%   4) The function returns only a STABLE ROM that satisfies the post-
%      stabilization frequency-weighted bound. It never returns a merely
%      stable-but-infeasible ROM.
%   5) Among stable weighted-feasible candidates, the returned ROM is the one
%      with the lowest relative unweighted grid error.
%
% Syntax:
%   [model,info] = mimoLoewnerWeightedAAA_stableStrict(s,G,Wo,Wi,opts)
%   [model,info] = mimoLoewnerWeightedAAA_stableStrict(s,G,Wo,Wi,opts,WoErrBound,WiErrBound)
%
% Transfer function:
%   Gr(s) = C*((s*E - A)\B) + D
%
% Required data:
%   s     N-by-1 complex samples, conjugate symmetric for forceReal=true.
%   G     p-by-m-by-N samples.
%   Wo,Wi interpolation-direction weights; [] means identity.
%
% Optional hard error weights:
%   WoErrBound, WiErrBound define the actual frequency-weighted error bound:
%       || WoErrBound(:,:,k)*(G-Gr)*WiErrBound(:,:,k) || <= specTol(k)
%   If omitted, Wo and Wi are used for the hard bound.
%
% Important options:
%   opts.specTol = 1 or N-by-1 vector
%   opts.maxSupportPairs = 200
%   opts.maxOrder = Inf
%   opts.svtol = 1e-10
%   opts.sideSelection = 'balanced' or 'alternate'
%
% Strict stability/feasibility options:
%   opts.enforceStableROM = true
%   opts.rawFeasibleStabilizationRatio = 1
%   opts.postStabilizationWeightedRatio = 1
%   opts.requireStableFeasibleReturn = true
%
% Accuracy-shaping options:
%   opts.weightedGuard = 1.2
%   opts.unweightedMixWeight = 0.2
%   opts.nearFeasibleDirectionMode = 'weighted'
%   opts.feasibleDirectionMode = 'unweighted'
%   opts.continueAfterStableFeasible = true
%   opts.extraBlocksAfterStableFeasible = 40
%   opts.unweightedPatience = 10
%
% Stabilization/refit:
%   opts.stableUpdateMatrix = 'C'    % currently implemented: 'C' or 'none'
%   opts.fitD = false
%   opts.refitAfterStabilization = true
%
% Notes:
%   This is a grid-based algorithm. Strict satisfaction is checked only on
%   the supplied frequency grid.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    s = s(:);
    [p,m,N] = size(G);
    if numel(s) ~= N
        error('numel(s) must equal size(G,3).');
    end
    if N < 2
        error('At least two samples are required.');
    end

    opts = localDefaults(opts,p,m,N);
    D = localFeedthrough(opts.D,p,m);
    if opts.forceReal
        D = real(D);
    end

    % Direction/interpolation weights.
    Wo = localExpandWeight(Wo,N,p,'Wo');
    Wi = localExpandInputWeight(Wi,N,m,'Wi');

    % Hard error-bound weights. If omitted, use Wo/Wi.
    if nargin < 6 || isempty(WoErrBound)
        WoErrBound = Wo;
    else
        WoErrBound = localExpandWeight(WoErrBound,N,p,'WoErrBound');
    end
    if nargin < 7 || isempty(WiErrBound)
        WiErrBound = Wi;
    else
        WiErrBound = localExpandInputWeight(WiErrBound,N,m,'WiErrBound');
    end

    pairTolAbs = opts.pairTol * max(1,max(abs(s)));
    partners = localConjPartners(s,pairTolAbs);
    if opts.forceReal && any(isnan(partners))
        error(['forceReal=true requires conjugate-symmetric samples. ', ...
               'Use mirrorMimoWeightedFrequencyData first.']);
    end

    [reps,repOfIndex] = localPairRepresentatives(s,partners);

    if opts.forceReal && opts.symmetrizeData
        [G,Wo,Wi,infoSym,WoErrBound,WiErrBound] = localSymmetrizeConjugateData( ...
            s,G,Wo,Wi,partners,pairTolAbs,WoErrBound,WiErrBound);
    else
        infoSym = struct('relG',NaN,'relWo',NaN,'relWi',NaN,'relWoErrBound',NaN,'relWiErrBound',NaN);
    end

    % Initial seeds: use absolute hard weighted error, not relative error.
    baseScores = zeros(numel(reps),1);
    for a = 1:numel(reps)
        k = reps(a);
        Rk = D - G(:,:,k);
        baseScores(a) = localMatNorm(WoErrBound(:,:,k)*Rk*WiErrBound(:,:,k),opts.normType);
    end

    [~,aRight] = max(baseScores);
    if numel(reps) < 2
        error('Need at least two distinct conjugate support blocks.');
    end

    sep = abs(s(reps) - s(reps(aRight)));
    if max(sep) > 0
        secondScore = baseScores .* (sep/max(sep));
    else
        secondScore = baseScores;
    end
    secondScore(aRight) = -Inf;
    [~,aLeft] = max(secondScore);
    if ~isfinite(secondScore(aLeft))
        error('Could not find two distinct initial support blocks.');
    end

    leftIdx = [];
    rightIdx = [];
    Ldir = zeros(p,0);
    Rdir = zeros(m,0);
    sideOfRep = zeros(numel(reps),1);

    % Seed right.
    kR = reps(aRight);
    RR = D - G(:,:,kR);
    [~,rR] = localWeightedDirections(RR,WoErrBound(:,:,kR),WiErrBound(:,:,kR));
    [rightIdx,Rdir] = localAddRightBlock(rightIdx,Rdir,kR,rR,s,partners,pairTolAbs,opts.forceReal);
    sideOfRep(aRight) = 2;

    % Seed left.
    kL = reps(aLeft);
    RL = D - G(:,:,kL);
    [ellL,~] = localWeightedDirections(RL,WoErrBound(:,:,kL),WiErrBound(:,:,kL));
    [leftIdx,Ldir] = localAddLeftBlock(leftIdx,Ldir,kL,ellL,s,partners,pairTolAbs,opts.forceReal);
    sideOfRep(aLeft) = 1;

    nextSide = 2;
    nBlocks = 2;

    info = struct();
    info.weightedRatioHistRaw = [];
    info.weightedRatioHistPost = [];
    info.weightedErrHistRaw = [];
    info.unweightedErrHistRaw = [];
    info.unweightedErrHistPost = [];
    info.orderHist = [];
    info.nLeftHist = [];
    info.nRightHist = [];
    info.selectedRepHist = [aRight; aLeft];
    info.selectedSideHist = [2; 1];
    info.selectedModeHist = {};
    info.greedySourceHist = {};
    info.postAcceptedHist = [];
    info.stabilizationAttemptHist = [];
    info.conjugateSymmetryErrorBeforeSymmetrization = infoSym;

    bestModel = [];
    bestMetrics = [];
    bestInfo = [];
    bestUnweighted = Inf;
    stableFeasibleFound = false;
    extraAfterStableFeasible = 0;
    unweightedNoImprove = 0;

    rawBuildOpts = opts;
    rawBuildOpts.stabilize = false;   % Never stabilize inside Loewner build.

    while true
        [rawModel,buildInfo] = localBuildReducedModelRaw( ...
            s,G,D,leftIdx,rightIdx,Ldir,Rdir,partners,pairTolAbs,rawBuildOpts);

        rawMetrics = localGridMetrics(rawModel,s,G,WoErrBound,WiErrBound,reps,partners,opts);

        postModel = [];
        postMetrics = [];
        postInfo = struct();
        postAttempted = false;
        postAccepted = false;

        % Attempt stabilization only when raw model is already feasible enough.
        if opts.enforceStableROM && ...
                rawMetrics.maxWeightedRatio <= opts.rawFeasibleStabilizationRatio

            postAttempted = true;
            [postModel,postInfo] = localStabilizeAndRefit( ...
                rawModel,s,G,WoErrBound,WiErrBound,opts);

            postMetrics = localGridMetrics(postModel,s,G,WoErrBound,WiErrBound,reps,partners,opts);

            postStable = localIsStable(postModel,opts);
            postAccepted = postStable && ...
                postMetrics.maxWeightedRatio <= opts.postStabilizationWeightedRatio;

            if postAccepted
                if ~stableFeasibleFound
                    stableFeasibleFound = true;
                    extraAfterStableFeasible = 0;
                    unweightedNoImprove = 0;
                end

                if postMetrics.maxUnweightedRel < bestUnweighted*(1 - opts.unweightedImproveTol)
                    bestModel = postModel;
                    bestMetrics = postMetrics;
                    bestInfo = postInfo;
                    bestUnweighted = postMetrics.maxUnweightedRel;
                    unweightedNoImprove = 0;
                else
                    unweightedNoImprove = unweightedNoImprove + 1;
                end
            end
        end

        info.weightedRatioHistRaw(end+1,1) = rawMetrics.maxWeightedRatio; %#ok<AGROW>
        info.weightedErrHistRaw(end+1,1) = rawMetrics.maxWeightedErr; %#ok<AGROW>
        info.unweightedErrHistRaw(end+1,1) = rawMetrics.maxUnweightedRel; %#ok<AGROW>
        info.orderHist(end+1,1) = rawModel.order; %#ok<AGROW>
        info.nLeftHist(end+1,1) = numel(leftIdx); %#ok<AGROW>
        info.nRightHist(end+1,1) = numel(rightIdx); %#ok<AGROW>
        info.stabilizationAttemptHist(end+1,1) = postAttempted; %#ok<AGROW>
        info.postAcceptedHist(end+1,1) = postAccepted; %#ok<AGROW>

        if postAttempted
            info.weightedRatioHistPost(end+1,1) = postMetrics.maxWeightedRatio; %#ok<AGROW>
            info.unweightedErrHistPost(end+1,1) = postMetrics.maxUnweightedRel; %#ok<AGROW>
        else
            info.weightedRatioHistPost(end+1,1) = NaN; %#ok<AGROW>
            info.unweightedErrHistPost(end+1,1) = NaN; %#ok<AGROW>
        end

        if opts.verbose
            if postAttempted
                fprintf(['blocks=%4d nL=%4d nR=%4d order=%4d ', ...
                         'rawRatio=%.4e rawUnw=%.4e postRatio=%.4e postUnw=%.4e postOK=%d\n'], ...
                    nBlocks,numel(leftIdx),numel(rightIdx),rawModel.order, ...
                    rawMetrics.maxWeightedRatio,rawMetrics.maxUnweightedRel, ...
                    postMetrics.maxWeightedRatio,postMetrics.maxUnweightedRel,postAccepted);
            else
                fprintf(['blocks=%4d nL=%4d nR=%4d order=%4d ', ...
                         'rawRatio=%.4e rawUnw=%.4e stab=0\n'], ...
                    nBlocks,numel(leftIdx),numel(rightIdx),rawModel.order, ...
                    rawMetrics.maxWeightedRatio,rawMetrics.maxUnweightedRel);
            end
        end

        % Stop conditions after a stable feasible model has been found.
        if stableFeasibleFound
            if ~opts.continueAfterStableFeasible
                break;
            end
            if extraAfterStableFeasible >= opts.extraBlocksAfterStableFeasible
                break;
            end
            if ~opts.exactExtraBlocksAfterStableFeasible && ...
                    unweightedNoImprove >= opts.unweightedPatience
                break;
            end
        end

        if nBlocks >= opts.maxSupportPairs
            break;
        end

        % Choose which model residual drives the next interpolation point.
        if postAttempted && strcmpi(opts.greedyModelSource,'postStabilizedWhenAttempted')
            greedyModel = postModel;
            greedyMetrics = postMetrics;
            greedySource = 'post';
        elseif postAttempted && strcmpi(opts.greedyModelSource,'postStabilizedIfAccepted') && postAccepted
            greedyModel = postModel;
            greedyMetrics = postMetrics;
            greedySource = 'post';
        else
            greedyModel = rawModel;
            greedyMetrics = rawMetrics;
            greedySource = 'raw';
        end

        [aNew,side,desiredSide,selectMode] = localChooseNextBlock( ...
            sideOfRep,errScoreForSelection(greedyMetrics),greedyMetrics, ...
            leftIdx,rightIdx,nextSide,opts);

        if isempty(aNew)
            if opts.verbose
                fprintf('No inactive conjugate support block remains. Stopping search.\n');
            end
            break;
        end

        kNew = reps(aNew);
        Rnew = localEvalOne(greedyModel,s(kNew)) - G(:,:,kNew);

        directionMode = localDirectionMode(greedyMetrics,selectMode,opts);
        if strcmpi(directionMode,'weighted')
            [ellNew,rNew] = localWeightedDirections(Rnew,WoErrBound(:,:,kNew),WiErrBound(:,:,kNew));
        elseif strcmpi(directionMode,'unweighted')
            [ellNew,rNew] = localUnweightedDirections(Rnew);
        else
            error('Unknown directionMode: %s',directionMode);
        end

        if side == 1
            [leftIdx,Ldir] = localAddLeftBlock(leftIdx,Ldir,kNew,ellNew,s,partners,pairTolAbs,opts.forceReal);
            sideOfRep(aNew) = 1;
        else
            [rightIdx,Rdir] = localAddRightBlock(rightIdx,Rdir,kNew,rNew,s,partners,pairTolAbs,opts.forceReal);
            sideOfRep(aNew) = 2;
        end

        nBlocks = nBlocks + 1;
        if stableFeasibleFound
            extraAfterStableFeasible = extraAfterStableFeasible + 1;
        end

        info.selectedRepHist(end+1,1) = aNew; %#ok<AGROW>
        info.selectedSideHist(end+1,1) = side; %#ok<AGROW>
        info.selectedModeHist{end+1,1} = sprintf('%s/%s',selectMode,directionMode); %#ok<AGROW>
        info.greedySourceHist{end+1,1} = greedySource; %#ok<AGROW>

        if strcmpi(opts.sideSelection,'alternate') || ...
                strcmpi(opts.sideSelection,'alternating') || ...
                strcmpi(opts.sideSelection,'strictalternate')
            nextSide = 3 - nextSide;
        end
    end

    if isempty(bestModel)
        msg = ['No stabilized ROM satisfying the frequency-weighted bound was found. ', ...
               'Try increasing opts.maxSupportPairs, opts.maxOrder, lowering opts.svtol, ', ...
               'or relaxing opts.postStabilizationWeightedRatio.'];
        if opts.requireStableFeasibleReturn
            error(msg);
        else
            warning(msg);
            model = rawModel;
            finalMetrics = rawMetrics;
        end
    else
        model = bestModel;
        finalMetrics = bestMetrics;
    end

    model.maxWeightedErr = finalMetrics.maxWeightedErr;
    model.maxWeightedRatio = finalMetrics.maxWeightedRatio;
    model.maxUnweightedRel = finalMetrics.maxUnweightedRel;
    model.maxUnweightedAbs = finalMetrics.maxUnweightedAbs;
    model.specTol = opts.specTol;
    model.weightedSpecSatisfiedOnGrid = finalMetrics.maxWeightedRatio <= opts.postStabilizationWeightedRatio;
    model.stableOnReturn = localIsStable(model,opts);

    info.finalMetrics = finalMetrics;
    info.bestStabilizationInfo = bestInfo;
    info.maxWeightedErr = finalMetrics.maxWeightedErr;
    info.maxWeightedRatio = finalMetrics.maxWeightedRatio;
    info.maxUnweightedRel = finalMetrics.maxUnweightedRel;
    info.maxUnweightedAbs = finalMetrics.maxUnweightedAbs;
    info.errAllSamples = finalMetrics.weightedErrAll;
    info.weightedRatioAllSamples = finalMetrics.weightedRatioAll;
    info.unweightedRelAllSamples = finalMetrics.unweightedRelAll;
    info.selectedBlocks = nBlocks;
    info.leftIdx = leftIdx;
    info.rightIdx = rightIdx;
    info.leftPoints = s(leftIdx);
    info.rightPoints = s(rightIdx);
    info.leftDirections = Ldir;
    info.rightDirections = Rdir;
    info.representativeIndices = reps;
    info.repOfIndex = repOfIndex;
    info.sideOfRep = sideOfRep;
    info.singularValuesHorizontal = buildInfo.svHorizontal;
    info.singularValuesVertical = buildInfo.svVertical;
    info.realificationImaginaryLeakage = buildInfo.imagLeak;

    model.leftPoints = info.leftPoints;
    model.rightPoints = info.rightPoints;
    model.leftDirections = Ldir;
    model.rightDirections = Rdir;

    if opts.makeDSS
        try
            model.sys = dss(model.A,model.B,model.C,model.D,model.E);
        catch
            model.sys = [];
        end
    end
end

function opts = localDefaults(opts,p,m,N)
    opts = localSetDefault(opts,'specTol',1);
    opts = localSetDefault(opts,'maxSupportPairs',200);
    opts = localSetDefault(opts,'maxOrder',Inf);
    opts = localSetDefault(opts,'order',[]);
    opts = localSetDefault(opts,'svtol',1e-10);
    opts = localSetDefault(opts,'D',zeros(p,m));
    opts = localSetDefault(opts,'normType','sigma');
    opts = localSetDefault(opts,'forceReal',true);
    opts = localSetDefault(opts,'symmetrizeData',true);
    opts = localSetDefault(opts,'pairTol',1e-10);
    opts = localSetDefault(opts,'verbose',false);
    opts = localSetDefault(opts,'makeDSS',true);
    opts = localSetDefault(opts,'sideSelection','balanced');

    % Strict stable/feasible return behavior.
    opts = localSetDefault(opts,'enforceStableROM',true);
    opts = localSetDefault(opts,'rawFeasibleStabilizationRatio',1);
    opts = localSetDefault(opts,'postStabilizationWeightedRatio',1);
    opts = localSetDefault(opts,'requireStableFeasibleReturn',true);
    opts = localSetDefault(opts,'greedyModelSource','postStabilizedWhenAttempted');

    % Stability enforcement.
    opts = localSetDefault(opts,'stabilityTol',1e-10);
    opts = localSetDefault(opts,'EcondTol',1e-12);
    opts = localSetDefault(opts,'regularizeSingularE',true);
    opts = localSetDefault(opts,'EsvdTol',1e-12);
    opts = localSetDefault(opts,'stableUpdateMatrix','C');
    opts = localSetDefault(opts,'fitD',false);
    opts = localSetDefault(opts,'refitAfterStabilization',true);
    opts = localSetDefault(opts,'lsRidge',1e-12);
    opts = localSetDefault(opts,'stableUpdateFitWeightMode','scalarWeighted');

    % Weighted/unweighted search policy.
    opts = localSetDefault(opts,'weightedGuard',1.2);
    opts = localSetDefault(opts,'unweightedMixWeight',0.2);
    opts = localSetDefault(opts,'unweightedMetric','relative');
    opts = localSetDefault(opts,'unweightedFloor',1e-12);
    opts = localSetDefault(opts,'nearFeasibleDirectionMode','weighted');
    opts = localSetDefault(opts,'feasibleDirectionMode','unweighted');

    % Continue after first stable weighted-feasible candidate to reduce
    % unweighted error.
    opts = localSetDefault(opts,'continueAfterStableFeasible',true);
    opts = localSetDefault(opts,'extraBlocksAfterStableFeasible',40);
    opts = localSetDefault(opts,'exactExtraBlocksAfterStableFeasible',false);
    opts = localSetDefault(opts,'unweightedPatience',10);
    opts = localSetDefault(opts,'unweightedImproveTol',1e-3);

    if ~isscalar(opts.specTol)
        if numel(opts.specTol) ~= N
            error('opts.specTol must be scalar or length N.');
        end
        opts.specTol = opts.specTol(:);
    end
end

function opts = localSetDefault(opts,name,value)
    if ~isfield(opts,name) || isempty(opts.(name))
        opts.(name) = value;
    end
end

function D = localFeedthrough(D,p,m)
    if isempty(D)
        D = zeros(p,m);
    elseif isscalar(D)
        D = D*ones(p,m);
    end
    if ~isequal(size(D),[p m])
        error('opts.D must be scalar or p-by-m.');
    end
end

function W = localExpandWeight(W,N,p,name)
    if isempty(W)
        W = repmat(eye(p),1,1,N);
        return;
    end
    if ndims(W) == 2
        if size(W,2) ~= p
            error('%s must have p columns.',name);
        end
        W = repmat(W,1,1,N);
    elseif ndims(W) == 3
        if size(W,2) ~= p || size(W,3) ~= N
            error('%s must be qo-by-p-by-N.',name);
        end
    else
        error('%s must be a matrix, a 3-D array, or [].',name);
    end
end

function W = localExpandInputWeight(W,N,m,name)
    if isempty(W)
        W = repmat(eye(m),1,1,N);
        return;
    end
    if ndims(W) == 2
        if size(W,1) ~= m
            error('%s must have m rows.',name);
        end
        W = repmat(W,1,1,N);
    elseif ndims(W) == 3
        if size(W,1) ~= m || size(W,3) ~= N
            error('%s must be m-by-ri-by-N.',name);
        end
    else
        error('%s must be a matrix, a 3-D array, or [].',name);
    end
end

function partners = localConjPartners(s,tolAbs)
    N = numel(s);
    partners = NaN(N,1);
    for k = 1:N
        [d,j] = min(abs(s - conj(s(k))));
        if d <= tolAbs
            partners(k) = j;
        end
    end
end

function [reps,repOfIndex] = localPairRepresentatives(s,partners)
    N = numel(s);
    seen = false(N,1);
    reps = [];
    repOfIndex = zeros(N,1);
    for k = 1:N
        if seen(k), continue; end
        kp = partners(k);
        if isnan(kp) || kp == k
            rep = k;
            pair = k;
        else
            pair = unique([k; kp]);
            if imag(s(k)) > imag(s(kp))
                rep = k;
            elseif imag(s(kp)) > imag(s(k))
                rep = kp;
            else
                rep = min(k,kp);
            end
        end
        reps(end+1,1) = rep; %#ok<AGROW>
        id = numel(reps);
        repOfIndex(pair) = id;
        seen(pair) = true;
    end
end

function [G,Wo,Wi,info,WoErrBound,WiErrBound] = localSymmetrizeConjugateData( ...
    s,G,Wo,Wi,partners,tolAbs,WoErrBound,WiErrBound)
    N = numel(s);
    seen = false(N,1);
    gNum=0; gDen=0; woNum=0; woDen=0; wiNum=0; wiDen=0; webNum=0; webDen=0; wibNum=0; wibDen=0;

    for k = 1:N
        if seen(k), continue; end
        kp = partners(k);
        if isnan(kp) || kp == k || abs(s(k)-conj(s(k))) <= tolAbs
            gDen = gDen + norm(G(:,:,k),'fro');
            woDen = woDen + norm(Wo(:,:,k),'fro');
            wiDen = wiDen + norm(Wi(:,:,k),'fro');
            webDen = webDen + norm(WoErrBound(:,:,k),'fro');
            wibDen = wibDen + norm(WiErrBound(:,:,k),'fro');

            gNum = gNum + norm(imag(G(:,:,k)),'fro');
            woNum = woNum + norm(imag(Wo(:,:,k)),'fro');
            wiNum = wiNum + norm(imag(Wi(:,:,k)),'fro');
            webNum = webNum + norm(imag(WoErrBound(:,:,k)),'fro');
            wibNum = wibNum + norm(imag(WiErrBound(:,:,k)),'fro');

            G(:,:,k) = real(G(:,:,k));
            Wo(:,:,k) = real(Wo(:,:,k));
            Wi(:,:,k) = real(Wi(:,:,k));
            WoErrBound(:,:,k) = real(WoErrBound(:,:,k));
            WiErrBound(:,:,k) = real(WiErrBound(:,:,k));
            seen(k) = true;
        else
            Gavg = 0.5*(G(:,:,k)+conj(G(:,:,kp)));
            Woavg = 0.5*(Wo(:,:,k)+conj(Wo(:,:,kp)));
            Wiavg = 0.5*(Wi(:,:,k)+conj(Wi(:,:,kp)));
            Webavg = 0.5*(WoErrBound(:,:,k)+conj(WoErrBound(:,:,kp)));
            Wibavg = 0.5*(WiErrBound(:,:,k)+conj(WiErrBound(:,:,kp)));

            gNum = gNum + norm(G(:,:,k)-conj(G(:,:,kp)),'fro');
            woNum = woNum + norm(Wo(:,:,k)-conj(Wo(:,:,kp)),'fro');
            wiNum = wiNum + norm(Wi(:,:,k)-conj(Wi(:,:,kp)),'fro');
            webNum = webNum + norm(WoErrBound(:,:,k)-conj(WoErrBound(:,:,kp)),'fro');
            wibNum = wibNum + norm(WiErrBound(:,:,k)-conj(WiErrBound(:,:,kp)),'fro');

            gDen = gDen + norm(G(:,:,k),'fro') + norm(G(:,:,kp),'fro');
            woDen = woDen + norm(Wo(:,:,k),'fro') + norm(Wo(:,:,kp),'fro');
            wiDen = wiDen + norm(Wi(:,:,k),'fro') + norm(Wi(:,:,kp),'fro');
            webDen = webDen + norm(WoErrBound(:,:,k),'fro') + norm(WoErrBound(:,:,kp),'fro');
            wibDen = wibDen + norm(WiErrBound(:,:,k),'fro') + norm(WiErrBound(:,:,kp),'fro');

            G(:,:,k) = Gavg; G(:,:,kp) = conj(Gavg);
            Wo(:,:,k) = Woavg; Wo(:,:,kp) = conj(Woavg);
            Wi(:,:,k) = Wiavg; Wi(:,:,kp) = conj(Wiavg);
            WoErrBound(:,:,k) = Webavg; WoErrBound(:,:,kp) = conj(Webavg);
            WiErrBound(:,:,k) = Wibavg; WiErrBound(:,:,kp) = conj(Wibavg);

            seen([k kp]) = true;
        end
    end

    info.relG = gNum/max(1,gDen);
    info.relWo = woNum/max(1,woDen);
    info.relWi = wiNum/max(1,wiDen);
    info.relWoErrBound = webNum/max(1,webDen);
    info.relWiErrBound = wibNum/max(1,wibDen);
end

function [leftIdx,Ldir] = localAddLeftBlock(leftIdx,Ldir,k,ell,s,partners,tolAbs,forceReal)
    [idxBlock,dirBlock] = localSupportBlock(k,ell,s,partners,tolAbs,forceReal);
    leftIdx = [leftIdx; idxBlock(:)];
    Ldir = [Ldir, dirBlock];
end

function [rightIdx,Rdir] = localAddRightBlock(rightIdx,Rdir,k,r,s,partners,tolAbs,forceReal)
    [idxBlock,dirBlock] = localSupportBlock(k,r,s,partners,tolAbs,forceReal);
    rightIdx = [rightIdx; idxBlock(:)];
    Rdir = [Rdir, dirBlock];
end

function [idxBlock,dirBlock] = localSupportBlock(k,d,s,partners,tolAbs,forceReal)
    d = localNormalizeDirection(d);
    kp = partners(k);
    if isnan(kp) || kp == k || abs(s(k)-conj(s(k))) <= tolAbs
        if forceReal
            d = real(d);
            d = localNormalizeDirection(d);
        end
        idxBlock = k;
        dirBlock = d;
    else
        if imag(s(k)) < imag(s(kp))
            tmp = k; k = kp; kp = tmp;
            d = conj(d);
        end
        idxBlock = [k; kp];
        dirBlock = [d, conj(d)];
    end
end

function d = localNormalizeDirection(d)
    nd = norm(d);
    if nd < 100*eps
        d = zeros(size(d));
        d(1) = 1;
    else
        d = d/nd;
    end
end

function [ell,r,sigma1] = localWeightedDirections(R,Wo,Wi)
    M = Wo*R*Wi;
    if isempty(M) || norm(M,'fro') < 100*eps
        [ell,r] = localUnweightedDirections(R);
        sigma1 = 0;
        return;
    end
    [U,S,V] = svd(M,'econ');
    sigma1 = S(1,1);
    ell = Wo' * U(:,1);
    r = Wi * V(:,1);
    if norm(ell) < 100*eps || norm(r) < 100*eps
        [ell,r] = localUnweightedDirections(R);
    else
        ell = ell/norm(ell);
        r = r/norm(r);
    end
end

function [ell,r] = localUnweightedDirections(R)
    if norm(R,'fro') < 100*eps
        ell = zeros(size(R,1),1); ell(1) = 1;
        r = zeros(size(R,2),1); r(1) = 1;
        return;
    end
    [U,~,V] = svd(R,'econ');
    ell = U(:,1);
    r = V(:,1);
end

function [model,info] = localBuildReducedModelRaw(s,G,D,leftIdx,rightIdx,Ldir,Rdir,partners,tolAbs,opts)
    [p,m,~] = size(G);
    nL = numel(leftIdx);
    nR = numel(rightIdx);
    if nL == 0 || nR == 0
        error('Both left and right interpolation sets must be nonempty.');
    end

    L = zeros(nL,nR);
    Ls = zeros(nL,nR);
    B = zeros(nL,m);
    C = zeros(p,nR);

    for j = 1:nR
        lamIdx = rightIdx(j);
        rj = Rdir(:,j);
        Glam = G(:,:,lamIdx) - D;
        C(:,j) = Glam*rj;
    end

    for i = 1:nL
        muIdx = leftIdx(i);
        mu = s(muIdx);
        ell = Ldir(:,i);
        Gmu = G(:,:,muIdx) - D;
        vi = ell' * Gmu;
        B(i,:) = vi;
        for j = 1:nR
            lamIdx = rightIdx(j);
            lam = s(lamIdx);
            rj = Rdir(:,j);
            den = mu - lam;
            if abs(den) <= 10*tolAbs
                error('Left and right support sets contain same/nearly same point.');
            end
            alpha = vi*rj;
            beta = ell'*C(:,j);
            L(i,j)  = (alpha - beta)/den;
            Ls(i,j) = (mu*alpha - lam*beta)/den;
        end
    end

    E0 = -L; A0 = -Ls; B0 = B; C0 = C; D0 = D;
    imagLeak = NaN;

    if opts.forceReal
        QL = localRealifyingTransform(leftIdx,s,partners,tolAbs);
        QR = localRealifyingTransform(rightIdx,s,partners,tolAbs);
        E0 = QL*E0*QR';
        A0 = QL*A0*QR';
        B0 = QL*B0;
        C0 = C0*QR';
        D0 = real(D0);
        imagLeak = max([norm(imag(E0),'fro'),norm(imag(A0),'fro'), ...
                        norm(imag(B0),'fro'),norm(imag(C0),'fro'),norm(imag(D0),'fro')]);
        E0 = real(E0); A0 = real(A0); B0 = real(B0); C0 = real(C0);
    end

    [Y,X,svH,svV,r] = localProjectionBases(E0,A0,opts);
    Er = Y'*E0*X;
    Ar = Y'*A0*X;
    Br = Y'*B0;
    Cr = C0*X;

    if opts.forceReal
        Er = real(Er); Ar = real(Ar); Br = real(Br); Cr = real(Cr); D0 = real(D0);
    end

    model = struct('E',Er,'A',Ar,'B',Br,'C',Cr,'D',D0,'order',r);
    info = struct('svHorizontal',svH,'svVertical',svV,'imagLeak',imagLeak,'rawSize',[nL nR]);
end

function Q = localRealifyingTransform(idx,s,partners,tolAbs)
    n = numel(idx);
    Q = zeros(n,n);
    q = 1;
    T = [1 1; -1i 1i]/sqrt(2);
    while q <= n
        k = idx(q);
        kp = partners(k);
        if ~(isnan(kp)) && kp ~= k && q < n && idx(q+1) == kp
            Q(q:q+1,q:q+1) = T;
            q = q + 2;
        elseif isnan(kp) || kp == k || abs(s(k)-conj(s(k))) <= tolAbs
            Q(q,q) = 1;
            q = q + 1;
        else
            error('Support indices are not ordered as conjugate blocks.');
        end
    end
end

function [Y,X,svH,svV,r] = localProjectionBases(E0,A0,opts)
    [Uh,Sh] = svd([E0 A0],'econ');
    [~,Sv,Vv] = svd([E0; A0],'econ');
    svH = diag(Sh);
    svV = diag(Sv);
    maxPossible = min([size(Uh,2), size(Vv,2), size(E0,1), size(E0,2)]);

    if isempty(opts.order)
        if isempty(svH), rH = 0; else, rH = sum(svH > opts.svtol*max(svH(1),eps)); end
        if isempty(svV), rV = 0; else, rV = sum(svV > opts.svtol*max(svV(1),eps)); end
        if isinf(opts.maxOrder), rMax = maxPossible; else, rMax = min(opts.maxOrder,maxPossible); end
        r = min([rH,rV,rMax]);
        if r < 1, r = min(1,maxPossible); end
    else
        r = min([opts.order,maxPossible]);
        if r < 1, error('opts.order must be at least 1.'); end
    end

    Y = Uh(:,1:r);
    X = Vv(:,1:r);
end

function metrics = localGridMetrics(model,s,G,Wo,Wi,reps,partners,opts)
    N = numel(s);
    weightedErrAll = zeros(N,1);
    weightedRatioAll = zeros(N,1);
    unweightedAbsAll = zeros(N,1);
    unweightedRelAll = zeros(N,1);

    for k = 1:N
        Rk = localEvalOne(model,s(k)) - G(:,:,k);
        weightedErrAll(k) = localMatNorm(Wo(:,:,k)*Rk*Wi(:,:,k),opts.normType);
        epsk = localSpecTolAt(opts.specTol,k);
        weightedRatioAll(k) = weightedErrAll(k)/max(opts.unweightedFloor,epsk);

        unweightedAbsAll(k) = localMatNorm(Rk,opts.normType);
        Gnorm = localMatNorm(G(:,:,k),opts.normType);
        switch lower(opts.unweightedMetric)
            case 'relative'
                unweightedRelAll(k) = unweightedAbsAll(k)/max(1,Gnorm);
            case 'absolute'
                unweightedRelAll(k) = unweightedAbsAll(k);
            otherwise
                error('opts.unweightedMetric must be relative or absolute.');
        end
    end

    nr = numel(reps);
    weightedErrByRep = zeros(nr,1);
    weightedRatioByRep = zeros(nr,1);
    unweightedRelByRep = zeros(nr,1);
    unweightedAbsByRep = zeros(nr,1);

    for a = 1:nr
        k = reps(a);
        kp = partners(k);
        if isnan(kp) || kp == k
            inds = k;
        else
            inds = [k kp];
        end
        weightedErrByRep(a) = max(weightedErrAll(inds));
        weightedRatioByRep(a) = max(weightedRatioAll(inds));
        unweightedRelByRep(a) = max(unweightedRelAll(inds));
        unweightedAbsByRep(a) = max(unweightedAbsAll(inds));
    end

    metrics = struct();
    metrics.weightedErrAll = weightedErrAll;
    metrics.weightedRatioAll = weightedRatioAll;
    metrics.unweightedAbsAll = unweightedAbsAll;
    metrics.unweightedRelAll = unweightedRelAll;
    metrics.weightedErrByRep = weightedErrByRep;
    metrics.weightedRatioByRep = weightedRatioByRep;
    metrics.unweightedRelByRep = unweightedRelByRep;
    metrics.unweightedAbsByRep = unweightedAbsByRep;
    metrics.maxWeightedErr = max(weightedErrByRep);
    metrics.maxWeightedRatio = max(weightedRatioByRep);
    metrics.maxUnweightedRel = max(unweightedRelByRep);
    metrics.maxUnweightedAbs = max(unweightedAbsByRep);
end

function epsk = localSpecTolAt(specTol,k)
    if isscalar(specTol)
        epsk = specTol;
    else
        epsk = specTol(k);
    end
end

function score = errScoreForSelection(metrics)
    % Placeholder wrapper so the call site is readable.
    score = metrics.weightedRatioByRep;
end

function [aNew,side,desiredSide,selectMode] = localChooseNextBlock( ...
    sideOfRep,~,metrics,leftIdx,rightIdx,nextSide,opts)

    feasible = find(sideOfRep == 0);
    if isempty(feasible)
        aNew = [];
        side = [];
        desiredSide = [];
        selectMode = '';
        return;
    end

    maxWR = metrics.maxWeightedRatio;
    if maxWR > opts.weightedGuard
        selectMode = 'weighted';
        score = metrics.weightedRatioByRep(feasible);
    elseif maxWR > 1
        selectMode = 'mixed';
        wu = metrics.unweightedRelByRep(feasible);
        wu = wu/max(opts.unweightedFloor,max(wu));
        score = metrics.weightedRatioByRep(feasible) + opts.unweightedMixWeight*wu;
    else
        selectMode = 'unweighted';
        score = metrics.unweightedRelByRep(feasible);
    end

    [~,ii] = max(score);
    aNew = feasible(ii);

    switch lower(opts.sideSelection)
        case {'balanced','smaller','minsize'}
            if numel(leftIdx) <= numel(rightIdx)
                side = 1;
            else
                side = 2;
            end
            desiredSide = side;
        case {'alternate','alternating','strictalternate'}
            side = nextSide;
            desiredSide = nextSide;
        otherwise
            error('Unknown opts.sideSelection. Use balanced or alternate.');
    end
end

function directionMode = localDirectionMode(metrics,selectMode,opts)
    if metrics.maxWeightedRatio > opts.weightedGuard
        directionMode = 'weighted';
    elseif metrics.maxWeightedRatio > 1
        directionMode = opts.nearFeasibleDirectionMode;
    else
        directionMode = opts.feasibleDirectionMode;
    end

    if strcmpi(directionMode,'matchfrequency')
        if strcmpi(selectMode,'unweighted')
            directionMode = 'unweighted';
        else
            directionMode = 'weighted';
        end
    end
end

function Hk = localEvalOne(model,sk)
    Hk = model.C*((sk*model.E - model.A)\model.B) + model.D;
end

function val = localMatNorm(M,normType)
    switch lower(normType)
        case {'sigma','spectral','2'}
            val = norm(M,2);
        case {'fro','frobenius'}
            val = norm(M,'fro');
        otherwise
            error('Unknown normType. Use sigma or fro.');
    end
end

function [model,stabInfo] = localStabilizeAndRefit(modelRaw,s,G,Wo,Wi,opts)
    stabInfo = struct();
    stabInfo.applied = false;
    stabInfo.message = '';
    stabInfo.polesBefore = [];
    stabInfo.polesAfter = [];
    stabInfo.regularizedE = false;

    E = modelRaw.E;
    A = modelRaw.A;
    B = modelRaw.B;
    C = modelRaw.C;
    D = modelRaw.D;

    if rcond(E) < opts.EcondTol
        if ~opts.regularizeSingularE
            error('Cannot stabilize: rcond(E)=%.3e < %.3e.',rcond(E),opts.EcondTol);
        end
        [Ue,Se,Ve] = svd(E,'econ');
        se = diag(Se);
        if isempty(se) || se(1) == 0
            error('Cannot regularize singular E: all singular values are zero.');
        end
        keep = se > opts.EsvdTol*se(1);
        if nnz(keep) < 1
            keep(1) = true;
        end
        U1 = Ue(:,keep);
        V1 = Ve(:,keep);
        E = U1'*E*V1;
        A = U1'*A*V1;
        B = U1'*B;
        C = C*V1;
        stabInfo.regularizedE = true;
    end

    As = E\A;
    Bs = E\B;
    Cs = C;

    [X,Lam] = eig(full(As));
    lam = diag(Lam);
    stabInfo.polesBefore = lam;

    unstable = real(lam) >= -opts.stabilityTol;
    lamNew = lam;
    for k = 1:numel(lam)
        if unstable(k)
            margin = max(abs(real(lam(k))),opts.stabilityTol);
            lamNew(k) = -margin + 1i*imag(lam(k));
        end
    end

    Bm = X\Bs;
    Cm = Cs*X;

    if opts.forceReal
        [Ar,Br,Cr] = localModalToReal(lamNew,Bm,Cm,opts);
        Dr = realIfSmallImag(D);
    else
        Ar = diag(lamNew);
        Br = Bm;
        Cr = Cm;
        Dr = D;
    end

    model = struct();
    model.E = eye(size(Ar));
    model.A = Ar;
    model.B = Br;
    model.C = Cr;
    model.D = Dr;
    model.order = size(Ar,1);

    if opts.refitAfterStabilization
        switch lower(opts.stableUpdateMatrix)
            case 'c'
                model = localRefitC(model,s,G,Wo,Wi,opts);
            case 'none'
                % Do nothing.
            otherwise
                error('Currently implemented opts.stableUpdateMatrix values: C or none.');
        end
    end

    if opts.forceReal
        model.E = realIfSmallImag(model.E);
        model.A = realIfSmallImag(model.A);
        model.B = realIfSmallImag(model.B);
        model.C = realIfSmallImag(model.C);
        model.D = realIfSmallImag(model.D);
    end

    stabInfo.applied = true;
    stabInfo.polesAfter = eig(model.A,model.E);
end

function [Ar,Br,Cr] = localModalToReal(lam,Bm,Cm,opts)
    n = numel(lam);
    m = size(Bm,2);
    p = size(Cm,1);
    used = false(n,1);
    Ar = [];
    Br = zeros(0,m);
    Cr = zeros(p,0);
    pairTol = 1e-8*max(1,max(abs(lam)));

    for k0 = 1:n
        if used(k0), continue; end
        lk = lam(k0);
        if abs(imag(lk)) <= pairTol
            Ablk = real(lk);
            Bblk = realIfSmallImag(Bm(k0,:));
            Cblk = realIfSmallImag(Cm(:,k0));
            used(k0) = true;
        else
            k = k0;
            if imag(lam(k)) < 0
                [~,kp] = min(abs(lam - conj(lam(k))));
                if ~used(kp), k = kp; end
            end
            lk = lam(k);
            [dist,j] = min(abs(lam - conj(lk)));
            a = real(lk);
            b = imag(lk);
            brow = Bm(k,:);
            ccol = Cm(:,k);
            Ablk = [a, -b; b, a];
            Bblk = [real(brow); imag(brow)];
            Cblk = [2*real(ccol), -2*imag(ccol)];
            used(k) = true;
            if dist <= pairTol
                used(j) = true;
            end
        end
        Ar = blkdiag(Ar,Ablk);
        Br = [Br; Bblk]; %#ok<AGROW>
        Cr = [Cr, Cblk]; %#ok<AGROW>
    end

    Ar = realIfSmallImag(Ar);
    Br = realIfSmallImag(Br);
    Cr = realIfSmallImag(Cr);
end

function model = localRefitC(model,s,G,Wo,Wi,opts)
    % Refit C in the stable standard realization:
    %   min_C sum_k alpha_k ||G_k - D - C*((s_k I-A)\B)||_F^2
    % with optional scalar weighting alpha_k from hard error weights.
    A = model.A;
    B = model.B;
    D = model.D;
    n = size(A,1);
    N = numel(s);
    [p,m,~] = size(G);

    Z = zeros(n,N*m);
    Y = zeros(p,N*m);

    col = 1;
    for k = 1:N
        Xk = (s(k)*eye(n) - A)\B;
        Yk = G(:,:,k) - D;

        switch lower(opts.stableUpdateFitWeightMode)
            case 'none'
                alpha = 1;
            case 'scalarweighted'
                alpha = sqrt(max(opts.unweightedFloor, ...
                    norm(Wo(:,:,k),'fro')*norm(Wi(:,:,k),'fro')));
            otherwise
                error('stableUpdateFitWeightMode must be none or scalarWeighted.');
        end

        idx = col:(col+m-1);
        Z(:,idx) = alpha*Xk;
        Y(:,idx) = alpha*Yk;
        col = col + m;
    end

    if opts.lsRidge > 0
        ZZ = Z*Z' + opts.lsRidge*eye(size(Z,1));
        Cnew = (Y*Z')/ZZ;
    else
        Cnew = Y/Z;
    end

    model.C = Cnew;
    if opts.fitD
        warning('fitD=true is reserved but not implemented in this C-refit routine; keeping D fixed.');
    end
end

function tf = localIsStable(model,opts)
    lam = eig(model.A,model.E);
    lam = lam(isfinite(lam));
    if isempty(lam)
        tf = false;
    else
        tf = max(real(lam)) < -0.5*opts.stabilityTol;
    end
end

function X = realIfSmallImag(X)
    if norm(imag(X),'fro') <= 1e-10*max(1,norm(real(X),'fro'))
        X = real(X);
    end
end
