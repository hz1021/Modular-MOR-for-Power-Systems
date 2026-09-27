function [s,G,Wo,Wi, WoErrBound, WiErrBound] = mirrorMimoWeightedFrequencyData(w,Gpos,WoPos,WiPos, WoErrBound, WiErrBound)
%MIRRORMIMOWEIGHTEDFREQUENCYDATA Add negative-frequency conjugate samples.
%
%   [S,G,WO,WI] = MIRRORMIMOWEIGHTEDFREQUENCYDATA(W,GPOS,WOPOS,WIPOS)
%   converts positive-frequency samples into conjugate-symmetric samples
%   suitable for a real continuous-time realization.
%
%   Inputs:
%       W          n-by-1 nonnegative frequencies in rad/s.
%       GPOS       p-by-m-by-n,      GPOS(:,:,k)  = G(1i*W(k)).
%       WOPOS      qo-by-p-by-n or [], output-weight samples at 1i*W(k).
%       WIPOS      m-by-ri-by-n or [], input-weight samples at 1i*W(k).
%
%   Outputs include positive frequencies, optional zero frequency, and the
%   conjugate negative-frequency data:
%
%       G(-1i*w)  = conj(G(1i*w)),
%       Wo(-1i*w) = conj(Wo(1i*w)),
%       Wi(-1i*w) = conj(Wi(1i*w)).

    w = w(:);
    if any(w < 0)
        error('w must contain nonnegative frequencies only.');
    end
    [p,m,n] = size(Gpos);
    if numel(w) ~= n
        error('numel(w) must equal size(Gpos,3).');
    end

    if isempty(WoPos)
        WoPos = repmat(eye(p),1,1,n);
    elseif ndims(WoPos) == 2
        WoPos = repmat(WoPos,1,1,n);
    end
    if isempty(WiPos)
        WiPos = repmat(eye(m),1,1,n);
    elseif ndims(WiPos) == 2
        WiPos = repmat(WiPos,1,1,n);
    end

    if size(WoPos,2) ~= p || size(WoPos,3) ~= n
        error('WoPos must be qo-by-p-by-n.');
    end
    if size(WiPos,1) ~= m || size(WiPos,3) ~= n
        error('WiPos must be m-by-ri-by-n.');
    end

    tol = 100*eps(max(1,max(w)));
    iz = find(abs(w) <= tol);
    ip = find(w > tol);

    % Sort positive frequencies for cleaner plots and support reporting.
    [~,ord] = sort(w(ip));
    ip = ip(ord);

    sPos = 1i*w(ip);
    Gp = Gpos(:,:,ip);
    Wop = WoPos(:,:,ip);
    Wip = WiPos(:,:,ip);

    Wopeb = WoErrBound(:,:,ip);
    Wipeb = WiErrBound(:,:,ip);

    if isempty(iz)
        s = [sPos; conj(sPos)];
        G = cat(3,Gp,conj(Gp));
        Wo = cat(3,Wop,conj(Wop));
        Wi = cat(3,Wip,conj(Wip));

        WoErrBound = cat(3,Wopeb,conj(Wopeb));
        WiErrBound = cat(3,Wipeb,conj(Wipeb));
    else
        z = iz(1);
        s = [sPos; 0; conj(sPos)];
        G0 = real(Gpos(:,:,z));
        Wo0 = real(WoPos(:,:,z));
        Wi0 = real(WiPos(:,:,z));

        Wo0eb = real(WoErrBound(:,:,z));
        Wi0eb = real(WiErrBound(:,:,z));

        G = cat(3,Gp,G0,conj(Gp));
        Wo = cat(3,Wop,Wo0,conj(Wop));
        Wi = cat(3,Wip,Wi0,conj(Wip));

        WoErrBound = cat(3,Wopeb,Wo0eb,conj(Wopeb));
        WiErrBound = cat(3,Wipeb,Wi0eb,conj(Wipeb));
    end
end
