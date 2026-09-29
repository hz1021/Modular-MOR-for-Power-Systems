function H = evalMimoLoewnerModel(model,s)
%   Evaluate descriptor Loewner model on sample points
%
%   H = EVALMIMOLOEWNERMODEL(MODEL,S) returns H(:,:,k) = Gr(S(k)), where
%   Gr(s) = C*((s*E - A)\B) + D

    s = s(:);
    p = size(model.C,1);
    m = size(model.B,2);
    H = zeros(p,m,numel(s));
    for k = 1:numel(s)
        H(:,:,k) = model.C*((s(k)*model.E - model.A)\model.B) + model.D;
    end
end
