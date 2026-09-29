function a = pageNorm(b)
% Compute the H_inf norm of b on each page. Similar to Matlab built-in
% pagenorm function, but also squeeze the result into a row.
b2 = pagesvd(b);
a = squeeze(b2(1,1,:)).';
end
