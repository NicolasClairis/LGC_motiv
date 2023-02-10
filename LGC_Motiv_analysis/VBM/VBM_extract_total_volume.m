%% select c1, c2 and c3 to create total volume
P=spm_select(Inf,'nifti');

N=nifti(P);
vol = zeros(numel(N),1);
for iS = 1:numel(N)
    dat = N(iS).dat(:,:,:);
    dat = dat(isfinite(dat));
    s = sum(dat(:));
    vv = abs(det(N(iS).mat(1:3,1:3)))/(100^3);
    % extract volume of c1 (grey matter), c2 (white matter) and c3 (CSF)
    vol(iS) = s*vv;
end

%% at the end you can use vol as a co-variate for your VBM analysis
