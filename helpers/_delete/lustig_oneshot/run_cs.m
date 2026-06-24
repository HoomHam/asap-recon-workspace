function run_cs(in_mat, out_mat, codes_root)
% run_cs(in_mat, out_mat, codes_root)
%
% Headless, exact reproduction of spiral3d_cs_3D_hoom.m (TV+L1, identity XFM).
% Reads ACR_test-format in_mat (ktrajs/kdatas/kcomps), runs 15 fnlCg iters,
% saves gas(15,100,100,100) to out_mat. No figures (batch-safe). Math is
% byte-for-byte the original: NUFFT3D(k,1,1,0,[100 100 100],2), w=kcomps/max,
% data=kdata/max(abs), TVWeight=xfmWeight=0.01, Itnlim=15, per-slice rot90 +
% per-iter max-normalization, exactly as the script that produced tv01g01.mat.
%
% codes_root = workspace/codes (holds 2025_CS/irt + 2025_CS/sparseMRI_v0.2)

% --- paths (IRT NUFFT + Lustig sparseMRI) ---
run(fullfile(codes_root, '2025_CS', 'irt', 'setup.m'));
addpath(genpath(fullfile(codes_root, '2025_CS', 'sparseMRI_v0.2')));

MS_recon_desired = 100;
msnow = MS_recon_desired;

data_frames = load(in_mat);
ktrajs = data_frames.ktrajs;
kdatas = data_frames.kdatas;
kcomps = data_frames.kcomps;

frame = 1;
dkx = ktrajs(frame, 1, :);
dky = ktrajs(frame, 2, :);
dkz = ktrajs(frame, 3, :);

% [-pi,pi] -> [-0.5,0.5]
dkx = dkx(:)/(2*pi);
dky = dky(:)/(2*pi);
dkz = dkz(:)/(2*pi);
k = [dkx, dky, dkz];

ph = 1;
slices = msnow;
imSize = [msnow, slices, msnow];
FT = NUFFT3D(k, 1, ph, 0, imSize, 2);

w = kcomps(frame,:)';
w = w/max(w(:));

data = kdatas(frame,:)';
data = data / max(abs(data));   % scale kdata, NOT the initial image

im_dc = FT'*(data.*w);

%%%% CS
XFM = 1;            % Identity transform
TVWeight = 0.01;    % TV penalty
xfmWeight = 0.01;   % L1 penalty
Itnlim = 15;

param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP3D;
param.data = data;
param.TVWeight = TVWeight;
param.xfmWeight = xfmWeight;
param.Itnlim = Itnlim;

res = XFM*im_dc;

gas = zeros(15, msnow, slices, msnow);
tic
for n = 1:15
    res = fnlCg(res, param);
    im_res = XFM'*res;
    img_gp = abs(im_res);
    for jj = 1:size(img_gp,1)
        tmp = squeeze(img_gp(jj,:,:));
        img_gp(jj,:,:) = rot90(tmp);
    end
    img_gp = img_gp/max(img_gp(:));
    gas(n,:,:,:) = img_gp;
end
toc

save(out_mat, 'gas');   % v7 (gas ~120MB < 2GB) — scipy.io.loadmat readable
fprintf('saved %s : gas(%d,%d,%d,%d)\n', out_mat, size(gas,1), size(gas,2), size(gas,3), size(gas,4));
end
