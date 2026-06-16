function run_cs_sweep(in_mat, out_mat, codes_root, tvweights)
% run_cs_sweep(in_mat, out_mat, codes_root, tvweights)
%
% NEW file — does NOT modify run_cs.m. Sweeps Lustig fnlCg over a vector of TV
% weights so the Lustig baseline matches BART's lambda sweep (handoff step 2).
% Pure TV: xfmWeight = 0 (the L1 term is OFF), so this is a clean single-
% regularizer TV run comparable to BART `-R T` and our PDHG-TV. Identity XFM is
% irrelevant here (only used by the disabled L1 term).
%
% Lustig's wavelet operator is FWT2_PO (2D only) — a true 3D-wavelet Lustig run
% is not available without a new operator, so the wavelet comparison stays
% ours-vs-BART. This sweep covers TV only, on purpose.
%
% Math otherwise byte-for-byte run_cs.m: NUFFT3D(k,1,1,0,[100 100 100],2),
% w=kcomps/max, data=kdata/max(abs), 15 fnlCg iters, per-slice rot90, per-iter
% max-norm. Saves only the FINAL iterate per weight to keep the file small.
%
%   tvweights : row vector, e.g. [1e-3 3e-3 1e-2 3e-2 1e-1]
%   out_mat   : saves gas_final(numW,100,100,100) and tvW(1,numW)

run(fullfile(codes_root, '2025_CS', 'irt', 'setup.m'));
addpath(genpath(fullfile(codes_root, '2025_CS', 'sparseMRI_v0.2')));

msnow = 100;
data_frames = load(in_mat);
ktrajs = data_frames.ktrajs;
kdatas = data_frames.kdatas;
kcomps = data_frames.kcomps;

frame = 1;
dkx = squeeze(ktrajs(frame,1,:))/(2*pi);   % [-pi,pi] -> [-0.5,0.5]
dky = squeeze(ktrajs(frame,2,:))/(2*pi);
dkz = squeeze(ktrajs(frame,3,:))/(2*pi);
k = [dkx, dky, dkz];

imSize = [msnow, msnow, msnow];
FT = NUFFT3D(k, 1, 1, 0, imSize, 2);

w = kcomps(frame,:)';  w = w/max(w(:));
data = kdatas(frame,:)';  data = data/max(abs(data));
im_dc = FT'*(data.*w);

numW = numel(tvweights);
gas_final = zeros(numW, msnow, msnow, msnow);
tvW = tvweights(:)';

for wi = 1:numW
    param = init;
    param.FT = FT;
    param.XFM = 1;              % identity (L1 term off, so unused)
    param.TV = TVOP3D;
    param.data = data;
    param.TVWeight = tvweights(wi);
    param.xfmWeight = 0;        % pure TV
    param.Itnlim = 15;

    res = im_dc;
    for n = 1:15
        res = fnlCg(res, param);
    end
    img = abs(res);
    for jj = 1:size(img,1)
        img(jj,:,:) = rot90(squeeze(img(jj,:,:)));
    end
    img = img/max(img(:));
    gas_final(wi,:,:,:) = img;
    fprintf('  TVWeight %.4g done\n', tvweights(wi));
end

save(out_mat, 'gas_final', 'tvW', '-v7');
fprintf('saved %s : gas_final(%d,%d,%d,%d)\n', out_mat, ...
        size(gas_final,1), size(gas_final,2), size(gas_final,3), size(gas_final,4));
end
