function [Ind,Dist,wi] = grid_lookup(KSpaceCoor,imgsize,fov,Options)
% load or save spiral grid parameters based on input trajectory/options

% parse arguments
arguments
    KSpaceCoor
    imgsize
    fov
    Options.os double = 3
    Options.zfill double = 1
    Options.kernelsize double = 5
    Options.beta double = []
    Options.wi double = []
    Options.parallel double = 1
    Options.filename char = 'grid_lookup_20220418.mat'
    Options.date char = datestr(now,'yyyy-mm-dd')
    Options.save double = 1
    Options.verbose double = 1
end

% initialize outputs
Ind = [];
Dist = [];
wi = [];

% check arguments
if isempty(Options.beta)
    Options.beta = pi*sqrt(Options.kernelsize^2/Options.os^2*(Options.os-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
end

% find file
fileloc = which(Options.filename);

% load grids
s = load(fileloc);
s = s.s;

% initial search
match1 = arrayfun(@(x) ...
    all(size(x.KSpaceCoor) == size(KSpaceCoor))         && ...
    abs(x.imgsize - imgsize)                    <= 0    && ...
    abs(x.fov - fov)                            <= 0    && ...
    abs(x.os - Options.os)                      <= 0    && ...
    abs(x.kernelsize - Options.kernelsize)      <= 0    && ...
    abs(x.beta - Options.beta)                  <= 0    && ...
    abs(x.zfill - Options.zfill)                <= 0    ...
    , s);

% filter out incorrect sized trajectories
s2 = s(match1);

% search for saved grid
match = arrayfun(@(x) ...
    abs(sum(x.KSpaceCoor(:) - single(KSpaceCoor(:)))) <= 0, s2);

% get index for last (most recent) match
i = find(match == 1,1,'last');

% get grid params if available or calculate
if isempty(i)
    
    % get last index
    i = 1+size(s,2);
    
    % check 2D vs 3D
    dimk = size(KSpaceCoor,2);
    
    % gridding
    if dimk == 2
        % create cartesian grid
        [x,y] = ndgrid(-ceil(Options.zfill*imgsize)/2/fov:1/(fov*Options.os):ceil(Options.zfill*imgsize)/2/fov);
        W = [x(:) y(:)];
        
        % find N^2 closest grid points for each sampled kspace point
        if Options.verbose, disp('Finding nearest grid points...'); tic; end
        try
            [Ind,Dist] = knnsearch(W,KSpaceCoor,'K',Options.kernelsize^2,'IncludeTies',true);
            Ind = cell2mat(Ind);
            Dist = cell2mat(Dist);
        catch
            [Ind,Dist] = knnsearch(W,KSpaceCoor,'K',Options.kernelsize^2);
        end
        if Options.verbose, fprintf('Nearest points found in %g seconds. \n',toc); end
        
        % create Kaiser-Bessel Kernel
%         beta = pi*sqrt(Options.kernelsize^2/Options.os^2*(Options.os-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
        width = Options.kernelsize/(fov*Options.os);
        klength = 10000;
        [kernel,u] = createKBkernel(width,Options.beta,klength);
        kernel2Dtable = interp1(u',kernel,Dist,'linear',0);
        
        % calculate density correction (if not provided)
        if isempty(Options.wi)
            if Options.verbose, tic; end
            iter = 5;
            dcftable = kernel2Dtable.^(1/2);
            wi = iterative_dcf_fa_20190910(iter,KSpaceCoor,dcftable,Ind,size(x,1),[],Options.verbose);
            if Options.verbose, fprintf('Density correction finished in %g seconds. \n',toc); end
        else
            wi = Options.wi;
        end
        
    elseif dimk == 3

        % create cartesian grid
        x = (-ceil(Options.zfill*imgsize)/2/fov:1/(fov*Options.os):ceil(Options.zfill*imgsize)/2/fov)';
        gridsz = size(x,1);

        % find N^2 closest grid points for each sampled kspace point
        if Options.verbose, disp('Finding nearest grid points...'); tic; end

        if gridsz > 200
            Indx = knnsearch(x,KSpaceCoor(:,1),'K',1);
            Indy = knnsearch(x,KSpaceCoor(:,2),'K',1);
            Indz = knnsearch(x,KSpaceCoor(:,3),'K',1);
            k = Options.kernelsize;
            Ind = zeros(size(KSpaceCoor,1),k^2,'single');
            Dist = zeros(size(KSpaceCoor,1),k^2,'single');
            if Options.parallel
                parfor ii = 1:size(KSpaceCoor,1)
                    [nx,ny,nz] = ndgrid(max(1,Indx(ii)-floor(k/2)):min(gridsz,Indx(ii)+floor(k/2)), ...
                                        max(1,Indy(ii)-floor(k/2)):min(gridsz,Indy(ii)+floor(k/2)), ...
                                        max(1,Indz(ii)-floor(k/2)):min(gridsz,Indz(ii)+floor(k/2)));
                    [nD,nI] = sort(sqrt(sum((KSpaceCoor(ii,:)-x([nx(:),ny(:),nz(:)])).^2,2)));
                    Dist(ii,:) = nD(1:k^2);
                    Ind(ii,:) = sub2ind([gridsz gridsz gridsz], ...
                        nx(nI(1:k^2)), ...
                        ny(nI(1:k^2)), ...
                        nz(nI(1:k^2)));
                end
            else
                for ii = 1:size(KSpaceCoor,1)
                    [nx,ny,nz] = ndgrid(max(1,Indx(ii)-floor(k/2)):min(gridsz,Indx(ii)+floor(k/2)), ...
                                        max(1,Indy(ii)-floor(k/2)):min(gridsz,Indy(ii)+floor(k/2)), ...
                                        max(1,Indz(ii)-floor(k/2)):min(gridsz,Indz(ii)+floor(k/2)));
                    [nD,nI] = sort(sqrt(sum((KSpaceCoor(ii,:)-x([nx(:),ny(:),nz(:)])).^2,2)));
                    Dist(ii,:) = nD(1:k^2);
                    Ind(ii,:) = sub2ind([gridsz gridsz gridsz], ...
                        nx(nI(1:k^2)), ...
                        ny(nI(1:k^2)), ...
                        nz(nI(1:k^2)));
                end
            end
        else
            [x,y,z] = ndgrid(-ceil(Options.zfill*imgsize)/2/fov:1/(fov*Options.os):ceil(Options.zfill*imgsize)/2/fov);
            try
                [Ind,Dist] = knnsearch([x(:) y(:) z(:)],KSpaceCoor,'K',Options.kernelsize^2,'IncludeTies',true,'NSMethod','kdtree');
                Ind = cell2mat(Ind);
                Dist = cell2mat(Dist);
            catch
                [Ind,Dist] = knnsearch([x(:) y(:) z(:)],KSpaceCoor,'K',Options.kernelsize^2);
            end
        end

        if Options.verbose, fprintf('Nearest points found in %g seconds. \n',toc); end

        % create Kaiser-Bessel Kernel
%         beta = pi*sqrt(Options.kernelsize^2/Options.os^2*(Options.os-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
        width = Options.kernelsize/(fov*Options.os);
        klength = 10000;
        [kernel,u] = createKBkernel(width,Options.beta,klength);
        kernel3Dtable = interp1(u',kernel,Dist,'linear',0);
        
        % calculate density correction (if not provided)
        if isempty(Options.wi)
            if Options.verbose, tic; end
            iter = 5;
            dcftable = kernel3Dtable.^(1/2);
            wi = iterative_dcf_fa_20190910(iter,KSpaceCoor,dcftable,Ind,size(x,1),[],Options.verbose);
            if Options.verbose, fprintf('Density correction finished in %g seconds. \n',toc); end
        else
            wi = Options.wi;
        end
    end
    
    % add new grid
    if Options.save
        s(i).KSpaceCoor = single(KSpaceCoor);
        s(i).fov = fov;
        s(i).imgsize = imgsize;
        s(i).os = Options.os;
        s(i).kernelsize = Options.kernelsize;
        s(i).beta = Options.beta;
        s(i).zfill = Options.zfill;
        s(i).wi = single(wi);
        s(i).Ind = single(Ind);
        s(i).Dist = single(Dist);
        s(i).date = Options.date;
        
        % update .mat file
        save(fileloc,'s');
        if Options.verbose
            disp(['Grid saved to ',fileloc]);
        end
    end
    
else
    % load grid
    Ind = s2(i).Ind;
    Dist = s2(i).Dist;
    wi = s2(i).wi;
    if Options.verbose
        disp(['Grid from ',s2(i).date,' loaded']);
    end
end

end