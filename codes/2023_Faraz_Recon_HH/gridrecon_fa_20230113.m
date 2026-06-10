function [Image_out,KSpace_out,wi,Ind,Dist] = gridrecon_fa_20230113(adKSpaceCoor,rawdata,NumK,fov,Options) 
% Look at each sampled kspace point, convolve with kernel, resample onto grid

% parse arguments
arguments
    adKSpaceCoor                                            % nx2 or nx3 kspace coordinates
    rawdata                                                 % raw kspace data (nCol*nLin,max([nSli,nPar]),nCha,nRep)
    NumK                                                    % desired matrix size
    fov                                                     % image FOV
    Options.os double = 3                                   % grid oversampling factor
    Options.zfill double = 1                                % grid zero-filling factor
    Options.k double = 5                                    % kernel size (along 1 dimension)
    Options.verbose double = 0                              % flag for printing progress
    Options.wi double = []                                  % nx1 density correction
    Options.Ind double = []                                 % indices for nearest neighbors
    Options.Dist double = []                                % distances for nearest neighbors
    Options.beta double = []                                % KB kernel beta
    Options.parallel double = 1                             % 1 to search neighbors in parallel (only for large grids)
    Options.lookup double = 1                               % flag for looking up grid
    Options.filename char = 'grid_lookup_20230113.mat'      % filename for saved grids
    Options.savegrid double = 1                             % flag for saving current grid
end

% initialize outputs
Image_out = [];
KSpace_out = [];

% check 2D vs 3D
dimk = size(adKSpaceCoor,2);

% find grid if available
if Options.lookup && (isempty(Options.Ind) && isempty(Options.Dist))
    [Ind,Dist,wi] = grid_lookup_20230113(adKSpaceCoor,NumK,fov,'os',Options.os,'zfill',Options.zfill, ...
        'kernelsize',Options.k,'beta',Options.beta,'wi',Options.wi,'parallel',Options.parallel, ...
        'filename',Options.filename,'save',Options.savegrid,'verbose',Options.verbose);
else
    Ind = Options.Ind;
    Dist = Options.Dist;
    wi = Options.wi;
end


if dimk == 2
    
    % Create cartesian grid
    [x,y] = ndgrid(-ceil(Options.zfill*NumK)/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK)/2/fov - 1/(fov*Options.os));
    W = [x(:) y(:)];
    
    % Find N^2 closest grid points for each sampled kspace point
    if isempty(Ind) || isempty(Dist)
        if Options.verbose, disp('Finding nearest grid points...'); tic; end
        try
            [Ind,Dist] = knnsearch(W,adKSpaceCoor,'K',Options.k^2,'IncludeTies',true);
            Ind = cell2mat(Ind);
            Dist = cell2mat(Dist);
        catch
            [Ind,Dist] = knnsearch(W,adKSpaceCoor,'K',Options.k^2);
        end
        if Options.verbose, fprintf('Nearest points found in %g seconds. \n',toc); end
    end
    
    % Create Kaiser-Bessel Kernel
    if isempty(Options.beta)
        Options.beta = pi*sqrt(Options.k^2/Options.os^2*(Options.os-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
    end
    width = Options.k/(fov*Options.os);
    klength = 10000;
    [kernel,u] = createKBkernel(width,Options.beta,klength);
    kernel2Dtable = interp1(u',kernel,Dist,'linear',0);
    
    % Calculate density correction (if not provided)
    if isempty(wi)
        if Options.verbose, tic; end
        iter = 5;
        dcftable = kernel2Dtable.^(1/2);
        wi = iterative_dcf_fa_20190910(iter,adKSpaceCoor,dcftable,Ind,size(x,1),[],Options.verbose);
        if Options.verbose, fprintf('Density correction finished in %g seconds. \n',toc); end
    end
    
    % Get raw data dimensions
    dim = size(rawdata);
    
    % Check dimensions
    if length(dim) < 4, dim(end+1:4) = 1; end
    
    % Pre-allocate memory
    KSpace_out = complex(zeros([size(W,1),dim(2),dim(3),dim(4)]),0);
    Auxiliary = zeros([size(W,1),dim(2),dim(3),dim(4)]);
    
    % Multiply sampled data with density correction
    M = rawdata.*wi;
    
    % Reshape density correction
    wr = repmat(wi,1,dim(2),dim(3),dim(4));
    
    % Transpose kernel table
    kernel2Dtable = kernel2Dtable';
    
    % Loop through each sampled data point
    if Options.verbose, disp('Adding points onto grid...'); tic; end
    for i = 1:size(adKSpaceCoor,1)
        % Convolve kernel with kspace point
        Mk = M(i,:,:,:).*kernel2Dtable(:,i);
        
        % Add convolved data to grid
        KSpace_out(Ind(i,:),:,:,:) = KSpace_out(Ind(i,:),:,:,:) + Mk;
        Auxiliary(Ind(i,:),:,:,:) = Auxiliary(Ind(i,:),:,:,:) + (wr(i,:,:,:).*kernel2Dtable(:,i));
    end
        
    % Get grid size
    tmpsize = size(x,1);
    
    % Partial de-apodization (ignore regions outside FOV)
    KSpace_out = KSpace_out ./ Auxiliary;
    KSpace_out(isnan(KSpace_out)) = 0;
    
    % FFT kspace to image
    KSpace_out = reshape(KSpace_out,[size(x),dim(2),dim(3),dim(4)]);
    Image_out = KSpace_out;
    for i = 1:2
        Image_out = fftshift(fft(fftshift(Image_out,i),[],i),i);
    end
    
    % Resize image
    if tmpsize >= ceil(NumK*Options.zfill)
        Image_out = Image_out(floor(1+(tmpsize-ceil(NumK*Options.zfill))/2):floor(1+(tmpsize-ceil(NumK*Options.zfill))/2)+(ceil(NumK*Options.zfill)-1),...
                              floor(1+(tmpsize-ceil(NumK*Options.zfill))/2):floor(1+(tmpsize-ceil(NumK*Options.zfill))/2)+(ceil(NumK*Options.zfill)-1),:,:,:);
    else
        Image_out = padarray(Image_out,[floor((ceil(NumK*Options.zfill)-tmpsize)/2),floor((ceil(NumK*Options.zfill)-tmpsize)/2)],0,'both');
    end
    
    % Print
    if Options.verbose, fprintf('Gridding completed in %g seconds. \n',toc); end
    
elseif dimk == 3

    % Check grid size
    if length(NumK) == 1
        NumK = repmat(NumK,[1 3]);
    end

    % Create cartesian grid
    x = (-ceil(Options.zfill*NumK(1))/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK(1))/2/fov - 1/(fov*Options.os))';
    y = (-ceil(Options.zfill*NumK(2))/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK(2))/2/fov - 1/(fov*Options.os))';
    z = (-ceil(Options.zfill*NumK(3))/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK(3))/2/fov - 1/(fov*Options.os))';
    gridszx = size(x,1);
    gridszy = size(y,1);
    gridszz = size(z,1);

    % Find N^2 closest grid points for each sampled kspace point
    if isempty(Ind) || isempty(Dist)
        if Options.verbose, disp('Finding nearest grid points...'); tic; end

        if gridszx > 200 || gridszy > 200 || gridszz > 200 || ~all(NumK==NumK(1))
            Indx = knnsearch(x,adKSpaceCoor(:,1),'K',1);
            Indy = knnsearch(y,adKSpaceCoor(:,2),'K',1);
            Indz = knnsearch(z,adKSpaceCoor(:,3),'K',1);
            k = Options.k;
            Ind = zeros(size(adKSpaceCoor,1),k^2,'single');
            Dist = zeros(size(adKSpaceCoor,1),k^2,'single');
            if Options.parallel
                parfor i = 1:size(adKSpaceCoor,1)
                    [nx,ny,nz] = ndgrid(max(1,Indx(i)-floor(k/2)):min(gridszx,Indx(i)+floor(k/2)), ...
                        max(1,Indy(i)-floor(k/2)):min(gridszy,Indy(i)+floor(k/2)), ...
                        max(1,Indz(i)-floor(k/2)):min(gridszz,Indz(i)+floor(k/2)));
                    [nD,nI] = sort(sqrt(sum((adKSpaceCoor(i,:)-[x(nx(:)),y(ny(:)),z(nz(:))]).^2,2)));
                    Dist(i,:) = nD(1:k^2);
                    Ind(i,:) = sub2ind([gridszx gridszy gridszz], ...
                        nx(nI(1:k^2)), ...
                        ny(nI(1:k^2)), ...
                        nz(nI(1:k^2)));
                end
            else
                for i = 1:size(adKSpaceCoor,1)
                    [nx,ny,nz] = ndgrid(max(1,Indx(i)-floor(k/2)):min(gridszx,Indx(i)+floor(k/2)), ...
                        max(1,Indy(i)-floor(k/2)):min(gridszy,Indy(i)+floor(k/2)), ...
                        max(1,Indz(i)-floor(k/2)):min(gridszz,Indz(i)+floor(k/2)));
                    [nD,nI] = sort(sqrt(sum((adKSpaceCoor(i,:)-[x(nx(:)),y(ny(:)),z(nz(:))]).^2,2)));
                    Dist(i,:) = nD(1:k^2);
                    Ind(i,:) = sub2ind([gridszx gridszy gridszz], ...
                        nx(nI(1:k^2)), ...
                        ny(nI(1:k^2)), ...
                        nz(nI(1:k^2)));
                end
            end
        else
            [x,y,z] = ndgrid(-ceil(Options.zfill*NumK(1))/2/fov:1/(fov*Options.os):ceil(Options.zfill*NumK(1))/2/fov - 1/(fov*Options.os));
            try
                [Ind,Dist] = knnsearch([x(:) y(:) z(:)],adKSpaceCoor,'K',Options.k^2,'IncludeTies',true);
                Ind = cell2mat(Ind);
                Dist = cell2mat(Dist);
            catch
                [Ind,Dist] = knnsearch([x(:) y(:) z(:)],adKSpaceCoor,'K',Options.k^2);
            end
        end

        if Options.verbose, fprintf('Nearest points found in %g seconds. \n',toc); end
    end


    
    % Create Kaiser-Bessel Kernel
    if isempty(Options.beta)
        Options.beta = pi*sqrt(Options.k^2/Options.os^2*(Options.os-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
    end
    width = Options.k/(fov*Options.os);
    klength = 10000;
    [kernel,u] = createKBkernel(width,Options.beta,klength);
    kernel3Dtable = interp1(u',kernel,Dist,'linear',0);
    
    % Calculate density correction (if not provided)
    if isempty(wi)
        if Options.verbose, tic; end
        iter = 5;
        dcftable = kernel3Dtable.^(1/2);
        wi = iterative_dcf_fa_20190910(iter,adKSpaceCoor,dcftable,Ind,[gridszx gridszy gridszz],[],Options.verbose);
        if Options.verbose, fprintf('Density correction finished in %g seconds. \n',toc); end
    end
    
    % Get raw data dimensions
    dim = size(rawdata);
    
    % Check dimensions
    if length(dim) < 3, dim(end+1:3) = 1; end
    
    % Pre-allocate memory
    k_real = zeros([gridszx*gridszy*gridszz,dim(2),dim(3)],'single');
    k_imag = zeros([gridszx*gridszy*gridszz,dim(2),dim(3)],'single');
    Auxiliary = zeros([gridszx*gridszy*gridszz,dim(2),dim(3)],'single');
    
    % Multiply sampled data with density correction
    M = single(rawdata.*wi);
    
    % Reshape density correction
    wr = single(repmat(wi,1,dim(2),dim(3)));
    
    % Transpose kernel table
    kernel3Dtable = single(kernel3Dtable');
    
    % Loop through each sampled data point
    if Options.verbose, disp('Adding points onto grid...'); tic; end
    for i = 1:size(adKSpaceCoor,1)
        % Convolve kernel with kspace point
        Mk = M(i,:,:).*kernel3Dtable(:,i);
        
        % Add convolved data to grid
        k_real(Ind(i,:),:,:) = k_real(Ind(i,:),:,:) + real(Mk);
        k_imag(Ind(i,:),:,:) = k_imag(Ind(i,:),:,:) + imag(Mk);
        Auxiliary(Ind(i,:),:,:) = Auxiliary(Ind(i,:),:,:) + (wr(i,:,:).*kernel3Dtable(:,i));
    end
    
    % Combine real and imag components
    KSpace_out = complex(k_real,k_imag);
    
    % Partial de-apodization (ignore regions outside FOV)
    KSpace_out = KSpace_out ./ Auxiliary;
    KSpace_out(isnan(KSpace_out)) = 0;
    
    % FFT kspace to image
    KSpace_out = reshape(KSpace_out,[gridszx,gridszy,gridszz,dim(2),dim(3)]);
    Image_out = KSpace_out;
    for i = 1:3
        Image_out = fftshift(fft(fftshift(Image_out,i),[],i),i);
    end
    
    % Resize image
    if Options.os > 1
        Image_out = Image_out(floor(1+(gridszx-ceil(NumK(1)*Options.zfill))/2):floor(1+(gridszx-ceil(NumK(1)*Options.zfill))/2)+(ceil(NumK(1)*Options.zfill)-1),...
                              floor(1+(gridszy-ceil(NumK(2)*Options.zfill))/2):floor(1+(gridszy-ceil(NumK(2)*Options.zfill))/2)+(ceil(NumK(2)*Options.zfill)-1),...
                              floor(1+(gridszz-ceil(NumK(3)*Options.zfill))/2):floor(1+(gridszz-ceil(NumK(3)*Options.zfill))/2)+(ceil(NumK(3)*Options.zfill)-1),:,:);
    elseif Options.os < 1
        Image_out = padarray(Image_out,...
            [floor((ceil(NumK(1)*Options.zfill)-gridszx)/2),...
            floor((ceil(NumK(2)*Options.zfill)-gridszy)/2),...
            floor((ceil(NumK(3)*Options.zfill)-gridszz)/2)],0,'both');
    end
    
    % Print
    if Options.verbose, fprintf('Gridding completed in %g seconds. \n',toc); end
    
end
end
