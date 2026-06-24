function [Image_out,KSpace_out,wi,Ind,Dist] = gridrecon_fa_20220210(adKSpaceCoor,rawdata,NumK,fov,zfill,GridOverSampleFactor,N,wi,verbose,Ind,Dist)
% Look at each sampled kspace point, convolve with kernel, resample onto grid
% adKSpaceCoor:             nx2 kspace coordinates
% rawdata:                  raw kspace data (nCol*nLin,max([nSli,nPar]),nCha,nRep)
% NumK:                     desired matrix size
% FOV:                      image FOV
% zfill                     zero filling factor for gridding
% GridOverSampleFactor:     oversample factor for gridding
% N:                        kernel size (along 1 dimension) for gridding
% wi:                       nx1 density correction (optional)
% verbose:                  print progress (optional)
% Ind:                      indices for nearest neighbors (optional)
% Dist:                     distances for nearest neighbors (optional)

if nargin<5, zfill = []; end
if nargin<6, GridOverSampleFactor = []; end
if nargin<7, N = []; end
if nargin<8, wi = []; end
if nargin<9, verbose = []; end
if nargin<10, Ind = []; end
if nargin<11, Dist = []; end
if isempty(zfill) || zfill < 1, zfill = 1; end
if isempty(GridOverSampleFactor), GridOverSampleFactor = 3; end
if isempty(N), N = 5; end
if isempty(verbose), verbose = 1; end

% initialize outputs
Image_out = [];
KSpace_out = [];

% check 2D vs 3D
dimk = size(adKSpaceCoor,2);

if dimk == 2
    
    % Create cartesian grid
    [x,y] = ndgrid(-ceil(zfill*NumK)/2/fov:1/(fov*GridOverSampleFactor):ceil(zfill*NumK)/2/fov);
    W = [x(:) y(:)];
    
    % Find N^2 closest grid points for each sampled kspace point
    if isempty(Ind) || isempty(Dist)
        if verbose, disp('Finding nearest grid points...'); end
        try
            [Ind,Dist] = knnsearch(W,adKSpaceCoor,'K',N^2,'IncludeTies',true);
            Ind = cell2mat(Ind);
            Dist = cell2mat(Dist);
        catch
            [Ind,Dist] = knnsearch(W,adKSpaceCoor,'K',N^2);
        end
    end
    
    % Create Kaiser-Bessel Kernel
    beta = pi*sqrt(N^2/GridOverSampleFactor^2*(GridOverSampleFactor-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
    width = 1*N/(fov*GridOverSampleFactor);
    klength = 10000;
    [kernel,u] = createKBkernel(width,beta,klength);
    kernel2Dtable = interp1(u',kernel,Dist,'linear',0);
    
    % Calculate density correction (if not provided)
    if isempty(wi)
        iter = 5;
        dcftable = kernel2Dtable.^(1/2);
        wi = iterative_dcf_fa_20190910(iter,adKSpaceCoor,dcftable,Ind,size(x,1),[],verbose);
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
    if verbose, tic; end
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
        Image_out = fftshift(ifft(Image_out,[],i),i);
    end
    
    % Resize image
    if tmpsize >= ceil(NumK*zfill)
        Image_out = Image_out(floor(1+(tmpsize-ceil(NumK*zfill))/2):floor(1+(tmpsize-ceil(NumK*zfill))/2)+(ceil(NumK*zfill)-1),...
            floor(1+(tmpsize-ceil(NumK*zfill))/2):floor(1+(tmpsize-ceil(NumK*zfill))/2)+(ceil(NumK*zfill)-1),:,:,:);
    else
        Image_out = padarray(Image_out,[floor((ceil(NumK*zfill)-tmpsize)/2),floor((ceil(NumK*zfill)-tmpsize)/2)],0,'both');
    end
    
    % Print
    if verbose, fprintf('Gridding completed in %g seconds. \n',toc); end
    
elseif dimk == 3
    % Create cartesian grid
    [x,y,z] = ndgrid(-ceil(zfill*NumK)/2/fov:1/(fov*GridOverSampleFactor):ceil(zfill*NumK)/2/fov);
    W = [x(:) y(:) z(:)];
    
    % Find N^3 closest grid points for each sampled kspace point
    if isempty(Ind) || isempty(Dist)
        if verbose, disp('Finding nearest grid points...'); end
        try
            [Ind,Dist] = knnsearch(W,adKSpaceCoor,'K',N^2,'IncludeTies',true);
            Ind = cell2mat(Ind);
            Dist = cell2mat(Dist);
        catch
            [Ind,Dist] = knnsearch(W,adKSpaceCoor,'K',N^2);
        end
    end
    
    % Create Kaiser-Bessel Kernel
    beta = pi*sqrt(N^2/GridOverSampleFactor^2*(GridOverSampleFactor-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
    width = 1*N/(fov*GridOverSampleFactor);
    klength = 10000;
    [kernel,u] = createKBkernel(width,beta,klength);
    kernel3Dtable = interp1(u',kernel,Dist,'linear',0);
    
    % Calculate density correction (if not provided)
    if isempty(wi)
        iter = 5;
        dcftable = kernel3Dtable.^(1/2);
        wi = iterative_dcf_fa_20190910(iter,adKSpaceCoor,dcftable,Ind,size(x,1),[],verbose);
    end
    
    % Get raw data dimensions
    dim = size(rawdata);
    
    % Check dimensions
    if length(dim) < 3, dim(end+1:3) = 1; end
    
    % Pre-allocate memory
    k_real = zeros([size(W,1),dim(2),dim(3)],'single');
    k_imag = zeros([size(W,1),dim(2),dim(3)],'single');
    Auxiliary = zeros([size(W,1),dim(2),dim(3)],'single');
    
    % Multiply sampled data with density correction
    M = single(rawdata.*wi);
    
    % Reshape density correction
    wr = single(repmat(wi,1,dim(2),dim(3)));
    
    % Transpose kernel table
    kernel3Dtable = single(kernel3Dtable');
    
    % Loop through each sampled data point
    if verbose, tic; end
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
    
    % Get grid size
    tmpsize = size(x,1);
    
    % Partial de-apodization (ignore regions outside FOV)
    KSpace_out = KSpace_out ./ Auxiliary;
    KSpace_out(isnan(KSpace_out)) = 0;
    
    % FFT kspace to image
    KSpace_out = reshape(KSpace_out,[size(x),dim(2),dim(3)]);
    Image_out = KSpace_out;
    for i = 1:3
        Image_out = fftshift(ifft(Image_out,[],i),i);
    end
    
    % Resize image
    if tmpsize >= ceil(NumK*zfill)
        Image_out = Image_out(floor(1+(tmpsize-ceil(NumK*zfill))/2):floor(1+(tmpsize-ceil(NumK*zfill))/2)+(ceil(NumK*zfill)-1),...
                              floor(1+(tmpsize-ceil(NumK*zfill))/2):floor(1+(tmpsize-ceil(NumK*zfill))/2)+(ceil(NumK*zfill)-1),...
                              floor(1+(tmpsize-ceil(NumK*zfill))/2):floor(1+(tmpsize-ceil(NumK*zfill))/2)+(ceil(NumK*zfill)-1),:,:);
    else
        Image_out = padarray(Image_out,...
            [floor((ceil(NumK*zfill)-tmpsize)/2),...
            floor((ceil(NumK*zfill)-tmpsize)/2),...
            floor((ceil(NumK*zfill)-tmpsize)/2)],0,'both');
    end
    
    % Print
    if verbose, fprintf('Gridding completed in %g seconds. \n',toc); end
    
end
end
