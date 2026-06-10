function [wi,error] = iterative_dcf_fa_20190910(iter,adKSpaceCoor,kerneltable,Ind,MatSize,w0,verbose)
% iter:                 # of iterations 
% adKSpaceCoor:         k space coordinates (nx3 or nx2)
% kerneltable:          kernel values (from gridding reconstruction)
% Ind:                  cartesian grid indices (from gridding reconstruction)
% MatSize:              length of cartesian matrix (int) or dimensions (array)
% w0:                   initial DCF estimate (optional)
% verbose:              print progress (optional)


% Adopted from https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.23041
% 3D was modified and needs to be double-checked

% Check for initial DCF
if nargin<6, w0 = []; end
if nargin<7, verbose = []; end
if isempty(w0), w0 = ones(size(adKSpaceCoor,1),1); end
if isempty(verbose), verbose = 1; end

% 2D or 3D?
dim = size(adKSpaceCoor,2);

% Round iter to nearest odd number
% iter = round((iter-1)/2)*2+1;

% Progress ticker
if verbose, fprintf(1,'dcf iter:  1 / %d',iter); end

if dim == 3   
    % Check grid size
    if length(MatSize) == 3
        msize = MatSize;
    else
        msize = [MatSize, MatSize, MatSize];
    end

    % Initial DCF Estimate
    wi = w0;
    kernelweights = zeros(msize);

    for zz = 1:iter           
        % Progress ticker
        if verbose
            fprintf(1,repmat('\b',1,3+numel([num2str(iter) num2str(zz)])));
            fprintf(1,'%d / %d',zz,iter);
        end
        
        % Allocate memory
        grid = zeros(msize);
        
      % % Gridding % %
      
        % Convolve kernel with kspace point
%         witmp = bsxfun(@times,kerneltable,repmat(wi,1,size(kerneltable,2))');
        witmp = (kerneltable .* wi);
        
        for i = 1:size(adKSpaceCoor,1)
            % Get grid indices for each kspace point
            tmpind = Ind(i,:);

            % Combine convolved data
            grid(tmpind) = grid(tmpind) + witmp(i,:);
            
            kernelweights(tmpind) = kernelweights(tmpind) + kerneltable(i,:);
        end
        
% % %         tmpimg = fftshift(ifftn(fftshift(grid)));
% % %         ifftweights = fftshift(ifftn(fftshift(kernelweights)));
% % % %         nonzeroweights = repmat(floor(1+(NumK*OverSampleFactor+1-NumK)/2):floor(1+(NumK*OverSampleFactor+1-NumK)/2)+(NumK-1),3,1)';
% % % %         tmpimg(nonzeroweights) = tmpimg(nonzeroweights) ./ (ifftweights(nonzeroweights).^2);
% % %         tmpimg = tmpimg ./ (ifftweights.^2);
% % %         grid = abs(fftshift(fftn(tmpimg)));

    % % Degridding % %        

        % Allocate memory
        w1i = zeros(size(wi));
        
        % Convolve grid with kernel
        w1itmp = (grid(Ind)) .* kerneltable;

        for i = 1:size(adKSpaceCoor,1)            
            % Combine convolved data
            w1i(i) = sum(w1itmp(i,:));
        end
        
        % Invert density to get weights
        nonzero = w1i ~= 0;
        wi(nonzero) = wi(nonzero) ./ w1i(nonzero);
        wi(~nonzero) = 0;
%         error = mean(1 - w1i.^-1); 
        error = max(w1i.^-1) - min(w1i.^-1);
    end
    
elseif dim == 2
    % Check grid size
    if length(MatSize) == 2
        msize = MatSize;
    else
        msize = [MatSize, MatSize];
    end

    % Initial DCF Estimate
    wi = w0;
    
    for zz = 1:iter
        % Progress ticker
        if verbose
            fprintf(1,repmat('\b',1,3+numel([num2str(iter) num2str(zz)])));
            fprintf(1,'%d / %d',zz,iter);
        end
        
        % Allocate memory
        grid = zeros(msize);
        
      % % Gridding % %
        
        for i = 1:size(adKSpaceCoor,1)
            % Get indices and kernel for each kspace point
            tmpind = Ind(i,:);
            tmpkernel = kerneltable(i,:);
            
            % Convolve kernel with kspace point
            witmp = wi(i) .* tmpkernel;
            
            % Combine convolved data
            grid(tmpind) = grid(tmpind) + witmp;
        end
        
      % % Degridding % %
        
        % Allocate memory
        w1i = zeros(size(wi));
        
        % Convolve grid with kernel
        w1itmp = grid(Ind) .* kerneltable;
%         w1itmp = convn(grid,kerneltable(1,:),'same');

        for i = 1:size(adKSpaceCoor,1)
            % Combine convolved data
            w1i(i) = sum(w1itmp(i,:));
%             w1i(i) = sum(w1itmp(Ind(i,:)));
        end
        
        % Invert density to get weights
        wi = wi .* w1i.^-1;

        %%%%%% JC suggestion to plot a section of wi
    end
end
if verbose, fprintf('\n'); end
end