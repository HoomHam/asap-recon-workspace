function [img, kspace, rawimg, rawkspace] = fa_spiral_dyn_recon(rawdata,KSpaceCoor,wi,acq_params,reco_params,progbar)
% rawdata               twix raw data (nsamples*nleaves,max([nslices,npartitions]),nchannels,nechoes,nrepetitions)
% KSpaceCoor            kspace coordines (nsamples*nleaves,2)
% wi                    density correction (nsamples*nleaves,1)
% acq_params            "allimages" structure with necessary parameters
% reco_params           structure with reconctruction/gridding parameters
% progbar               structure with progressbar info

% unpack required acq_param values
protocol = acq_params.protocol;
nsamples = acq_params.nCol;
nleaves = acq_params.nLin;
nslices = acq_params.nSli;
npartitions = acq_params.nPar;
nchannels = acq_params.nCha;
nrepetitions = acq_params.nRep;
acqOrder = acq_params.acqOrder;
numRes = acq_params.numRes;
if isfield(acq_params,'GA')
    GA_flag = acq_params.GA;
    if GA_flag == 1
        GAperiod = acq_params.GAperiod;
    end
else
    if isfield(acq_params,'GAperiod')
        GAperiod = acq_params.GAperiod;
    end
end

% unpack reco_params
imgsize = reco_params.imgsize;
fovPE = reco_params.fovPE;
zfill = reco_params.zfill;
GridOSFactor = reco_params.GridOSFactor;
KernelSize = reco_params.KernelSize;
savekspace = reco_params.savekspace;
saveraw = reco_params.saveraw;
if isfield(reco_params,'nreps')
    nreps = reco_params.nreps;
end

% unpack progbar
if nargin < 6, progbar = []; end
if ~isempty(progbar)
    h = progbar.h;
    enum = progbar.enum;
    N = progbar.N;
    progressbar = 1;
else
    progressbar = 0;
end

% verbose flag
verbose = 0;

% initialize
kspace = [];
rawimg = [];
rawkspace = [];

% recon based on SOS vs Fancy spirals
if contains(protocol,'fancy_v2') || contains(protocol,'fancy_v3')

    % reshape
    rawdata = permute(reshape(rawdata,nsamples,nleaves,max([nslices,npartitions]),nchannels,nrepetitions),[1 4 2 5 3]);
    rawdata = reshape(rawdata,nsamples,nchannels,[]);

    % get indices for each interleave
    intind = [];
    for j = 1:ceil(size(rawdata,3)/numRes)
        intind = [intind; repmat(j,numRes,1)];
    end
    intind(size(rawdata,3)+1:end) = [];

    % initialize combined rawdata
    rawdata2 = zeros(nsamples,nleaves,nreps,nchannels);
    weights = zeros(nsamples,nleaves,nreps,nchannels);

    % combine rawdata
    for i = 1:size(rawdata,3)
        rep = 0;
        while intind(i) - (rep+1)*nleaves > 0
            rep = rep + 1;
        end
        int = intind(i) - (rep)*nleaves;
        nrep = mod(rep,nreps)+1;

        rawdata2(:,int,nrep,:) = rawdata2(:,int,nrep,:) + reshape(rawdata(:,:,i),nsamples,1,1,nchannels);
        weights(:,int,nrep,:) = weights(:,int,nrep,:) + 1;
    end
    rawdata2 = reshape(rawdata2./weights,nsamples*nleaves*nreps,nchannels);

    % update progressbar
    if progressbar, waitbar((0.25+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',25)); end

    % get grid
    [Ind,Dist,wi] = grid_lookup_20230113(KSpaceCoor,imgsize,fovPE,'os',GridOSFactor,'kernelsize',KernelSize,'verbose',verbose);

    % update progressbar
    if progressbar, waitbar((0.5+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',50)); end

    % reconstruct
    if savekspace
        [rawimg,rawkspace] = gridrecon_fa_20230113(KSpaceCoor,rawdata2,imgsize,fovPE,'wi',wi,'Ind',Ind,'Dist',Dist,'os',GridOSFactor,'k',KernelSize,'verbose',verbose);
    else
        rawimg = gridrecon_fa_20230113(KSpaceCoor,rawdata2,imgsize,fovPE,'wi',wi,'Ind',Ind,'Dist',Dist,'os',GridOSFactor,'k',KernelSize,'verbose',verbose);
    end

    % update progressbar
    if progressbar, waitbar((0.75+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',75)); end

    % combine coils
    if savekspace
        [img,kspace] = combinecoils_fa(rawimg,rawkspace);
    else
        img = combinecoils_fa(rawimg);
    end

    % update progressbar
    if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end  

elseif contains(protocol,'fancy')

    if strcmpi(acqOrder,'ADC')

        % create cartesian grid
        [x,y,z] = ndgrid(-ceil(zfill*imgsize)/2/fovPE:1/(fovPE*GridOSFactor):ceil(zfill*imgsize)/2/fovPE);

        % find N^2 closest grid points for each sampled kspace point
        try
            [Ind,Dist] = knnsearch([x(:) y(:) z(:)],KSpaceCoor,'K',KernelSize^2,'IncludeTies',true);
            Ind = cell2mat(Ind);
            Dist = cell2mat(Dist);
        catch
            [Ind,Dist] = knnsearch([x(:) y(:) z(:)],KSpaceCoor,'K',KernelSize^2);
        end

        % create Kaiser-Bessel Kernel
        beta = pi*sqrt(KernelSize^2/GridOSFactor^2*(GridOSFactor-0.5)^2-0.8); % Eq(5) in https://ieeexplore.ieee.org/document/1435541/
        width = 1*KernelSize/(fovPE*GridOSFactor);
        klength = 10000;
        [kernel,u] = createKBkernel(width,beta,klength);
        kernel3Dtable = interp1(u',kernel,Dist,'linear',0);

        % Calculate density correction (if not provided)
        if isempty(wi)
            iter = 5;
            dcftable = kernel3Dtable.^(1/2);
            wi = iterative_dcf_fa_20190910(iter,KSpaceCoor,dcftable,Ind,size(x,1),[],0);
        end


        % reshape and split acqs with and without diffusion gradients
        rawdata = permute(reshape(rawdata,nsamples,nleaves,nchannels,nrepetitions),[1 2 4 3]);
        rawdata = reshape(rawdata,nsamples,nleaves*nrepetitions,nchannels);
        rawdata1 = reshape(rawdata(:,1:2:end,:),nsamples,nleaves,nrepetitions/2,nchannels);
        rawdata2 = reshape(rawdata(:,2:2:end,:),nsamples,nleaves,nrepetitions/2,nchannels);
        
        % reshape raw data for recon
        rawdata1 = permute(reshape(rawdata1,nsamples*nleaves,nrepetitions/2,nchannels),[1 3 2]);
        rawdata2 = permute(reshape(rawdata2,nsamples*nleaves,nrepetitions/2,nchannels),[1 3 2]);

        % get indices for different rotations
        kinds = 1:nsamples*nleaves:nsamples*nleaves*nreps+1;
        ninds = mod(0:nrepetitions/2-1,nreps)+1;
        
        % initialize outputs
        rawimg1 = zeros(imgsize,imgsize,imgsize,nchannels,nrepetitions/2);
        rawimg2 = zeros(imgsize,imgsize,imgsize,nchannels,nrepetitions/2);
        
        % reconstruct
        for i = 1:nreps
            rawimg1(:,:,:,:,ninds==i) = gridrecon_fa_20220506(KSpaceCoor(kinds(i):kinds(i+1)-1,:),...
                rawdata1(:,:,ninds==i),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi(kinds(i):kinds(i+1)-1,:),verbose,...
                Ind(kinds(i):kinds(i+1)-1,:),Dist(kinds(i):kinds(i+1)-1,:));
            
            rawimg2(:,:,:,:,ninds==i) = gridrecon_fa_20220506(KSpaceCoor(kinds(i):kinds(i+1)-1,:),...
                rawdata2(:,:,ninds==i),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi(kinds(i):kinds(i+1)-1,:),verbose,...
                Ind(kinds(i):kinds(i+1)-1,:),Dist(kinds(i):kinds(i+1)-1,:));
            
            % update progressbar
            if progressbar, waitbar(((i/nreps*0.9)+0.05+enum-1)/N,h,...
                    sprintf('Reconstructing...%3.0f%%',(95*i/nreps))); end
        end
         
        % combine channels
        img1 = combinecoils_fa(rawimg1);
        img2 = combinecoils_fa(rawimg2);
        
        % update progressbar
        if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
        
        % concatenate image sets
        img = cat(4,img1,img2);
        if saveraw
            rawimg = cat(5,rawimg1,rawimg2);
        else 
            rawimg = [];
        end  
        
    else
        
        % reshape raw data
        rawdata = permute(reshape(rawdata,nsamples,nleaves,nchannels,nrepetitions),[1 2 4 3]);

        % get indices for different rotations
        ninds = mod(0:nrepetitions-1,nreps)+1;

        % initialize output
        rawdata2 = zeros(nsamples,nleaves*nreps,nchannels);

        % reconstruct
        for i = 1:nreps
            rawdata2(:,1+(i-1)*nleaves:i*nleaves,:) = squeeze(mean(rawdata(:,:,ninds==i,:),3));
        end
        
        % reshape output
        rawdata2 = reshape(rawdata2,nsamples*nleaves*nreps,nchannels);

        % recon
        if savekspace
            [rawimg, rawkspace] = gridrecon_fa_20220520(KSpaceCoor,rawdata2,imgsize,fovPE,...
                'zfill',zfill,'os',GridOSFactor,'k',KernelSize,'wi',wi,'verbose',verbose);
        else
            rawimg = gridrecon_fa_20220520(KSpaceCoor,rawdata2,imgsize,fovPE,...
                'zfill',zfill,'os',GridOSFactor,'k',KernelSize,'wi',wi,'verbose',verbose);
        end

        % combine channels
        if savekspace
            [img,kspace] = combinecoils_fa(rawimg,rawkspace);
        else
            img = combinecoils_fa(rawimg);
        end
        
        % update progressbar
        if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end       
    end
    
else
    % recon based on acqOrder
    if numRes > 1
        if strcmpi(acqOrder,'Inter-Interleaved')
            
            % reshape raw data (assume only 1 echo)
            rawdata = reshape(rawdata,nsamples,nleaves,max([nslices,npartitions]),nchannels,nrepetitions);
            
            % recon based on correct protocol version
            if contains(protocol,'fa_spiral_dyn_20210622')
                
                % initialize
                KSpaceCoor_gp = [];
                KSpaceCoor_dp = [];
                wi_gp = [];
                wi_dp = [];
                rawdata_gp = [];
                rawdata_dp = [];
                cgp = 1;
                cdp = 1;
                
                % separate gp and dp
                for i = 1:nleaves
                    ind = 1+(i-1)*nsamples:i*nsamples;
                    if mod(i,2) == 1
                        KSpaceCoor_gp = [KSpaceCoor_gp; KSpaceCoor(ind,:)];
                        wi_gp = [wi_gp; wi(ind)];
                        rawdata_gp(:,cgp,:,:,:) = rawdata(:,i,:,:,:);
                        cgp = cgp + 1;
                    else
                        KSpaceCoor_dp = [KSpaceCoor_dp; KSpaceCoor(ind,:)];
                        wi_dp = [wi_dp; wi(ind)];
                        rawdata_dp(:,cdp,:,:,:) = rawdata(:,i,:,:,:);
                        cdp = cdp + 1;
                    end
                end
                
                % reshape (assume only 1 echo)
                rawdata_gp = reshape(rawdata_gp,nsamples*(cgp-1),max([nslices,npartitions]),nchannels,nrepetitions);
                rawdata_dp = reshape(rawdata_dp,nsamples*(cdp-1),max([nslices,npartitions]),nchannels,nrepetitions);
                
                % reconstruct images in batches
                rawimg_gp = [];
                rawimg_dp = [];
                rawkspace_gp = [];
                rawkspace_dp = [];
                for i = 1:ceil(100/(zfill*2)):size(rawdata_gp,4)
                    ind = i:i+ceil(100/(zfill*2))-1;
                    ind(ind > size(rawdata_gp,4)) = [];
                    [temp_gp,tempk_gp] = gridrecon_fa_20210701(KSpaceCoor_gp,rawdata_gp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi_gp,verbose);
                    [temp_dp,tempk_dp] = gridrecon_fa_20210701(KSpaceCoor_dp,rawdata_dp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi_dp,verbose);
                    
                    % 3D FFT
                    if npartitions > 1
                        temp_gp = fftshift(fft(temp_gp,[],3),3);
                        temp_dp = fftshift(fft(temp_dp,[],3),3);
                    end
                    
                    % concatenate
                    rawimg_gp = cat(5, rawimg_gp, temp_gp);
                    rawimg_dp = cat(5, rawimg_dp, temp_dp);
                    rawkspace_gp = cat(5, rawkspace_gp, tempk_gp);
                    rawkspace_dp = cat(5, rawkspace_dp, tempk_dp);
                    
                    % update progressbar
                    if progressbar, waitbar(((ind(end)/size(rawdata_gp,4)*0.9)+0.05+enum-1)/N,h,...
                            sprintf('Reconstructing...%3.0f%%',(95*ind(end)/size(rawdata_gp,4)))); end
                end
                
                % combine channels
                img_gp = combinecoils_fa(rawimg_gp);
                img_dp = combinecoils_fa(rawimg_dp);
                
                % update progressbar
                if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
                
                % concatenate interleaved image sets
                img = cat(4,img_gp,img_dp);
                kspace = cat(5,rawkspace_gp,rawkspace_dp);
                rawimg = cat(5,rawimg_gp,rawimg_dp);
                
            else
                %         elseif  contains(protocol,'fa_spiral_dyn_20210624') ||...
                %                 contains(protocol,'fa_spiral_dyn_20210628') ||...
                %                 contains(protocol,'fa_spiral_dyn_20210709')
                
                % initialize
                rawdata_gp = zeros(size(rawdata));
                rawdata_dp = zeros(size(rawdata));
                
                % separate gp and dp
                if contains(protocol,'fa_spiral_dyn_20210624')
                    c = 0;
                else
                    c = 1;
                end
                for i = 1:nrepetitions
                    for j = 1:max([nslices,npartitions])
                        for k = 1:nleaves
                            if mod(c,2) == 1
                                rawdata_gp(:,k,j,:,i) = rawdata(:,k,j,:,i);
                            else
                                rawdata_dp(:,k,j,:,i) = rawdata(:,k,j,:,i);
                            end
                            c = c + 1;
                        end
                    end
                end
                
                % reshape (assume only 1 echo)
                rawdata_gp = reshape(rawdata_gp,nsamples*nleaves,max([nslices,npartitions]),nchannels,nrepetitions);
                rawdata_dp = reshape(rawdata_dp,nsamples*nleaves,max([nslices,npartitions]),nchannels,nrepetitions);
                
                % reconstruct images in batches
                rawimg_gp = [];
                rawimg_dp = [];
                rawkspace_gp = [];
                rawkspace_dp = [];
                
                % rotate by golden angle
                if GA_flag == 1
                    
                    % initialize
                    rawimg = [];
                    kspace = [];
                    
                    % reshape
                    rawdata = reshape(rawdata,nsamples*nleaves,max([nslices,npartitions]),nchannels,nrepetitions);
                    
                    % pseudo golden angle
                    GA = 2*pi * 2/(sqrt(5)+1);
                    
                    % rotate each repetition
                    for j = 1:nrepetitions
                        rotphiGA = - mod(j-1,GAperiod)*GA;
                        
                        RGA = [cos(rotphiGA) -sin(rotphiGA);
                            sin(rotphiGA)  cos(rotphiGA)];
                        
                        g = RGA^-1 * [KSpaceCoor(:,1)'; KSpaceCoor(:,2)'];
                        
                        [temp,tempk] = gridrecon_fa_20210701(g',rawdata(:,:,:,j),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                        
                        % 3D FFT
                        if npartitions > 1
                            temp = fftshift(fft(temp,[],3),3);
                        end
                        
                        % concatenate
                        rawimg = cat(5, rawimg, temp);
                        kspace = cat(5,kspace,tempk);
                        
                        % update progressbar
                        if progressbar, waitbar(((j/size(rawdata,4)*0.9)+0.05+enum-1)/N,h,...
                                sprintf('Reconstructing...%3.0f%%',(95*j/size(rawdata,4)))); end
                    end
                    
                    % combine channels
                    img = combinecoils_fa(rawimg);
                    
                    % update progressbar
                    if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
                    
                else
                    for i = 1:ceil(100/(zfill*2)):size(rawdata_gp,4)
                        ind = i:i+ceil(100/(zfill*2))-1;
                        ind(ind > size(rawdata_gp,4)) = [];
                        if savekspace
                            [temp_gp,tempk_gp] = gridrecon_fa_20210701(KSpaceCoor,rawdata_gp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                            [temp_dp,tempk_dp] = gridrecon_fa_20210701(KSpaceCoor,rawdata_dp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                        else
                            temp_gp = gridrecon_fa_20210701(KSpaceCoor,rawdata_gp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                            temp_dp = gridrecon_fa_20210701(KSpaceCoor,rawdata_dp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                            tempk_gp = [];
                            tempk_dp = [];
                        end
                        
                        % 3D FFT
                        if npartitions > 1
                            temp_gp = fftshift(fft(temp_gp,[],3),3);
                            temp_dp = fftshift(fft(temp_dp,[],3),3);
                        end
                        
                        % concatenate
                        rawimg_gp = cat(5, rawimg_gp, temp_gp);
                        rawimg_dp = cat(5, rawimg_dp, temp_dp);
                        rawkspace_gp = cat(5, rawkspace_gp, tempk_gp);
                        rawkspace_dp = cat(5, rawkspace_dp, tempk_dp);
                        
                        % update progressbar
                        if progressbar, waitbar(((ind(end)/size(rawdata_gp,4)*0.9)+0.05+enum-1)/N,h,...
                                sprintf('Reconstructing...%3.0f%%',(95*ind(end)/size(rawdata_gp,4)))); end
                    end
                    
                    % combine channels
                    img_gp = combinecoils_fa(rawimg_gp);
                    img_dp = combinecoils_fa(rawimg_dp);
                    
                    % update progressbar
                    if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
                    
                    % concatenate interleaved image sets
                    img = cat(4,img_gp,img_dp);
                    kspace = cat(5,rawkspace_gp,rawkspace_dp);
                    rawimg = cat(5,rawimg_gp,rawimg_dp);
                    
                end
            end
            
        elseif strcmpi(acqOrder,'Radial-Interleaved')
            
            % get raw echo data
            rawechoes = reco_params.VD13_Data(:,1:end-1);
            
            % reshape raw data (assume 2 echoes)
            rawdata = reshape(rawdata,nsamples,nleaves,max([nslices,npartitions]),nchannels,2,nrepetitions);
            
            if contains(protocol,'fa_spiral_dyn_20210722')
                % get params
                freq = [];%acq_params.frequency/1.498*1e-6;
                orient = acq_params.orientation;
                c = 0;
                ncols = 64;
                ncolskip = 0;
                ncols2 = ncols - ncolskip;
                
                % initialize
                rawdata_gp = zeros(nsamples,nleaves,max([nslices,npartitions]),nchannels,nrepetitions);
                raw_out = zeros(ncols2,nleaves,nrepetitions,max([nslices,npartitions]),nchannels);
                raw_rl = zeros(ncols*2,nleaves,nrepetitions,max([nslices,npartitions]),nchannels);
                raw_lr = zeros(ncols*2,nleaves,nrepetitions,max([nslices,npartitions]),nchannels);
                kspace_dp_out = zeros(2,ncols2,nleaves,nrepetitions,max([nslices,npartitions]));
                kspace_dp_rl = zeros(2,ncols*2,nleaves,nrepetitions,max([nslices,npartitions]));
                kspace_dp_lr = zeros(2,ncols*2,nleaves,nrepetitions,max([nslices,npartitions]));
                
                % loop through raw data
                for i = 1:nrepetitions
                    for n = 1:max([nslices,npartitions])
                        for k = 1:nleaves
                            if mod(c,2)~=0
                                for j = 1:nchannels
                                    ind = (i-1)*max([nslices,npartitions])*nleaves*nchannels+(n-1)*nleaves*nchannels+(k-1)*nchannels+j;
                                    if ind <= size(rawechoes,2)
                                        raw_out(:,k,i,n,j) = rawechoes(1+ncolskip+1:ncols+1,ind);
                                        raw_rl(:,k,i,n,j) = rawechoes(1+ncols+1:ncols*3+1,ind);
                                        raw_lr(:,k,i,n,j) = rawechoes(1+ncols*3+1:ncols*5+1,ind);
                                    end
                                end
                                spoke_out = calcradialspoke(ncols,c,freq,'spoke',0,orient);
                                kspace_dp_out(:,:,k,i,n) = spoke_out(1+ncolskip:end,:)';
                                kspace_dp_rl(:,:,k,i,n) = ...
                                    calcradialspoke(ncols*2,c,freq,'traverserl',norm(spoke_out(end,:)),orient)';
                                kspace_dp_lr(:,:,k,i,n) = ...
                                    calcradialspoke(ncols*2,c,freq,'traverselr',norm(spoke_out(end,:)),orient)';
                            else
                                rawdata_gp(:,k,n,:,i) = rawdata(:,k,n,:,1,i);
                            end
                            c = c + 1;
                        end
                    end
                end
                
            else %if contains(protocol,'fa_spiral_dyn_20210809')
                % get params
                freq = [];%acq_params.frequency/1.498*1e-6;
                orient = acq_params.orientation;
                c = 0;
                ncols = 64;
                ncolskip = 7;
                ncolshift = 8;
                ncols2 = 66;
                
                % initialize
                rawdata_gp = zeros(nsamples,nleaves,max([nslices,npartitions]),nchannels,nrepetitions);
                raw_out = zeros(ncols2,nleaves,nrepetitions,max([nslices,npartitions]),nchannels);
                raw_rl = zeros(ncols*2,nleaves,nrepetitions,max([nslices,npartitions]),nchannels);
                raw_lr = zeros(ncols*2,nleaves,nrepetitions,max([nslices,npartitions]),nchannels);
                kspace_dp_out = zeros(2,ncols2,nleaves,nrepetitions,max([nslices,npartitions]));
                kspace_dp_rl = zeros(2,ncols*2,nleaves,nrepetitions,max([nslices,npartitions]));
                kspace_dp_lr = zeros(2,ncols*2,nleaves,nrepetitions,max([nslices,npartitions]));
                
                % loop through raw data
                for i = 1:nrepetitions
                    for n = 1:max([nslices,npartitions])
                        for k = 1:nleaves
                            if mod(c,2)~=0
                                for j = 1:nchannels
                                    ind = (i-1)*max([nslices,npartitions])*nleaves*nchannels+(n-1)*nleaves*nchannels+(k-1)*nchannels+j;
                                    if ind <= size(rawechoes,2)
                                        raw_out(:,k,i,n,j) = rawechoes(1+ncolskip:ncols2+ncolskip,ind);
                                        raw_rl(:,k,i,n,j) = rawechoes(1+ncols+ncolshift:ncols*3+ncolshift,ind);
                                        raw_lr(:,k,i,n,j) = rawechoes(1+ncols*3+ncolshift:ncols*5+ncolshift,ind);
                                    end
                                end
                                spoke_out = calcradialspoke_20210809(ncols,c-1,freq,'spoke',-6e-4,orient);
                                kspace_dp_out(:,:,k,i,n) = spoke_out(1:end,:)';
                                kspace_dp_rl(:,:,k,i,n) = ...
                                    calcradialspoke_20210809(ncols*2,c-1,freq,'traverserl',norm(spoke_out(end,:))+6e-4,orient)';
                                kspace_dp_lr(:,:,k,i,n) = ...
                                    calcradialspoke_20210809(ncols*2,c-1,freq,'traverselr',norm(spoke_out(end,:))+6e-4,orient)';
                            else
                                rawdata_gp(:,k,n,:,i) = rawdata(:,k,n,:,1,i);
                            end
                            c = c + 1;
                        end
                    end
                end
            end
            
            % reformat echo data
            raw_out = reshape(raw_out,ncols2*nleaves*nrepetitions,max([nslices,npartitions]),nchannels);
            kspace_dp_out = permute(reshape(kspace_dp_out,2,ncols2*nleaves*nrepetitions,max([nslices,npartitions])),[2 1 3]);
            raw_rl = reshape(raw_rl,ncols*2*nleaves*nrepetitions,max([nslices,npartitions]),nchannels);
            kspace_dp_rl = permute(reshape(kspace_dp_rl,2,ncols*2*nleaves*nrepetitions,max([nslices,npartitions])),[2 1 3]);
            raw_lr = reshape(raw_lr,ncols*2*nleaves*nrepetitions,max([nslices,npartitions]),nchannels);
            kspace_dp_lr = permute(reshape(kspace_dp_lr,2,ncols*2*nleaves*nrepetitions,max([nslices,npartitions])),[2 1 3]);
            
            % initialize images
            fovPE = fovPE * 1.4;
            rawimg_dp = zeros(ncols*zfill,ncols*zfill,max([nslices,npartitions]),nchannels,3);
            if savekspace
                klength = length(-ceil(zfill*ncols)/2/fovPE:1/(fovPE*GridOSFactor):ceil(zfill*ncols)/2/fovPE);
                rawkspace_dp = zeros(klength,klength,max([nslices,npartitions]),nchannels,3);
            else
                rawkspace_dp = [];
            end
            
            % reconstruct DP echoes
            for i = 1:max([nslices,npartitions])
                if savekspace
                    tmpind = find(raw_out(:,i,1)~=0);
                    [rawimg_dp(:,:,i,:,1),rawkspace_dp(:,:,i,:,1)] = gridrecon_fa_20210701(kspace_dp_out(tmpind,:,i),raw_out(tmpind,i,:),ncols,fovPE,zfill,GridOSFactor,KernelSize,[],verbose);
                    tmpind = find(raw_rl(:,i,1)~=0);
                    [rawimg_dp(:,:,i,:,2),rawkspace_dp(:,:,i,:,2)] = gridrecon_fa_20210701(kspace_dp_rl (tmpind,:,i),raw_rl (tmpind,i,:),ncols,fovPE,zfill,GridOSFactor,KernelSize,[],verbose);
                    tmpind = find(raw_lr(:,i,1)~=0);
                    [rawimg_dp(:,:,i,:,3),rawkspace_dp(:,:,i,:,3)] = gridrecon_fa_20210701(kspace_dp_lr (tmpind,:,i),raw_lr (tmpind,i,:),ncols,fovPE,zfill,GridOSFactor,KernelSize,[],verbose);
                else
                    tmpind = find(raw_out(:,i,1)~=0);
                    rawimg_dp(:,:,i,:,1) = gridrecon_fa_20210701(kspace_dp_out(tmpind,:,i),raw_out(tmpind,i,:),ncols,fovPE,zfill,GridOSFactor,KernelSize,[],verbose);
                    tmpind = find(raw_rl(:,i,1)~=0);
                    rawimg_dp(:,:,i,:,2) = gridrecon_fa_20210701(kspace_dp_rl (tmpind,:,i),raw_rl (tmpind,i,:),ncols,fovPE,zfill,GridOSFactor,KernelSize,[],verbose);
                    tmpind = find(raw_lr(:,i,1)~=0);
                    rawimg_dp(:,:,i,:,3) = gridrecon_fa_20210701(kspace_dp_lr (tmpind,:,i),raw_lr (tmpind,i,:),ncols,fovPE,zfill,GridOSFactor,KernelSize,[],verbose);
                end
                % update progressbar
                if progressbar, waitbar(((i/max([nslices,npartitions])*0.45)+0.05+enum-1)/N,h,...
                        sprintf('Reconstructing...%3.0f%%',(i/max([nslices,npartitions])*50))); end
                
            end
            
            fovPE = fovPE / 1.4;
            
            % 3D FFT
            if npartitions > 1
                rawimg_dp = fftshift(fft(rawimg_dp,[],3),3);
            end
            
            % reshape GP data
            rawdata_gp = reshape(rawdata_gp,nsamples*nleaves,max([nslices,npartitions]),nchannels,nrepetitions);
            rawimg_gp = [];
            rawkspace_gp = [];
            
            % reconstruct GP images in batches
            for i = 1:ceil(100/(zfill*2)):size(rawdata_gp,4)
                ind = i:i+ceil(100/(zfill*2))-1;
                ind(ind > size(rawdata_gp,4)) = [];
                if savekspace
                    [temp_gp,tempk_gp] = gridrecon_fa_20210701(KSpaceCoor,rawdata_gp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                else
                    temp_gp = gridrecon_fa_20210701(KSpaceCoor,rawdata_gp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                    tempk_gp = [];
                end
                
                % 3D FFT
                if npartitions > 1
                    temp_gp = fftshift(fft(temp_gp,[],3),3);
                end
                
                % concatenate
                rawimg_gp = cat(5, rawimg_gp, temp_gp);
                rawkspace_gp = cat(5, rawkspace_gp, tempk_gp);
                
                % update progressbar
                if progressbar, waitbar(((ind(end)/size(rawdata_gp,4)*0.45)+0.5+enum-1)/N,h,...
                        sprintf('Reconstructing...%3.0f%%',(50+45*ind(end)/size(rawdata_gp,4)))); end
            end
            
            % interpolate to same size
            if size(rawimg_dp,1) ~= size(rawimg_gp,1)
                rawimg_dp = resize_fa(rawimg_dp,[size(rawimg_gp,1) size(rawimg_gp,2)],'nearest');
            end
            if savekspace && (size(rawkspace_dp,1) ~= size(rawkspace_gp,1))
                rawkspace_dp = resize_fa(rawkspace_dp,[size(rawkspace_gp,1) size(rawkspace_gp,2)],'nearest');
            end
            
            % combine channels
            img_gp = combinecoils_fa(rawimg_gp);
            img_dp = zeros([size(img_gp,1:3),size(rawimg_dp,5)]);
            for i = 1:size(rawimg_dp,5)
                img_dp(:,:,:,i) = combinecoils_fa(rawimg_dp(:,:,:,:,i));
            end
            
            % update progressbar
            if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
            
            % concatenate interleaved image sets
            img = cat(4,img_gp,img_dp);
            kspace = cat(5,rawkspace_gp,rawkspace_dp);
            rawimg = cat(5,rawimg_gp,rawimg_dp);
            
        elseif strcmpi(acqOrder,'Sequential')
            
            % separate resonances
            rawdata_gp = rawdata(:,:,:,1:nrepetitions/2);
            rawdata_dp = rawdata(:,:,:,1+nrepetitions/2:end);
            
            % reconstruct images in batches
            rawimg_gp = [];
            rawimg_dp = [];
            rawkspace_gp = [];
            rawkspace_dp = [];
            for i = 1:ceil(100/(zfill*2)):size(rawdata_gp,4)
                ind = i:i+ceil(100/(zfill*2))-1;
                ind(ind > size(rawdata_gp,4)) = [];
                [temp_gp,tempk_gp] = gridrecon_fa_20210701(KSpaceCoor,rawdata_gp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                [temp_dp,tempk_dp] = gridrecon_fa_20210701(KSpaceCoor,rawdata_dp(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                
                % 3D FFT
                if npartitions > 1
                    temp_gp = fftshift(fft(temp_gp,[],3),3);
                    temp_dp = fftshift(fft(temp_dp,[],3),3);
                end
                
                % concatenate
                rawimg_gp = cat(5, rawimg_gp, temp_gp);
                rawimg_dp = cat(5, rawimg_dp, temp_dp);
                rawkspace_gp = cat(5, rawkspace_gp, tempk_gp);
                rawkspace_dp = cat(5, rawkspace_dp, tempk_dp);
                
                % update progressbar
                if progressbar, waitbar(((ind(end)/size(rawdata_gp,4)*0.9)+0.05+enum-1)/N,h,...
                        sprintf('Reconstructing...%3.0f%%',(95*ind(end)/size(rawdata_gp,4)))); end
            end
            
            % combine channels
            img_gp = combinecoils_fa(rawimg_gp);
            img_dp = combinecoils_fa(rawimg_dp);
            
            % update progressbar
            if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
            
            % concatenate interleaved image sets
            img = cat(4,img_gp,img_dp);
            kspace = cat(5,rawkspace_gp,rawkspace_dp);
            rawimg = cat(5,rawimg_gp,rawimg_dp);
            
        elseif strcmpi(acqOrder,'Saturation')
            % reconstruct images in batches
            rawimg = [];
            kspace = [];
            for i = 1:ceil(100/(zfill*2)):size(rawdata,4)
                ind = i:i+ceil(100/(zfill*2))-1;
                ind(ind > size(rawdata,4)) = [];
                if savekspace
                    [temp,tempk] = gridrecon_fa_20210701(KSpaceCoor,rawdata(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                else
                    temp = gridrecon_fa_20210701(KSpaceCoor,rawdata(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                    tempk = [];
                end
                
                % 3D FFT
                if npartitions > 1
                    temp = fftshift(fft(temp,[],3),3);
                end
                
                % concatenate
                rawimg = cat(5, rawimg, temp);
                kspace = cat(5, kspace, tempk);
                
                % update progressbar
                if progressbar, waitbar(((ind(end)/size(rawdata,4)*0.9)+0.05+enum-1)/N,h,...
                        sprintf('Reconstructing...%3.0f%%',(95*ind(end)/size(rawdata,4)))); end
            end
            
            % combine channels
            img = combinecoils_fa(rawimg);
            
            % update progressbar
            if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
        end
    else
        % rotate by golden angle
        if GA_flag == 1
            
            % initialize
            rawimg = [];
            kspace = [];
            
            % pseudo golden angle
            GA = 2*pi * 2/(sqrt(5)+1);
            
            % rotate each repetition
            for j = 1:nrepetitions
                rotphiGA = - mod(j-1,GAperiod)*GA;
                
                RGA = [cos(rotphiGA) -sin(rotphiGA);
                    sin(rotphiGA)  cos(rotphiGA)];
                
                g = RGA^-1 * [KSpaceCoor(:,1)'; KSpaceCoor(:,2)'];
                
                [temp,tempk] = gridrecon_fa_20210701(g',rawdata(:,:,:,j),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                
                % 3D FFT
                if npartitions > 1
                    temp = fftshift(fft(temp,[],3),3);
                end
                
                % concatenate
                rawimg = cat(5, rawimg, temp);
                kspace = cat(5, kspace, tempk);
                
                % update progressbar
                if progressbar, waitbar(((j/size(rawdata,4)*0.9)+0.05+enum-1)/N,h,...
                        sprintf('Reconstructing...%3.0f%%',(95*j/size(rawdata,4)))); end
            end
            
            % combine channels
            img = combinecoils_fa(rawimg);
            
            % update progressbar
            if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
            
        else
            % reconstruct images in batches
            rawimg = [];
            kspace = [];
            for i = 1:ceil(100/(zfill*2)):size(rawdata,4)
                ind = i:i+ceil(100/(zfill*2))-1;
                ind(ind > size(rawdata,4)) = [];
                if savekspace
                    [temp,tempk] = gridrecon_fa_20210701(KSpaceCoor,rawdata(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                else
                    temp = gridrecon_fa_20210701(KSpaceCoor,rawdata(:,:,:,ind),imgsize,fovPE,zfill,GridOSFactor,KernelSize,wi,verbose);
                    tempk = [];
                end
                
                % 3D FFT
                if npartitions > 1
                    temp = fftshift(fft(temp,[],3),3);
                end
                
                % concatenate
                rawimg = cat(5, rawimg, temp);
                kspace = cat(5, kspace, tempk);
                
                % update progressbar
                if progressbar, waitbar(((ind(end)/size(rawdata,4)*0.9)+0.05+enum-1)/N,h,...
                        sprintf('Reconstructing...%3.0f%%',(95*ind(end)/size(rawdata,4)))); end
            end
            
            % combine channels
            img = combinecoils_fa(rawimg);
            
            % update progressbar
            if progressbar, waitbar((1.0+enum-1)/N,h,sprintf('Reconstructing...%3.0f%%',100)); end
            
        end
    end
end

end
