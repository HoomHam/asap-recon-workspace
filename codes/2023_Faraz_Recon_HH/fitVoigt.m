function [fitParams, err] = fitVoigt(freq, spectra, options)

% parse inputs
arguments
    freq {mustBeNumeric}
    spectra {mustBeNumeric}
    options.peaks double = 2
    options.peaks_locs double = []
    options.MinPeakDistance double = 2
    options.MinPeakHeight double = 0
    options.lb double = [-Inf, -Inf, 0.1, 0, -pi, -pi, -Inf, -Inf]
    options.ub double = [Inf, Inf, 6, 1, pi, pi, Inf, Inf]
    options.startPoint double = []
    options.warnings double = 0
    options.opts = []
    options.solver = []
    options.view double = 0
    options.normalize double = 1
end

% check solver
if isempty(options.solver)
    if isempty(options.opts)
        options.solver = 'lsqnonlin';
    else
        if strcmpi(class(options.opts),'optim.options.Lsqnonlin')
            options.solver = 'lsqnonlin';
        elseif strcmpi(class(options.opts),'optim.options.Lsqcurvefit')
            options.solver = 'lsqcurvefit';
        end
    end
end

% define fit options
if isempty(options.opts)
    if strcmpi(options.solver,'lsqnonlin')
        lsq_opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
            'MaxFunctionEvaluations',5e4,'MaxIterations',5e4,'Display','off',...
            'FunctionTolerance',1e-7,'StepTolerance',1e-8);
    elseif strcmpi(options.solver,'lsqcurvefit')
        lsq_opts = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt',...
            'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'Display','off',...
            'FunctionTolerance',1e-9,'StepTolerance',1e-8);
    end
else
    lsq_opts = options.opts;
end

% suppress warnings
if options.warnings == 0
    warning('off','optimlib:levenbergMarquardt:InfeasibleX0');
end

% create figure 
if options.view == 1
    f1 = figure;
end

% reshape
freq = reshape(freq,[],1);

% check if double
spectra = double(spectra);

% normalize
if options.normalize == 1
    scale = 10 / max(abs(spectra));
    spectra = spectra .* scale;
else
    scale = 1;
end

% figure(300)
% plot(abs(spectra))
% beep
% pause

% initialize
fitParams = [];
err = [];

% find peaks of each resonance
if isempty(options.peaks_locs)
    if options.peaks == 2
        MaxAmps = zeros(2,1);
        MaxInds = zeros(2,1);
        [MaxAmps(1), MaxInds(1)] = findpeaks(abs(spectra),'MinPeakDistance',options.MinPeakDistance,...
            'MinPeakHeight',options.MinPeakHeight,...
            'NPeaks',1,'SortStr','descend');
        [MaxAmps(2), MaxInds(2)] = findpeaks(abs(spectra(MaxInds(1)+1:end)),'MinPeakDistance',options.MinPeakDistance,...
            'MinPeakHeight',options.MinPeakHeight,...
            'NPeaks',1,'SortStr','descend');
        MaxInds(2) = MaxInds(2) + MaxInds(1);
    else
        [MaxAmps, MaxInds] = findpeaks(abs(spectra),'MinPeakDistance',options.MinPeakDistance,...
            'MinPeakHeight',options.MinPeakHeight,...
            'NPeaks',options.peaks,'SortStr','descend');
    end
else
    MaxInds = options.peaks_locs;
    MaxAmps = abs(spectra(MaxInds));
    options.peaks = length(MaxInds);
end

% define start point
if isempty(options.startPoint)
    options.startPoint = zeros(options.peaks,8);
    for j = 1:options.peaks
        % find signal less than half the maximum value
        %     LowInd = find(abs(spectra(MaxInds(j):end)) <= abs(MaxAmps(j)/2));
        LowInds = abs(spectra) <= abs(MaxAmps(j)/2);

        %     StartFWHM = abs(freq(MaxInds(j)) - freq(MaxInds(j)+LowInd(1)-1));
        StartFWHM = min(abs(freq(MaxInds(j))-freq(LowInds)));
        
        if isempty(StartFWHM)
            StartFWHM = NaN;
        end

        % define starting point
        options.startPoint(j,:) = [MaxAmps(j), freq(MaxInds(j)), StartFWHM, 0.5, 0, 0, 0, 0];
    end
end

% check bounds

if sum(isnan(options.startPoint(:,3))),
    options.startPoint(:,3)
    
    err = NaN;
    fitParams(:,9) = NaN;

else

    if size(options.lb,1) == 1 && options.peaks > 1
        options.lb = repmat(options.lb,options.peaks,1);
    end
    if size(options.ub,1) == 1 && options.peaks > 1
        options.ub = repmat(options.ub,options.peaks,1);
    end

    % fit
    if strcmpi(options.solver,'lsqnonlin')
        [x, y] = lsqnonlin(@voigt, options.startPoint,...
            options.lb, options.ub, lsq_opts);
    elseif strcmpi(options.solver,'lsqcurvefit')
        [x, y] = lsqcurvefit(@(params,x)[real(voigt2(params,x)); imag(voigt2(params,x))], ...
            options.startPoint, [freq; freq], ...
            [real(spectra); imag(spectra)],...
            options.lb, options.ub, lsq_opts);
    elseif strcmpi(options.solver,'fminsearch')
        [x, y] = lsqnonlin(@voigt, options.startPoint);
    end

    fitParams = cat(3, fitParams, x);
    err = [err; y];
    fitParams(:,9) = scale;

end


% % close figure 
% if options.view == 1
%     close(f1);
% end

% un-suppress warnings
if options.warnings == 0
    warning('on','optimlib:levenbergMarquardt:InfeasibleX0');
end

    function y = voigt2(params,x)
        voigt = 0;
        for i = 1:size(params,1)
            amp = params(i,1);
            shift = params(i,2);
            fwhmL = params(i,3);
            fwhmG = params(i,3);
            eta = params(i,4);
            phase = params(i,5);
            phaserec = params(i,6);
            offsetR = params(i,7);
            offsetI = params(i,8);
            
            df = x(1:end/2) - shift;
            offset = complex(offsetR,offsetI);

            lorR = fwhmL ./ (df.^2 + fwhmL.^2);
            lorI = -df ./ (df.^2 + fwhmL.^2);
            gauss = exp(-(df .* sqrt(log(2)) / fwhmG).^2);

            voigt_tmp = amp .* exp(1i*phase) .* ...
                (eta * (complex(lorR,0) + complex(0,lorI) * exp(1i*phaserec)) +...
                (1 - eta) * gauss) + offset;
            
            voigt = voigt + voigt_tmp;

        end
        y = voigt;
    end

    function err = voigt(params)
        voigt = 0;
        for i = 1:size(params,1)
            amp = params(i,1);
            shift = params(i,2);
            fwhmL = params(i,3);
            fwhmG = params(i,3);
            eta = params(i,4);
            phase = params(i,5);
            phaserec = params(i,6);
            offsetR = params(i,7);
            offsetI = params(i,8);
            
            df = freq - shift;
            offset = complex(offsetR,offsetI);

            lorR = fwhmL ./ (df.^2 + fwhmL.^2);
            lorI = -df ./ (df.^2 + fwhmL.^2);
            gauss = exp(-(df .* sqrt(log(2)) / fwhmG).^2);

            voigt_tmp = amp .* exp(1i*phase) .* ...
                (eta * (complex(lorR,0) + complex(0,lorI) * exp(1i*phaserec)) +...
                (1 - eta) * gauss) + offset;
            
            voigt = voigt + voigt_tmp;

        end
        err = sum(real(voigt-spectra).^2) + sum(imag(voigt-spectra).^2);
        
        if exist('f1','var')
            if randi(10000) == 25
                figure(f1);
                subplot(2,2,1);
                plot(abs(spectra)); hold on;
                plot(abs(voigt)); hold off;
                subplot(2,2,2);
                plot(real(spectra)); hold on;
                plot(real(voigt)); hold off;
                subplot(2,2,3);
                plot(imag(spectra)); hold on;
                plot(imag(voigt)); hold off;
                subplot(2,2,4);
                hold on;
                if isempty(f1.Children(1).Children)
                    scatter(1,err);
                else
                    scatter(f1.Children(1).Children(1).XData+1,err);
                end
                hold off;
            end
        end
    end

end
