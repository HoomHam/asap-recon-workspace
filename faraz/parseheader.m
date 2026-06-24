function special = parseheader(twix,protocol)
% Parses through the variables of the "Special" tab for Siemens scanners

% check input
if nargin < 2, protocol = ''; end

% initalize output structure
special = struct;

if isfield(twix.hdr.MeasYaps,'sWipMemBlock')
    if contains(protocol,'fa_spiral_dyn')
        if contains(protocol,'fancy_v2') || contains(protocol,'fancy_v3') || contains(protocol,'fancy_v4')
            % version number
            if contains(protocol,'fancy_v4')
                            version = 20250815
            else
            version = extract(protocol,digitsPattern(8));
            version = str2double(version{1});
            end

            % tx factor
            special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};

            % number of images
            special.numImg = twix.hdr.MeasYaps.sWipMemBlock.alFree{2};

            % number of resonances
            special.numRes = 1 + twix.hdr.MeasYaps.sWipMemBlock.alFree{9};
            if isempty(special.numRes), special.numRes = 1; end

            % pulse shape for additional resonances
            special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
            if isempty(special.RFshape)
                special.RFshape = 'aSinc';
            elseif special.RFshape == 1
                special.RFshape = 'Gauss';
            elseif special.RFshape == 2
                special.RFshape = 'Rect';
            elseif special.RFshape == 3
                special.RFshape = 'Kai';
            end

            % acquisition order
            special.acqOrder = twix.hdr.MeasYaps.sWipMemBlock.alFree{10};
            if isempty(special.acqOrder)
                special.acqOrder = 'Sequential';
            elseif special.acqOrder == 1
                special.acqOrder = 'GP-DP';
            elseif special.acqOrder == 2
                special.acqOrder = 'GP-DP-DP';
            else
                special.acqOrder = 'who fucking knows';
            end

            % frequency, flip angle, and pulse duration for other resonances
            special.RFdur = 1.0240; % ms
            if special.numRes >= 2
                special.freq2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{3}; % ppm
                special.FA2 = twix.hdr.MeasYaps.sWipMemBlock.alFree{4}; % deg
                special.RFdur2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{5}; % ms
            end

            % matrix size
            special.MatSize = twix.hdr.MeasYaps.sWipMemBlock.alFree{6};

            % golden angle periodicity
            special.GAperiod = twix.hdr.MeasYaps.sWipMemBlock.alFree{7};

            % check if spectra acquired
            if twix.hdr.MeasYaps.sWipMemBlock.alFree{11} == 1
                special.bSpectra = true;
            else
                special.bSpectra = false;
            end

            % spectra parameters
            if special.bSpectra
                special.numSpec = twix.hdr.MeasYaps.sWipMemBlock.alFree{12}; % number of spectra
                special.dtSpec = twix.hdr.MeasYaps.sWipMemBlock.alFree{13}; % dwelltime (us/point)
                special.FASpec = twix.hdr.MeasYaps.sWipMemBlock.alFree{14}; % deg
            end

        elseif contains(protocol,'fancy')
            % version number
            version = extract(protocol,digitsPattern(8));
            version = str2double(version{1});

            % tx factor
            special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};

            % number of images
            special.numImg = twix.hdr.MeasYaps.sWipMemBlock.alFree{2};

            % number of resonances
            special.numRes = 1 + twix.hdr.MeasYaps.sWipMemBlock.alFree{9};
            if isempty(special.numRes), special.numRes = 1; end

            if version < 20220316
                % delay between images (ms)
                special.repDelay = twix.hdr.MeasYaps.sWipMemBlock.adFree{3};
            else
                % number of images before alternating saturation
                if special.numRes > 1
                    special.numAlt = twix.hdr.MeasYaps.sWipMemBlock.alFree{3};
                end
            end

            % pulse shape for additional resonances
            special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
            if isempty(special.RFshape)
                special.RFshape = 'aSinc';
            elseif special.RFshape == 1
                special.RFshape = 'Gauss';
            elseif special.RFshape == 2
                special.RFshape = 'Rect';
            elseif special.RFshape == 3
                special.RFshape = 'Kai';
            end
            
            % acquisition order
            special.acqOrder = twix.hdr.MeasYaps.sWipMemBlock.alFree{10};
            if isempty(special.acqOrder)
                special.acqOrder = 'Sequential';
            elseif special.acqOrder == 1
                special.acqOrder = 'Interleaved';
            elseif special.acqOrder == 2
                special.acqOrder = 'Inter-Interleaved';
            elseif special.acqOrder == 3
                special.acqOrder = 'Saturation';
            elseif special.acqOrder == 4
                special.acqOrder = 'ADC';
                special.bvalue = twix.hdr.MeasYaps.sWipMemBlock.adFree{11};
            elseif special.acqOrder == 5
                special.acqOrder = 'Forced DP Traj';
            elseif special.acqOrder == 6
                special.acqOrder = 'GP-DP-DP';
            elseif special.acqOrder == 7
                special.acqOrder = 'GP-DP-DP-Spectra';
            else
                special.acqOrder = 'who fucking knows';
            end
            
            % frequency, flip angle, and pulse duration for other resonances
            special.RFdur = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms
            if special.numRes >= 2
                special.freq2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{4}; % ppm
                special.FA2 = twix.hdr.MeasYaps.sWipMemBlock.alFree{5}; % deg
                special.RFdur2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms
                
                if special.numRes >= 3
                    special.freq3 = twix.hdr.MeasYaps.sWipMemBlock.adFree{11}; % ppm
                    special.FA3 = twix.hdr.MeasYaps.sWipMemBlock.alFree{12}; % deg
                    special.RFdur3 = twix.hdr.MeasYaps.sWipMemBlock.adFree{13}; % ms
                end
            end
            
            % matrix size 
            special.MatSize = twix.hdr.MeasYaps.sWipMemBlock.alFree{7};
            
            % golden angle periodicity
            special.GAperiod = twix.hdr.MeasYaps.sWipMemBlock.alFree{14};
        else
            % tx factor
            special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};
            
            % number of images
            special.numImg = twix.hdr.MeasYaps.sWipMemBlock.alFree{2};
            
            % number of resonances
            special.numRes = 1 + twix.hdr.MeasYaps.sWipMemBlock.alFree{9};
            if isempty(special.numRes), special.numRes = 1; end
            
            % delay between images (ms)
            special.repDelay = twix.hdr.MeasYaps.sWipMemBlock.adFree{3};
            
            % pulse shape for additional resonances
            special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
            if isempty(special.RFshape)
                special.RFshape = 'aSinc';
            elseif special.RFshape == 1
                special.RFshape = 'Gauss';
            elseif special.RFshape == 2
                special.RFshape = 'Rect';
            elseif special.RFshape == 3
                special.RFshape = 'Kai';
            end
            
            % acquisition order
            special.acqOrder = twix.hdr.MeasYaps.sWipMemBlock.alFree{10};
            if isempty(special.acqOrder)
                special.acqOrder = 'Sequential';
            elseif special.acqOrder == 1
                special.acqOrder = 'Interleaved';
            elseif special.acqOrder == 2
                special.acqOrder = 'Inter-Interleaved';
            elseif special.acqOrder == 3
                special.acqOrder = 'Radial-Interleaved';
            elseif special.acqOrder == 4
                special.acqOrder = 'Saturation';
            end
            
            % frequency, flip angle, and pulse duration for other resonances
            special.RFdur = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms
            if special.numRes >= 2
                special.freq2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{4}; % ppm
                special.FA2 = twix.hdr.MeasYaps.sWipMemBlock.alFree{5}; % deg
                special.RFdur2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms
                
                if special.numRes >= 3
                    special.freq3 = twix.hdr.MeasYaps.sWipMemBlock.adFree{11}; % ppm
                    special.FA3 = twix.hdr.MeasYaps.sWipMemBlock.alFree{12}; % deg
                    special.RFdur3 = twix.hdr.MeasYaps.sWipMemBlock.adFree{13}; % ms
                end
            end
            
            % golden angle trajectory used?
            special.GA = twix.hdr.MeasYaps.sWipMemBlock.alFree{7};
            
            % golden angle periodicity
            if special.GA == 1
                special.GAperiod = twix.hdr.MeasYaps.sWipMemBlock.alFree{14};
            end
        end
        
    elseif contains(protocol,'fa_spiral_hybridXTC')
        
        % tx factor
        special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};
        
        % default saturation method
        special.SatMethod = 'none';
        
        % number of images per set
        special.NumImg = twix.hdr.MeasYaps.sWipMemBlock.alFree{3};
        
        % GP saturation after 1st image?
        if twix.hdr.MeasYaps.sWipMemBlock.alFree{7} == 1
            special.GPsat = true;
        else
            special.GPsat = false;
        end
        
        if ~isempty(twix.hdr.MeasYaps.sWipMemBlock.alFree{2})
            % saturation pulse flip angle
            special.SatFA = twix.hdr.MeasYaps.sWipMemBlock.alFree{8}; % deg
            
            % saturation pulse width
            special.SatDuration = twix.hdr.MeasYaps.sWipMemBlock.adFree{9}; % ms
            
            % saturation pulse type
            if twix.hdr.MeasYaps.sWipMemBlock.alFree{10} == 1
                special.SatPulse = 'Rect';
            else
                special.SatPulse = 'Gauss';
            end
            
            % number frequency, and spacing of RBC saturations
            if ismember(twix.hdr.MeasYaps.sWipMemBlock.alFree{2}, [1 3])
                special.SatNum_RBC = twix.hdr.MeasYaps.sWipMemBlock.alFree{4};
                special.SatFreq_RBC = twix.hdr.MeasYaps.sWipMemBlock.adFree{5}; % ppm
                try
                    special.SatSpacing_RBC = twix.hdr.MeasYaps.sWipMemBlock.adFree{6}; % ms
                catch
                    special.SatSpacing_RBC = 0;
                end
            end
            
            % number frequency, and spacing of RBC saturations
            if ismember(twix.hdr.MeasYaps.sWipMemBlock.alFree{2}, [2 3])
                special.SatNum_TP = twix.hdr.MeasYaps.sWipMemBlock.alFree{11};
                special.SatFreq_TP = twix.hdr.MeasYaps.sWipMemBlock.adFree{12}; %  ppm
                try
                    special.SatSpacing_TP = twix.hdr.MeasYaps.sWipMemBlock.adFree{13}; % ms
                catch
                    special.SatSpacing_TP = 0;
                end
            end
            
            % saturation method
            if twix.hdr.MeasYaps.sWipMemBlock.alFree{2} == 1
                special.SatMethod = 'RBC';
            elseif twix.hdr.MeasYaps.sWipMemBlock.alFree{2} == 2
                special.SatMethod = 'TP';
            elseif twix.hdr.MeasYaps.sWipMemBlock.alFree{2} == 3
                special.SatMethod = 'RBC+TP';
            end
        end
    elseif contains(protocol,'fa_spec_')
        % tx factor
        special.txfactor = twix.hdr.MeasYaps.sWipMemBlock.adFree{1};

        % number of reps
        special.numReps = twix.hdr.MeasYaps.sWipMemBlock.alFree{2};
        
        % delay between images (ms)
        special.repDelay = twix.hdr.MeasYaps.sWipMemBlock.adFree{3};

        % number of resonances
        special.numRes = 1 + twix.hdr.MeasYaps.sWipMemBlock.alFree{4};      
        if isempty(special.numRes), special.numRes = 1; end

        % pulse shape for additional resonances
        special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
        if isempty(special.RFshape)
            special.RFshape = 'aSinc';
        elseif special.RFshape == 1
            special.RFshape = 'Gauss';
        elseif special.RFshape == 2
            special.RFshape = 'Rect';
        elseif special.RFshape == 3
            special.RFshape = 'Kai';
        end

        % acquisition order
        special.acqOrder = twix.hdr.MeasYaps.sWipMemBlock.alFree{10};
        if isempty(special.acqOrder)
            special.acqOrder = 'Sequential';
        elseif special.acqOrder == 1
            special.acqOrder = 'Interleaved';
        elseif special.acqOrder == 2
            special.acqOrder = 'Inter-Interleaved';
        end

        % frequency, flip angle, and pulse duration for other resonances
        special.RFdur = twix.hdr.MeasYaps.sWipMemBlock.adFree{9}; % ms
        if special.numRes >= 2
            special.freq2 = twix.hdr.MeasYaps.sWipMemBlock.adFree{5}; % ppm
            special.FA2 = twix.hdr.MeasYaps.sWipMemBlock.alFree{6}; % deg
        end

    elseif contains(protocol,'fa_radial_')
        % RF duration
        special.RFdur = twix.hdr.MeasYaps.sWipMemBlock.adFree{1}; % ms

        % RF pulse shape 
        special.RFshape = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
        if isempty(special.RFshape)
            special.RFshape = 'aSinc';
        elseif special.RFshape == 1
            special.RFshape = 'Gauss';
        elseif special.RFshape == 2
            special.RFshape = 'Rect';
        elseif special.RFshape == 3
            special.RFshape = 'Kai';
        end


    else
        for i = 1:14
            if isfield(twix.hdr.MeasYaps.sWipMemBlock,'alFree') &&...
                i <= length(twix.hdr.MeasYaps.sWipMemBlock.alFree) && ~isempty(twix.hdr.MeasYaps.sWipMemBlock.alFree{i})
            
                special.(['param',num2str(i)]) = twix.hdr.MeasYaps.sWipMemBlock.alFree{i};
            else
                if isfield(twix.hdr.MeasYaps.sWipMemBlock,'adFree') &&...
                        i <= length(twix.hdr.MeasYaps.sWipMemBlock.adFree) && ~isempty(twix.hdr.MeasYaps.sWipMemBlock.adFree{i})
                    
                    special.(['param',num2str(i)]) = twix.hdr.MeasYaps.sWipMemBlock.adFree{i};
                else
                    special.(['param',num2str(i)]) = 0;
                end
            end
        end
        
    end
    
    % zero any empty fields
    fnames = fieldnames(special);
    for i = 1:length(fnames)
        if isempty(special.(fnames{i}))
            special.(fnames{i}) = 0;
        end
    end
end
end
