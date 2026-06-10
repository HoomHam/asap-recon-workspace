function adKCoor = calcradialtraj_20220713(points,spokes,res,flatTime_spoke,flatTime_traverse)
% calculate trajectory for fa_radial sequences
% INPUTS: 
%    points             total number of readout points
%    spokes             total number of spokes
%    res                image resolution (1/mm)
%    flatTime_spoke     duration of plateau for out spoke
%    flatTime_traverse  duration of plateau for traverse spoke

%% Check inputs
if nargin < 4, flatTime_spoke = []; end
if nargin < 5, flatTime_traverse = []; end
if isempty(flatTime_spoke), flatTime_spoke = 800; end
if isempty(flatTime_traverse), flatTime_traverse = 1780; end

%% Control Parameters
dPhi1                = 0.4656;                  % Golden Mean 1
dPhi2                = 0.6823;                  % Golden Mean 2
PreSampleDwellTimes  = 20;                      % Number of dwell times during which spokes gradient is already ramping up before data sampling begins
lRuTimePreph_us      = 30;                      % Duration of ramp-up time for spoke prephaser
lRdTimePreph_us      = 30;                      % Duration of ramp-down time for spoke prephaser
lRuTimeSpoke_us      = 110;                     % Duration of ramp-up time for spoke
lFlatTimeSpoke_us    = flatTime_spoke;          % Duration of flat-top time for spoke
lRdTimeSpoke_us      = 110;                     % Duration of ramp-down time for spoke
lRuTimeTraverse_us   = 110;                     % Duration of ramp-up time for traverse
lFlatTimeTraverse_us = flatTime_traverse;       % Duration of flat-top time for traverse
lRdTimeTraverse_us   = 110;                     % Duration of ramp-down time for traverse
dPrephAmpl_mT_m      = -2.018408;               % Prephaser amp for 400/256 res
dSpokeAmpl_mT_m      = 8.325933;                % Spoke amp for 400/256 res
dTraverseAmpl_mT_m   = -7.953488;               % Traverse amp for 400/256 res
dGamma_MHz_T         = 42.5756;                 % Gyromagnetic ratio for 1H

%% Calculate center-out and right-left trajectories

% Calculate gradient amplitudes
dPrephAmpl_mT_m = dPrephAmpl_mT_m * ( (400/256) / res );
dSpokeAmpl_mT_m = dSpokeAmpl_mT_m * ( (400/256) / res );
dTraverseAmpl_mT_m = dTraverseAmpl_mT_m * ( (400/256) / res );

% number of sample points for each spoke 
samples_out = ceil(points/3);
samples_tra = floor(points/3)*2;


% Allocate memory
adKCoorSpokes = zeros(spokes,samples_out+PreSampleDwellTimes,3,'single');   % K-space coordinates for spokes
adKCoorTraverses = zeros(spokes,samples_tra,3,'single');                  % K-space coordinates for right-left traverses

% Calculate prephaser moment
dPrephMoment = 0.5 * double(lRuTimePreph_us + lRdTimePreph_us) * dPrephAmpl_mT_m;

% Calculate k-space coordinates of center-out spoke
dSpokeDuration_us = double(lRuTimeSpoke_us + lFlatTimeSpoke_us + lRdTimeSpoke_us);
dSpokeDurInc_us = dSpokeDuration_us / double(samples_out+PreSampleDwellTimes);
for i = 1:samples_out+PreSampleDwellTimes
    dTime_us = i * dSpokeDurInc_us;

    % Precalculate gradient ramps shorter than the dwell time
    if (i == 1 && dSpokeDurInc_us > double(lRuTimeSpoke_us))
        adKCoorSpokes(1,1,1) = 0.5 * dSpokeAmpl_mT_m * double(lRuTimeSpoke_us) + (dSpokeDurInc_us - double(lRuTimeSpoke_us)) * dSpokeAmpl_mT_m;
    else
        if (dTime_us <= lRuTimeSpoke_us)
            adKCoorSpokes(1,i,1) = (dTime_us * dSpokeAmpl_mT_m / double(lRuTimeSpoke_us)) * dTime_us * 0.5;
        elseif (dTime_us > lRuTimeSpoke_us && dTime_us <= lRuTimeSpoke_us + lFlatTimeSpoke_us)
            adKCoorSpokes(1,i,1) = adKCoorSpokes(1,i-1,1) + dSpokeDurInc_us * dSpokeAmpl_mT_m;
            lKCoorEndFlatTop = adKCoorSpokes(1,i,1);     % Remember gradient moment at the end of flat top
        else
            dTimeOnRamp_us = double(samples_out+PreSampleDwellTimes-i) * dSpokeDurInc_us;    % Time since end of flat top

            % Integrate down ramp by subtracting ramping-up momentum from total ramp momentum
            adKCoorSpokes(1,i,1) = lKCoorEndFlatTop + 0.5 * double(lRdTimeSpoke_us) * dSpokeAmpl_mT_m - (dTimeOnRamp_us * dSpokeAmpl_mT_m / double(lRuTimeSpoke_us)) * dTimeOnRamp_us * 0.5;
        end
    end
end
% Scale coordinates to mm^-1
adKCoorSpokes(1,:,1) = dGamma_MHz_T * (dPrephMoment + adKCoorSpokes(1,:,1)) / 1e6;
dOffset = adKCoorSpokes(1,end,1);

dTraverseDuration_us = double(lRuTimeTraverse_us + lFlatTimeTraverse_us + lRdTimeTraverse_us);
dTraverseDurInc_us = dTraverseDuration_us / double(samples_tra);
for i = 1:samples_tra
    dTime_us = i * dTraverseDurInc_us;

    % Precalculate gradient ramps shorter than the dwell time
    if (i == 1 && dTraverseDurInc_us > double(lRuTimeTraverse_us))
        adKCoorTraverses(1,1,1) = 0.5 * dTraverseAmpl_mT_m * double(lRuTimeTraverse_us) + (dTraverseDurInc_us - double(lRuTimeTraverse_us)) * dTraverseAmpl_mT_m;
    else
        if (dTime_us <= lRuTimeTraverse_us)
            adKCoorTraverses(1,i,1) = (dTime_us * dTraverseAmpl_mT_m / double(lRuTimeTraverse_us)) * dTime_us * 0.5;
        elseif (dTime_us > lRuTimeTraverse_us && dTime_us <= lRuTimeTraverse_us + lFlatTimeTraverse_us)
            adKCoorTraverses(1,i,1) = adKCoorTraverses(1,i-1,1) + dTraverseDurInc_us * dTraverseAmpl_mT_m;
            lKCoorEndFlatTop = adKCoorTraverses(1,i,1);        % Remember gradient moment at the end of flat top
        else
            dTimeOnRamp_us = double(samples_tra-i) * dTraverseDurInc_us;    % Time since end of flat top

            % Integrate down ramp by subtracting ramping-up momentum from total ramp momentum
            adKCoorTraverses(1,i,1) = lKCoorEndFlatTop + 0.5 * double(lRdTimeTraverse_us) * dTraverseAmpl_mT_m - (dTimeOnRamp_us * dTraverseAmpl_mT_m / double(lRuTimeTraverse_us)) * dTimeOnRamp_us * 0.5;
        end
    end
end
% Scale coordinates to mm^-1 and add offset momentum from outwards spoke
adKCoorTraverses(1,:,1) = dGamma_MHz_T * adKCoorTraverses(1,:,1) / 1e6 + dOffset;

%% Rotate spokes
for lSpoke = 2:spokes
    % Calculate new azimuth angle
    dAzimuth = mod(double(lSpoke-1) * 2.0 * pi * dPhi2, 2.0 * pi);
    
    % Calculate new z-coordinate
    zCoor = mod(double(lSpoke-1) * dPhi1 + 1.0, 2.0) - 1.0;
    
    % Calculate new polar angle
    dPolar = acos(zCoor);

    % Get rotations about Alpha and Beta
    Ry = [cos(dPolar) 0 -sin(dPolar); 
              0       1      0; 
          sin(dPolar) 0  cos(dPolar)];

    Rz = [cos(dAzimuth) -sin(dAzimuth) 0; 
          sin(dAzimuth)  cos(dAzimuth) 0; 
              0              0         1];

    rotmat = Ry * Rz;

    % Rotate
    adKCoorSpokes(lSpoke,:,:) = (rotmat * [adKCoorSpokes(lSpoke-1,:,1);
                                           adKCoorSpokes(lSpoke-1,:,2);
                                           adKCoorSpokes(lSpoke-1,:,3)])';

    adKCoorTraverses(lSpoke,:,:) = (rotmat * [adKCoorTraverses(lSpoke-1,:,1);
                                              adKCoorTraverses(lSpoke-1,:,2);
                                              adKCoorTraverses(lSpoke-1,:,3)])';

end
adKCoorSpokes = adKCoorSpokes(:,1+PreSampleDwellTimes:end,:);
adKCoor = reshape(cat(2,adKCoorSpokes,adKCoorTraverses),[],3);

end