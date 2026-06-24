function [KSpaceCoor,KSpaceDiameter,DownRampLength] = calcspiralgrad_20190918(nsamples,nleaves,fov,orientation,larmor,traj)
% calculates the trajectory (1/mm) for the fa_spiral pulse sequence
% nsamples:                 number of sample points per interleave
% nleaves:                  number of interleaves
% orientation:              'axial','coronal', or 'sagittal'
% larmor:                   larmor frequency (MHz/T)
% traj:                     'arch' or 'var'

if nargin < 3, fov = []; end
if nargin < 4, orientation = []; end
if nargin < 5, larmor = []; end
if nargin < 6, traj = []; end
if isempty(fov), fov = 400; end
if isempty(orientation), orientation = 'axial'; end
if isempty(larmor), larmor = 42.58; end
if isempty(traj), traj = 'var'; end

%% Constant Parameters
risetime = 370;
dwelltime = 10;
fsgcm = 4; % fullscale G/cm
GAM = 100 * larmor; % Hz/G
A = 1;
ac = 0.94; 

%% Prepare gradient parameters
Gmax = 0;
S = dwelltime*A/risetime;
om = (2*pi/nleaves) * (fov/10)/(1/(GAM*fsgcm*dwelltime*1e-6));
npts = ceil(nsamples);

%% Generate gradient trajectory

switch lower(traj)
    case 'var'
        % Initial calculation
        gx = vsp(npts, nleaves, om, Gmax, S, A, ac);
        
        % Recalculate gradients
        DownRampLength = length(gx)-npts;
        [gx,gy,KSpaceDiameter] = vsp(npts-DownRampLength, nleaves, om, Gmax, S, A, ac);
        DownRampLength = DownRampLength + 1;
        
        % Pads gradient if too short
        gx(npts) = 0;
        gy(npts) = 0;
        
    case 'arch'
        dramplen = 50;
        [gx,gy] = archspiral(npts,nleaves,fov,larmor,dwelltime,dramplen);
        gx = gx(1:end-dramplen);
        gy = gy(1:end-dramplen);
        KSpaceDiameter = [];
        DownRampLength = 0;
end

% Scale gradients to mT/m
gx = gx * fsgcm * 10;
gy = gy * fsgcm * 10;


%% Create spiral interleaves
% Create rotation matrices
rotmat = zeros(3,3,nleaves);
switch lower(orientation)
    case 'axial'
        initrotmat =    [ 0   -1    0;
                         -1    0    0;
                          0    0   -1];                                   
    case 'coronal'
        initrotmat =    [ 1    0    0;
                          0    0   -1;
                          0    1    0];   
    case 'sagittal'
        initrotmat =    [ 0    0    1;
                         -1    0    0;
                          0   -1    0];
end

for i = 1:nleaves
    rotphi = -2 * pi * (i-1)/nleaves;
    rotcos = cos(rotphi);
    rotsin = sin(rotphi);
    
    R = [rotcos -rotsin 0;
         rotsin  rotcos 0;
         0       0      1];
        
    rotmat(:,:,i) = initrotmat * R;
end

ngx = zeros(nleaves, length(gx));
ngy = zeros(nleaves, length(gy));

% Rotate spirals
for ii = 1:nleaves
    g = rotmat(:,:,ii) * [gy; gx; ones(1,length(gx))];
    switch lower(orientation)
        case 'axial'
            ngx(ii,:) = g(1,:);
            ngy(ii,:) = g(2,:);
        case 'coronal'
            ngx(ii,:) = g(1,:);
            ngy(ii,:) = g(3,:);  
        case 'sagittal'
            ngx(ii,:) = g(2,:);
            ngy(ii,:) = g(3,:);
    end
end

% Get trajectory coordinates from integral of gradients
[kx_mat,ky_mat] = integral_fa(ngx, ngy, larmor);

% Concatenate interleaves
kx = reshape(kx_mat',1,[]);
ky = reshape(ky_mat',1,[]);

KSpaceCoor = [kx',ky'];

if isempty(KSpaceDiameter)
    KSpaceDiameter = max(sqrt(kx_mat(1,:).^2 + ky_mat(1,:).^2)) * 2 * fov;
end

end

%% Functions

function [gx,gy] = archspiral(npts,nleaves,fov,larmor,dwelltime,dramplen)

    lambda = nleaves/(2*pi*(fov)); % 1/mm
    theta = sqrt((0:npts-1).*dwelltime); % us
    
    gx = 1/larmor .* lambda .* (cos(theta)-theta.*sin(theta)) .* 100;
    gy = 1/larmor .* lambda .* (sin(theta)+theta.*cos(theta)) .* 100;
    
    gx(1) = 0;
    gy(1) = 0;
    
    stepx = -gx(npts)/(dramplen);
    stepy = -gy(npts)/(dramplen);
    for i = 1:dramplen
        gx(npts+i) = gx(npts) + (i*stepx);
        gy(npts+i) = gy(npts) + (i*stepy);
    end
    
end

function [gx,gy,kspaceDiameter] = vsp(npts, nleaves, OM, Gmax, S, A, ac)

maxdecratio = 32;
theta0 = 0;

if Gmax == 0
    vfactor = 1;
else
    vfactor = 1 - S/Gmax;
end

loop = 1;
decratio = 1;
while loop
    loop = 0;
    dnpts = npts*decratio;
    om = OM/decratio;
    s = S/decratio;
    g0 = 0;
    gx(1) = g0;
    gy(1) = 0;
    absg = hyp(g0,0);
    oldkx = 0;
    oldky = 0;
    kx = gx(1);
    ky = gy(1);
    thetan_1 = theta0;
    taun = 0;
    n = 0;
    
    while n < (dnpts-1)
        taun_1 = taun;
        taun = hyp(kx,ky)/A;
        tauhat = taun;
        theta = theta0 + atan2(om*tauhat,1) + (om*tauhat);
        
        if absg < ac
            deltheta = theta - thetan_1;
            B = 1/(1+tan(deltheta)*tan(deltheta));
            gtilde = absg*vfactor;
            t1 = s.^2;
            t2 = gtilde.*gtilde.*(1-B);
            if t2 > t1
                decratio = decratio*2;
                if (decratio > maxdecratio)
                    printf('Iteration failed');
                end
                loop = 1;
                break;
            end
            t3 = sqrt(t1-t2);
            absg = sqrt(B) * gtilde + t3;
        elseif absg > ac
            absg = ac;
        end
        tgx = absg*cos(theta);
        tgy = absg*sin(theta);
        kx = kx + tgx;
        ky = ky + tgy;
        thetan_1 = theta;
        
        if ~mod(n,round(decratio))
            m = n/round(decratio);
            if ((m+1)>(npts-1))
                break;
            end
            gx(m+2) = (kx-oldkx)/decratio;
            gy(m+2) = (ky-oldky)/decratio;
            oldkx = kx;
            oldky = ky;
        end
        n = n+1;
    end
end

kspaceDiameter = 2*nleaves*om*taun/(2*pi);
% disp(['K-space diameter = ',num2str(kspaceDiameter)]);

% downramp
t1 = atan2(gy(npts),gx(npts)) + pi;
delx = S*cos(t1) * .8;
dely = S*sin(t1) * .8;

if delx ~= 0
    m = round(abs(gx(npts)/delx));
else
    m = 0;
end

if dely ~= 0
    n = round(abs(gy(npts)/dely));
else
    n = 0;
end

if m > n
    dnpts = npts+n;
else
    dnpts = npts+m;
end

tgx = gx(npts);
tgy = gy(npts);

for n = npts+1:dnpts+1
    tgx = tgx + delx;
    tgy = tgy + dely;
    gx(n) = tgx;
    gy(n) = tgy;
end

n = n - 1;

while gx(n) ~=0 || gy(n) ~= 0 
    if abs(gx(n)) < abs(delx) || delx == 0
        gx(n+1) = 0;
    else
        gx(n+1) = gx(n) + delx;
    end
    if abs(gy(n)) < abs(dely) || dely == 0
        gy(n+1) = 0;
    else
        gy(n+1) = gy(n) + dely;  
    end
    n = n + 1;
end

dnpts = n + 1;

end

function z = hyp(x,y)
    z = sqrt(x.^2 + y.^2);
end

function [kx,ky] = integral_fa(gx, gy, larmor, dwelltime)
% gx,gy         nleaves x nsamples
% larmor        MHz/T
% dwelltime     us
% integral assumes equal spacing
if nargin < 3, larmor = []; end
if nargin < 4, dwelltime = []; end
if isempty(larmor), larmor = 42.58; end
if isempty(dwelltime), dwelltime = 10; end

kx = zeros(size(gx));
ky = zeros(size(gy));
for n = 1:size(gx,1)
    % mT/m  * 10G/1mT * 1m/1000mm * Hz/G * s
    kx(n,:) = cumsum(gx(n,:)).*(10/1000).*(larmor*100).*(dwelltime*1e-6);
    ky(n,:) = cumsum(gy(n,:)).*(10/1000).*(larmor*100).*(dwelltime*1e-6);
end

end