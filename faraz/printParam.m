function printParam(data,stitle)
% data          struct with necessary parameters
% stitle        title for acquisition (optional)

if nargin < 2, stitle = []; end

% print title (if provided)
if ~isempty(stitle)
    fprintf('\n%s \n',stitle);
end

% check if spiral or cartesian
if isfield(data,'trajectory')
    if strcmpi(data.trajectory,'Spiral')
        fprintf('\nSamples/Interleaves = %g/%g\n',data.nCol,data.nLin);
    else
        fprintf('\n');
    end
else
    fprintf('\n');
end

% print imaging parameters
fprintf('TR/TE = %g/%g ms\nFA = %g%c\n',data.TR,data.TE,data.FA,176);

% get dimensions
if isfield(data,'images')
    dim = size(data.images);         
else
    dim = [];
end

% print fov/res
if length(dim) > 2 && strcmp(data.mode,'2D')
    fov = [data.fovRO,data.fovPE,data.thickness*dim(3)];
    tmpstring2 = 'FOV = %g x %g x %g mm3\n';
    res = [data.fovRO/dim(1),data.fovPE/dim(2),data.thickness];
    tmpstring = 'Res = %.3g x %.3g x %.3g mm3\n';
    fprintf(tmpstring2,fov);
    fprintf(tmpstring,res);
    fprintf('%s \n',data.orientation);
elseif length(dim) > 2 && strcmp(data.mode,'3D')
    fov = [data.fovRO,data.fovPE,data.thickness*dim(3)];
    tmpstring2 = 'FOV = %g x %g x %g mm3\n';
    res = [data.fovRO/dim(1),data.fovPE/dim(2),data.thickness];
    tmpstring = 'Res = %.3g x %.3g x %.3g mm3\n';
    fprintf(tmpstring2,fov);
    fprintf(tmpstring,res);
    fprintf('3D %s \n',data.orientation);
elseif length(dim) == 2
    fov = [data.fovRO,data.fovPE];
    tmpstring2 = 'FOV = %g x %g mm2\n';
    res = [data.fovRO/dim(1),data.fovPE/dim(2)];
    tmpstring = 'Res = %.3g x %.3g mm2\n';
    fprintf(tmpstring2,fov);
    fprintf(tmpstring,res);
    if data.thickness == 500
        fprintf('%s Projection\n',data.orientation)
    else
        fprintf('Slice Thickness = %g mm\n',data.thickness);
        fprintf('%s \n',data.orientation);
    end 
elseif isempty(dim)
    if isfield(data,'MatSize') && ~isempty(data.MatSize)
        fov = [data.fovRO,data.fovPE,data.thickness*data.MatSize];
        res = fov ./ data.MatSize;
        tmpstring2 = 'FOV = %g x %g x %g mm3\n';
        tmpstring = 'Res = %.3g x %.3g x %.3g mm3\n';
        fprintf(tmpstring2,fov);
        fprintf(tmpstring,res);
    else
        fov = [data.fovRO,data.fovPE,data.thickness*max(data.nSli,data.nPar)];
        tmpstring2 = 'FOV = %g x %g x %g mm3\n';
        fprintf(tmpstring2,fov);
    end
    fprintf('%s %s \n',data.mode,data.orientation);
end

% extra space
fprintf('\n');

end