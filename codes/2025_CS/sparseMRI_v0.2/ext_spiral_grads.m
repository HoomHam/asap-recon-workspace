% rf pulse;
filename = 'C:\Users\P53-LOCAL\Desktop\AI\compressed_sensing\codes\sparseMRI_v0.2\VD_spiral_grad.txt';

% gradient waveform;
filename = 'C:\Users\P53-LOCAL\Desktop\SpatialSpect\3T_extrf\MRD1\c13grad_Pyr0_ss0.txt';

g_scale = max(max(g), abs(min(g))); % max(abs(g))
fid = fopen(filename, 'wt');
if fid ~= -1
    fprintf(fid, 'Gradient waveform for selective exciation\n');
    fprintf(fid, '%d\n', length(g));
    for idx=1:length(g)
        fprintf(fid, '%f\n', g(idx)/g_scale);
    end
    fprintf(fid,'\n');
    fprintf(fid, '%f\n', g_scale);
    fprintf(fid, '%f\n', z_thk);

    fclose(fid);
else
    warningMessage = sprintf('Cannot open file %s', filename);
    uiwait(warndlg(warningMessage));
end