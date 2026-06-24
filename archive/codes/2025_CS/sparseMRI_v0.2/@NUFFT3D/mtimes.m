function ress = mtimes(a,bb)
% performs the normal nufft

% for n=1:size(bb,3) %% for 2d muti-slices??
% b = bb(:,:,n);  %% 

b = bb;
if a.adjoint  
	b = b(:).*a.w(:);
	res = nufft_adj(b, a.st)/sqrt(prod(a.imSize)); % adjoint : kspace to image
    % 	res = reshape(res, a.imSize(1), a.imSize(2));
	res = reshape(res, a.imSize(1), a.imSize(2), a.imSize(3)); %%
	res = res.*conj(a.phase);
	if a.mode==1
		res = real(res);
	end

else
	b = reshape(b,a.imSize(1),a.imSize(2), a.imSize(3));
	if a.mode==1
		b = real(b);
	end
	b = b.*a.phase;
	res = nufft(b, a.st)/sqrt(prod(a.imSize)); % image to k-space
	res = reshape(res,a.dataSize(1),a.dataSize(2));
end
% ress(:,:,n) = res;
ress = res; %%
    
% end

