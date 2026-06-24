function res = adjD(y)

res = zeros(size(y,1),size(y,2),size(y,3),size(y,4));

res = adjDx(y(:,:,:,:,1)) + adjDy(y(:,:,:,:,2)) + adjDz(y(:,:,:,:,3)) + adjDt(y(:,:,:,:,4));

return;


function res = adjDx(x)
res = x([1,1:end-1],:,:,:) - x;
res(1,:,:,:) = -x(1,:,:,:);
res(end,:,:,:) = x(end-1,:,:,:);

function res = adjDy(x)
res = x(:,[1,1:end-1],:,:) - x;
res(:,1,:,:) = -x(:,1,:,:);
res(:,end,:,:) = x(:,end-1,:,:);

function res = adjDz(x)
res = x(:,:,[1,1:end-1],:) - x;
res(:,:,1,:) = -x(:,:,1,:);
res(:,:,end,:) = x(:,:,end-1,:);

function res = adjDt(x)
res = x(:,:,:,[1,1:end-1]) - x;
res(:,:,:,1) = -x(:,:,:,1);
res(:,:,:,end) = x(:,:,:,end-1);

% % 
% res = zeros(size(y,1),size(y,2),size(y,3),size(y,5));
% 
% for fr=1:size(y,5)
%     res(:,:,:,fr) = adjDx(y(:,:,:,1,fr)) + adjDy(y(:,:,:,2,fr)) + adjDz(y(:,:,:,3,fr));
% end
% 
% return;
% 
% 
% function res = adjDx(x)
% res = x([1,1:end-1],:,:) - x;
% res(1,:,:) = -x(1,:,:);
% res(end,:,:) = x(end-1,:,:);
% 
% function res = adjDy(x)
% res = x(:,[1,1:end-1],:) - x;
% res(:,1,:) = -x(:,1,:);
% res(:,end,:) = x(:,end-1,:);
% 
% function res = adjDz(x)
% res = x(:,:,[1,1:end-1]) - x;
% res(:,:,1) = -x(:,:,1);
% res(:,:,end) = x(:,:,end-1);
% 
% function res = adjDt(x)
% res = x(:,:,[1,1:end-1]) - x;
% res(:,:,1) = -x(:,:,1);
% res(:,:,end) = x(:,:,end-1);
