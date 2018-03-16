function c = myconv2(a,b)
% Two dimensional convolution.

vpad = ceil((size(b,1) - 1)/2);
hpad = ceil((size(b,2) - 1)/2);
u = repmat(a(1,:),[vpad 1]);
d = repmat(a(end,:),[vpad 1]);
l = repmat(a(:,1),[1 hpad]);
r = repmat(a(:,end),[1 hpad]);
ul = repmat(a(1,1),[vpad hpad]);
ur = repmat(a(1,end),[vpad hpad]);
dl = repmat(a(end,1),[vpad hpad]);
dr = repmat(a(end,end),[vpad hpad]);
ap = [ul u ur; l  a r; dl d dr];
cp = conv2(ap,b,'same');
c = cp(vpad+1:vpad+size(a,1),hpad+1:hpad+size(a,2));