function [hd,hd2] = plotmeansem(xax, mn, sm1, sm2,lineSpec1, lineSpec2, multiplier, rectify)
% use as plotmeansem(xax, mn, sm1, sm2,lineSpec1, lineSpec2, multiplier)

if nargin<7
  multiplier = 1; rectify = 0;
end
u = mn(:)+sm1(:);
l = mn(:)-sm2(:);
if rectify==1
l(l<0) = 0;
u(u<0) = 0;
mn(mn<0) = 0;
end
[X,Y] = polygon(xax,u,l,multiplier);

hd = plot(xax,mn,lineSpec1);
hold on
hd2 = plot(X,Y,lineSpec2);

function [X,Y] = polygon(x,sm1,sm2,multiplier)

x = x(:)';
if nargin < 4
    multiplier = 1;
end
X = [x x(end:-1:1) x(1)];
up = sm1(:)';
down = sm2(:)';
Y = [down up(end:-1:1) down(1)];
