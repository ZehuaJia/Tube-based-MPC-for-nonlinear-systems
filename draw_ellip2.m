function  h = draw_ellip2(P,alpha,x,color)
%% Draw the ellipse centered at origin
% The ellipsoid£ºa * x^2 + c * y^2 + b * x * y = f 
% P = [a, b/2; b/2, c]
a = P(1, 1);
b = 2 * P(1, 2);
c = P(2, 2);
d = 0;
e = 0;
f = alpha;
delta = b^2-4*a*c;
if delta >= 0
    warning('This is not an ellipse')
    return;
end
x0 = (b*e-2*c*d)/delta;
y0 = (b*d-2*a*e)/delta;
r = a*x0^2 + b*x0*y0 +c*y0^2 + f;
if r <= 0
    warning('This is not an ellipse')
    return;
end

aa = sqrt(r/a); 
bb = sqrt(-4*a*r/delta);
t = linspace(0, 2*pi, 60);
xy = [1 -b/(2*a);0 1]*[aa*cos(t);bb*sin(t)];
plot(xy(1,:)+x(1),xy(2,:)+x(2), color, 'linewidth', 2);