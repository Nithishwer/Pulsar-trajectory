clc;

E = 1.2;

plot_param(1.5);

function func = plot_param(E)
%Initialize all constants
M0 = 1.989e30;
Mp = 1.4*M0;
Mc = 10e-8*M0;
i = 60.0*pi/180.0;
G = 6.67e-11;
v_inf = 5000.0;
omega = 0;
t_init = 0;

%Calculate all basic parameters
mu = G * (Mp + Mc);
a = mu ./ v_inf.^2;
ap = a * sin(i);
nu = (mu ./ a.^3).^0.5;
eterm = E.^2 - 1;
fmax = acos(-1/E);
%Calculate time array from the true anomaly
k=1;
t = zeros;
for f=-fmax+fmax/100:fmax/100.0:fmax-fmax/100
    H = 2* atan(((E-1)/(E+1))^0.5 * tan(f/2.0));
    M = E*sinh(H) - H;
    t(k) = M/nu + t_init;
    k=k+1;
end

f=(-fmax+fmax/100:fmax/100.0:fmax-fmax/100);

%Calculate all motion parameters from the true anomaly
rl = ap * eterm * sin(f + omega) ./ (1+E*cos(f));
vl = nu * ap / eterm^0.5 * (cos(f + omega) + E * cos(omega));
al = - nu^2 * ap * sin(f + omega) .* (1 + E * cos(f)).^2 / eterm^2;
jl = - nu^3 * ap / eterm^(7.0/2.0) * (1 + E * cos(f)).^3 .* (cos(f + omega) + E * cos(omega) - 3.0 * E * sin(f + omega) .* sin(f));

%Plot line-of-sight radius, velocity, acceleration and jerk wrt time
figure(1);
plot(t, rl);
hold on
title('Graph of Line-of-sight radius with time for hyperbola');
xlabel('Time (s)');
ylabel('Radius');
hold off

figure(2);
plot(t, vl);
hold on
title('Graph of Line-of-sight velocity with time for hyperbola');
xlabel('Time (s)');
ylabel('Velocity');
hold off

figure(3);
plot(t, al);
hold on
title('Graph of Line-of-sight acceleration with time for hyperbola');
xlabel('Time (s)');
ylabel('Acceleration');
hold off

figure(4);
plot(t, jl);
hold on
title('Graph of Line-of-sight jerk with time for hyperbola');
xlabel('Time (s)');
ylabel('Jerk');
hold off

end