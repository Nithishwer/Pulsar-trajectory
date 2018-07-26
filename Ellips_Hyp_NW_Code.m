figure;
hold on

hypw0_e15      = hyperbol(0,1.5,"jl","--b");
hypw0_e2       = hyperbol(0,2,"jl","b");

legend(" e = 1.5"," e = 2")
xlabel(" Time - Time during Perihelion (s) ")
ylabel(" j_{pl} ( ms^{-3} ) ")

%function to return the properties of a hyperbolic trajectory for a linear f

function plotob = hyperbol(w,e,type,style)

vinfinity = 7.6217e+05;         %velocity at infinity
Mo        = 1.989*10^30;        %Solar Mass  
Mp        = 1.4*Mo;             %Pulsar Mass     
Mc        = 0.3*Mo;             %Companion Mass
i         = pi/3;               %Orbital inclination angle
Tp        = 0;                  %Epoch of the periastron Passage
G         = 6.674*10^(-11);     %Gravitational Constant
mu        = G*(Mc+Mp);          %mu
ap        = 388478173.5251782;  %semi-major axis length
n         = (mu/ap^3)^0.5;      %Mean Motion
nu        = (mu)^0.5/ap;        %nu
fmax      = acos(-1/e);         %maximum and minimum true anomaly   
apd       = ap*sin(i);          %Line of sight semi-major axis length
k=0;

syms f

r   =  ap*(e^2-1)/(1+e*cos(f));
rl  =  apd*(e^2-1)*(1+e*cos(f))^(-1)*sin(f+w);
v   =  n*ap*((e^2+1+2*e*cos(f))/(e^2-1))^0.5;
vl  =  n*apd/(e^2-1)^0.5*( cos(f+w) + e*cos(w) );
a   =  -n^2*ap/(e^2-1)^2*(1+e*cos(f))^2;
al  =  -n^2*apd/(e^2-1)^2*sin(f+w)*(1+e*cos(f))^2;
jl  =  -n^3*apd*(1+e*cos(f))^2/(e^2-1)^(7/2)*(cos(f+w)+e*cos(w)-3*e*sin(f)*sin(f+w));

range     =-fmax+fmax/200:fmax/200:fmax-fmax/200;
range_deg = range*180/pi

for fi=range
    
    k=k+1;
    F=2*atanh(((e-1)/(e+1))^0.5*tan(fi/2));
    M=e*sinh(F)-F;
    t= M/nu+Tp;
    Time(k)=t;
    
% Calculating los array for type property    
    
    if type == "r" 
            y(k)  = subs(r,f,fi);
    elseif type == "v"
            y(k) = subs(v,f,fi);
    elseif type == "a"
            y(k) = subs(a,f,fi);
    elseif type == "j"
            y(k) = subs(j,f,fi);
    elseif type == "rl"   
            y(k) = subs(rl,f,fi);
    elseif type == "vl"   
            y(k) = subs(vl,f,fi);
    elseif type == "al"
            y(k) = subs(al,f,fi);
    elseif type == "jl"
            y(k) = subs(jl,f,fi);
    end
    
end

plotob=plot(Time,y,style);

end



function plotob = ellip(w,e,style)

a=0;
Mo=1.989*10^30;       %Solar Mass
Mp=1.4*Mo;            %Pulsar Mass     
Mc=0.3*Mo;            %Companion Mass
Po=0.1*24*60*60;      %Pulsar Period
i=pi/3;               %Orbital inclination angle
Tp=0;                 %Epoch of the periastron Passage

if nargin<2
e=0.5;                %Orbital eccentricity
end

G=6.674*10^(-11);     %Gravitational Constant
wo=2*pi/Po;           %Orbital angular frequency =2*pi/Po
%w=0 ;                %Longitude of Periastron


syms f

apd=((Po/(2*pi))^2*G*(Mp+Mc))^(1/3)*(Mc*sin(i))/(Mp+Mc); 
ap=apd/sin(i)
r(f) =ap*(1-e^2)/(1+e*cos(f));
rl(f)=apd*(1-e^2)*(1+e*cos(f))^(-1)*sin(f+w);
vl(f)=(2*pi*apd*(cos(f+w)+e*cos(w)))/(Po*(1-e^2)^0.5);
al(f)=((-((2)*pi/Po)^2)*apd*sin(f+w)*(1+e*cos(f))^2)/(1-e^2)^2;
jl(f)=(((-2)*pi/Po)^3*apd/(1-e^2)^7/2)*(1+e*cos(f))^3*[cos(f+w)+e*cos(w)-3*e*sin(f+w)*sin(f)];


for fi=-pi:2*pi/100:pi
    
    a=a+1;
    
% Calculating los array for r, v, and j
       ra(a)=subs(r,f,fi);   
%     rla(a)=subs(rl,f,fi);
%     vla(a)=subs(vl,f,fi);
%     ala(a)=subs(al,f,fi);
%     jla(a)=subs(jl,f,fi);

end
plotob=plot(-pi:2*pi/100:pi,ra,style)

end

