%Code for a single particle in a central force field:

m = 1;
dtau = 0.1;%0.01
r0 = 1;
v0 = 1;
k = 1;
i = 1;
E = [];
Fr = [];
Fx = [];
Fy = [];
R = [];
x = [];
y = [];
u = [];
ux = [];
uy = [];
L = [];
midKinEnergy = 0;
midPotEnergy = 0;
midEnergy = zeros(1, 30001);
T = 13;
%theta = [];


x(end+1) = 4;%"Perfect" ellipse for x = 4, y = 6, ux = 0.3 and uy = -0.1
y(end+1) = 6;%"Perfect" circle for same x and y (or x = 4 = y), ux = 0.3 = uy
R(end+1) = sqrt(x(1)^2 + y(1)^2);
ux(end+1) = 0.3/v0;
uy(end+1) = 0/v0;
u(end+1) = sqrt(ux(1)^2 + uy(1)^2);%Dimensionless
%Fr(end+1) = -k/(v0^2*r0*R(1)^2);%Dimensionless
Fx(end+1) = -k*x(1)/(v0^2*r0^2*R(1)^3);%Dimensionless
Fy(end+1) = -k*y(1)/(v0^2*r0^2*R(1)^3);%Dimensionless
%theta(end+1) = 0;

for t = 0:dtau:300
   
   uxmid = ux(i) + Fx(i)*dtau*r0/2*m*v0;
   uymid = uy(i) + Fy(i)*dtau*r0/2*m*v0;
   x(end+1) = x(i) + r0*dtau*ux(i)/v0 + Fx(i)*(dtau^2)*r0^2/2*m*v0^2;
   y(end+1) = y(i) + r0*dtau*uy(i)/v0 + Fy(i)*(dtau^2)*r0^2/2*m*v0^2;
   
   figure(1)
   plot(x, y);
   xlabel('x');
   ylabel('y');
   axis([-3 9 -6 8]);
   
   i = i + 1;
   R(end+1) = sqrt(x(i)^2+y(i)^2);
   Fx(end+1) = -k*x(i)/(v0^2*r0*R(i)^3);
   Fy(end+1) = -k*y(i)/(v0^2*r0*R(i)^3);
   ux(end+1) = uxmid + Fx(i)*r0*dtau/2*m*v0;
   uy(end+1) = uymid + Fy(i)*r0*dtau/2*m*v0;
   u(end+1) = sqrt(ux(i)^2 + uy(i)^2);
   E(end+1) = 0.5*u(i)*u(i) - k/(r0*v0^2*m*R(i));
   midKinEnergy = midKinEnergy + 0.5*u(i)*u(i);
   midPotEnergy = midPotEnergy - k/(r0*v0^2*m*R(i));
   L(end+1) = x(i)*m*uy(i) - m*y(i)*ux(i);%Only contribution from this component of L
   
end

t = 0:dtau:300;

figure(1)
plot(x, y);
xlabel('x');
ylabel('y');

midKinEnergy = midKinEnergy/T;
midPotEnergy = midPotEnergy/T;
midEnergy(:) = midKinEnergy/midPotEnergy;

figure(2)
plot(t, E, t, L, 'g', t, midEnergy, 'r');
xlabel('Time, t');
ylabel('Energy, E(u, R)');
axis([0 300 -4 1])
legend('Energy', 'Angular momentum', 'Kinetic/potential ratio');
