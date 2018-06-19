%Two-body problem

m1 = 1;%1
m2 = 1;%1
M = 1000;%1000
r01 = 1;
r02 = 1;
v01 = 1;
v02 = 1;
Fx1 = [];
Fx2 = [];
Fy1 = [];
Fy2 = [];
G = 0.001;%m^3 kg^-1 s^-2
ux1 = [];
ux2 = [];
uy1 = [];
uy2 = [];
u1 = [];
u2 = [];
x1 = [];
x2 = [];
y1 = [];
y2 = [];
R1 = [];
R2 = [];
R12 = [];
E = [];
L = [];
i = 1;
dtau = 0.1;%0.01

x1(end+1) = 5.5/r01;% x1 = 5.5
x2(end+1) = 5/r02;% x2 = -3.5, set this to 5 ...
y1(end+1) = 0/r01;% y1 = 0
y2(end+1) = 0/r02;% y2 = 0
ux1(end+1) = 0/v01;% ux1 = 0
ux2(end+1) = 0/v02;% ux2 = 0
uy1(end+1) = 0.65/v01;% uy1 = 0.55, this to 0.65 ...
uy2(end+1) = 0.63/v02;% uy2 = -0.65, and this to 0.63 to see how the planets affect each other!
u1(end+1) = sqrt(ux1(i)^2 + uy1(i)^2);
u2(end+1) = sqrt(ux2(i)^2 + uy2(i)^2);
R1(end+1) = sqrt(x1(1)^2 + y1(1)^2);
R2(end+1) = sqrt(x2(1)^2 + y2(1)^2);
R12(end+1) = sqrt((x1(1)-x2(1))^2 + (y1(1) - y2(1))^2);
Fx1(end+1) = -G*M*(x1(1))/(r01*v01*v01*(R1(i))^3) ...
    - G*m2*(x1(1) - x2(1))/(r01*v01*v01*(R12(i))^3);
Fx2(end+1) = -G*M*(x2(1))/(r02*v02*v02*(R2(i))^3) ...
    + G*m1*(x1(1) - x2(1))/(r02*v02*v02*(R12(i))^3);
Fy1(end+1) = -G*M*(y1(1))/(r01*v01*v01*R1(i)^3) ...
    - G*m2*(y1(1) - y2(1))/(r01*v01*v01*R12(i)^3);
Fy2(end+1) = -G*M*(y2(1))/(r02*v02*v02*R2(i)^3) ...
    + G*m1*(y1(1) - y2(1))/(r02*v02*v02*R12(i)^3);


for t = 0:dtau:400
    ux1mid = ux1(i) + 0.5*r01*dtau*Fx1(i)/v01;
    ux2mid = ux2(i) + 0.5*r02*dtau*Fx2(i)/v02;
    uy1mid = uy1(i) + 0.5*r01*dtau*Fy1(i)/v01;
    uy2mid = uy2(i) + 0.5*r02*dtau*Fy2(i)/v02;
    
    x1(end+1) = x1(i) + ux1mid*dtau*r01/v01;
    x2(end+1) = x2(i) + ux2mid*dtau*r02/v02;
    y1(end+1) = y1(i) + uy1mid*dtau*r01/v01;
    y2(end+1) = y2(i) + uy2mid*dtau*r02/v02;
    
    %figure(1)
    %plot(x1, y1, 'b', x2, y2, 'g');
    %xlabel('x');
    %ylabel('y');
    %axis([-20 10 -10 10]);
    
    i = i + 1;
    
    R1(end+1) = sqrt(x1(i)^2 + y1(i)^2);
    R2(end+1) = sqrt(x2(i)^2 + y2(i)^2);
    R12(end+1) = sqrt((x1(i)-x2(i))^2 + (y1(i) - y2(i))^2);
    Fx1(end+1) = -G*M*(x1(i))/(r01*v01*v01*(R1(i))^3) ...
    - G*m2*(x1(i) - x2(i))/(r01*v01*v01*(R12(i))^3);
    Fx2(end+1) = -G*M*(x2(i))/(r02*v02*v02*(R2(i))^3) ...
    + G*m1*(x1(i) - x2(i))/(r02*v02*v02*(R12(i))^3);
    Fy1(end+1) = -G*M*(y1(i))/(r01*v01*v01*(R1(i)^3)) ...
    - G*m2*(y1(i) - y2(i))/(R12(i)^3);
    Fy2(end+1) = -G*M*(y2(i))/(r02*v02*v02*(R2(i)^3)) ...
    + G*m1*(y1(i) - y2(i))/(R12(i)^3);
    ux1(end+1) = ux1mid + dtau*Fx1(i);
    ux2(end+1) = ux2mid + dtau*Fx2(i);
    uy1(end+1) = uy1mid + dtau*Fy1(i);
    uy2(end+1) = uy2mid + dtau*Fy2(i);
    u1(end+1) = sqrt(ux1(i)^2 + uy1(i)^2);
    u2(end+1) = sqrt(ux2(i)^2 + uy2(i)^2);
    
    E(end+1) = 0.5*u1(i)*u1(i) + 0.5*u2(i)*u2(i) ...
        - (G*M*m1)/(R1(i)) - (G*M*m2)/(R2(i)) - (G*m1*m2)/(R12(i));
    L(end+1) = m1*x1(i)*uy1(i) - m1*y1(i)*ux1(i) ...
        + m2*x2(i)*uy2(i) - m2*y2(i)*ux2(i);
end


t = 0:dtau:400;

figure(1)
plot(x1, y1, 'b', x2, y2, 'g');
xlabel('x');
ylabel('y');
title('Planetary interaction');

figure(2)
plot(t, E, 'b', t, L, 'g');
xlabel('Time, t');
legend('Energy', 'Angular momentum');

