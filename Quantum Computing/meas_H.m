function E = meas_H(theta)
X = [0,1;1,0];
Y = [0,-1i;1i,0];
Z = [1,0;0,-1];
I = [1,0;0,1];
M = cos(theta)*kron(I,I) - sin(theta)*1i*kron(X,Y);
phi = M*[0;1;0;0];

c0 = -1.47*kron(I,I);
c1 = 0.442*kron(Z,I);
c2 = -0.246*kron(I,Z);
c3 = 0.373*kron(Z,Z);
c4 = 0.191*kron(Y,Y);
c5 = 0.191*kron(X,X);
c6 = 0.01*kron(X,Y);

E = phi'*(c0+c1+c2+c3+c4+c5+c6)*phi;
end