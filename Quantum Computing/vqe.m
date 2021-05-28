X = 0:2*pi/1000:2*pi;
Y = zeros(1, 1001);
E = meas_H(pi);

for i = 1:1001
    Y(i)=meas_H(X(i));
end
    
plot(X,Y)
ylabel('Energy')
xlabel('Theta')
savefig('output/energyvtheta.fig')

[m,i] = max(Y);
Y(i) = 0;
[m,j] = max(Y);
disp('Minimum Energy:')
disp(m)
disp('Optimal Theta(s):')
disp(X(j))
disp(X(i))
