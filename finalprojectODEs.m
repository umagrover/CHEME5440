function dxdt = finalprojectODEs(t,x,k,kappa);

dxdt = zeros(6,1);

% free Ab  
dxdt(1) = -k(1)*x(1)*x(2) +k(2)*x(4) -k(3)*x(1)*x(3) +k(4)*x(5); 

% free CD4 receptor 
dxdt(2) = -k(1)*x(1)*x(2) +k(2)*x(4) -kappa*k(1)*x(5)*x(2) +k(2)*x(6);  

% free CD70 receptor
dxdt(3) = -k(3)*x(1)*x(3) +k(4)*x(5) -kappa*k(3)*x(4)*x(3) +k(4)*x(6);

% Ab:CD4
dxdt(4) = k(1)*x(1)*x(2) -k(2)*x(4) -kappa*k(3)*x(4)*x(3) +k(4)*x(6);

% Ab:CD70
dxdt(5) = k(3)*x(1)*x(3) -k(4)*x(5) -kappa*k(1)*x(5)*x(2) +k(2)*x(6);

% CD4:Ab:CD70
dxdt(6) = kappa*k(1)*x(5)*x(2) -k(2)*x(6) +kappa*k(3)*x(4)*x(3) -k(4)*x(6);  