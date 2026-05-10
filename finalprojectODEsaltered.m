function dxdt = finalprojectODEsaltered(t,x,k,Aeff,pcd4,pcd70);

dxdt = zeros(6,1);

% free Ab  
dxdt(1) = -k(1)*x(1)*x(2) +k(2)*x(4) -k(3)*x(1)*x(3) +k(4)*x(5); 

% free CD4 receptor 
dxdt(2) = -k(1)*x(1)*x(2) +k(2)*x(4) -Aeff*pcd4*k(1)*x(5) +k(2)*x(6);  

% free CD70 receptor
dxdt(3) = -k(3)*x(1)*x(3) +k(4)*x(5) -Aeff*pcd70*k(3)*x(4) +k(4)*x(6);

% Ab:CD4a
dxdt(4) = k(1)*x(1)*x(2) -k(2)*x(4) -Aeff*pcd70*k(3)*x(4) +k(4)*x(6);

% Ab:CD70
dxdt(5) = k(3)*x(1)*x(3) -k(4)*x(5) -Aeff*pcd4*k(1)*x(5) +k(2)*x(6);

% CD4:Ab:CD70/second binding event
dxdt(6) = Aeff*pcd4*k(1)*x(5) -k(2)*x(6) +Aeff*pcd70*k(3)*x(4) -k(4)*x(6);  