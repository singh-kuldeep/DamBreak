close all
clc
clear
% One gas plots 
i = 0;
while(1)
	i = i+1;
	Array=csvread('premetive.csv');
	x = Array(:, 1);
	rho = Array(:, 2);
	u = Array(:, 3);
	p = Array(:, 4);
	T = p./(287.15*rho);

	% Plots
	handle = figure('units','normalized','outerposition',[0 0 1 1])
	subplot(2,2,1)
	plot(x, rho, 'b--o')
	title('Density v/s x')

	subplot(2,2,2)
	plot(x, u, 'r--o')
	title('Velocity v/s x')

	subplot(2,2,3)
	plot(x, p, 'g--o')
	title('Pressure v/s x')

	subplot(2,2,4)
	plot(x, T, 'b--o')
	title('Temperature v/s x')

	pause(3.5)

	% saving the plots
	saveas(handle,['filename' num2str(i) '.jpg']);
	close all
end

