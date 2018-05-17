% plot graphs to compare +P and Wigner results

namefmt = 'posp/%s_M%d_gamma%.1g.mat';
offset = 0.1;

gs = [0.01 0.03 0.1];
pths = {'DD/coht_0200_1884.mat' ...
	'H/coht_0152_1960.mat' 'CC/coht_0100_1936.mat'};
abc = 'abc';

for k = 1:3
	g = gs(k);
	pth = ['~/Desktop/gstar/' pths{k}];
	xgraph(pth)
	figure(9), A = gca;  C = A.Children;
	load(pth)
	N = log2(input{1}.points(2));
	
	figure
	for j = 1:2
		% factor of 2 for paperunits
		plot(C(j).XData/2, C(j).YData - (N-5)*offset, '--k', 'LineWidth', 2)
		hold on
	end
	
	for n = 5:8
		load(sprintf(namefmt, 'wigner', 2^n, g), 'time', 'g2_mean', 'g2_err')
		g2_mean = g2_mean - (n-5)*offset;
		plot(time, g2_mean - g2_err, ':k', time, g2_mean + g2_err, ':k')
		load(sprintf(namefmt, 'posp', 2^n, g), 'time', 'g2_mean', 'g2_err')
		g2_mean = g2_mean - (n-5)*offset;
		plot(time, g2_mean - g2_err, '-k', time, g2_mean + g2_err, '-k')
	end
	A = gca;
	xlim([min(time) max(time)]), ylim([0.55 1.05])
	A.YTick = [0.6 0.7 0.73 0.8 0.83 0.9 0.93 1];
	A.YTickLabel = {'0.9' '1' '0.93' '1' '0.93' '1' '0.93' '1'};
	xlabel s, ylabel g^{(2)}(0)
	
	for i = 1:26, close(i), end
	saveTightFigure(['resp190418' abc(k) '.pdf'])
end