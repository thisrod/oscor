% plot graphs to compare +P and Wigner results

namefmt = 'posp/%s_M%d_gamma%.1g.mat';
offset = 0.1;

g = 0.03;
xgraph ~/Desktop/gstar/H/coht_0152_1960.mat
figure(9), A = gca;  C = A.Children;
load ~/Desktop/gstar/H/coht_0152_1960.mat
N = log2(input{1}.points(2));

figure
plot(C(1).XData, C(1).YData - (N-5)*offset, ':k', C(2).XData, C(2).YData - (N-5)*offset, ':k')
hold on

for n = 5:8
	load(sprintf(namefmt, 'wigner', 2^n, g), 'time', 'g2_mean', 'g2_err')
	g2_mean = g2_mean - (n-5)*offset;
	plot(time, g2_mean - g2_err, ':k', time, g2_mean + g2_err, ':k')
	load(sprintf(namefmt, 'posp', 2^n, g), 'time', 'g2_mean', 'g2_err')
	g2_mean = g2_mean - (n-5)*offset;
	plot(time, g2_mean - g2_err, '-k', time, g2_mean + g2_err, '-k')
end
xlim([min(time) max(time)])

for i = 1:26, close(i), end