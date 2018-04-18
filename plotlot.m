function plotlot(s, sys)
%PLOTLOT generate the desired graphs from XSPDE oscor graphs
%  PLOTLOT(filename, label)

if nargin < 2
	sys = '';
end

if ~ismember(s(1), '/.')
	s = ['/Users/rpolkinghorne/Desktop/gstar/' s];
end
load(s),  xgraph(s), input = input{:};

reimplot(23,25,input)
xlabel \omega, ylabel k
h = colorbar;  set(get(h,'Label'),'String','S_T')
if sys, title(['KE distribution, ' sys]), end

reimplot(19,21,input)
xlabel \omega, ylabel k
h = colorbar;  set(get(h,'Label'),'String','S_n')
if sys, title(['phonon excitation, ' sys]), end

figure(17), wfplot, ylim([min(ylim),0])
xlabel s, ylabel k, zlabel 'n(k)', title(sys)

tote(15,16), xlabel s, ylabel 'T/L and V/L'
if sys, title(['CoE in ' sys]), end

figure(9), A = gca;  C = A.Children;
figure, plot(C(1).XData, C(1).YData, '-k', C(2).XData, C(2).YData, '-k')
xlabel s, ylabel 'g^{(2)}(0)',title(sys)

for i = 1:26, close(i), end

end