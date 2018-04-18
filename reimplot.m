function reimplot(reh, imh, in)
%REIMPLOT greyscale plot of modulus from real and imaginary parts

R = groot;
B = R.Children([R.Children.Number]==reh).Children;
Sre = B.Children;
Sim = R.Children([R.Children.Number]==imh).Children.Children;

% keyboard
data = abs(Sre.ZData+1i*Sim.ZData);
data(:, Sre.XData==0) = 0;
figure, imagesc(Sre.XData, in.kc{2}, data), hold on, colormap gray

kk = num2cell(ylim);  kk = linspace(kk{:});
% two squeezing cycles per sound wave cycle
ww = abs(kk).*sqrt(kk.^2+2/in.c.healing^2);
plot(ww,kk,'-w', 'LineWidth', 2)
plot(xlim, sqrt(2)/in.c.healing*[1 -1; 1 -1], ':w', 'LineWidth', 2)
xlim([0 max(xlim)])

% text(0, 0, ...
% 	sprintf('range %.2e to %.2e', min(data(:)), max(data(:))), ...
% 	'Color', 'w')
% xlabel(B.XLabel.String), ylabel k

end