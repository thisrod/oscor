%TPLOT greyscale plot of kinetic energy distribution over momentum modes

B = gca;  A = B.Children;
figure, imagesc(A.XData, A.YData, abs(A.ZData)), colormap gray
text(min(A.XData(:)), 0, ...
	sprintf('range %.2e to %.2e', min(A.ZData(:)), max(A.ZData(:))), ...
	'Color', 'w')
xlabel(B.XLabel.String), ylabel(B.YLabel.String)
fprintf('Remember to replace YData with momentum grid\n')