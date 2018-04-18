% adjust time and frequency axes for factor of 2 between thesis and paper

figure(27); A = gca;
A.XTickLabel = arrayfun(@(n) num2str(n/2), A.XTick, 'UniformOutput', false);

figure(28); A = gca;
A.XTickLabel = arrayfun(@(n) num2str(n/2), A.XTick, 'UniformOutput', false);

figure(29); A = gca;
A.XTickLabel = arrayfun(@(n) num2str(2*n), A.XTick, 'UniformOutput', false);

figure(30); A = gca;
A.XTickLabel = arrayfun(@(n) num2str(2*n), A.XTick, 'UniformOutput', false);

figure(31); A = gca;
A.XTickLabel = arrayfun(@(n) num2str(2*n), A.XTick, 'UniformOutput', false);
