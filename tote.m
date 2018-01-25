function tote(varargin)
%TOTE plot total of xspde graphs
%
%  TOTE(H, L, ...)  Plot the total of graphs H, L, ... 

S = 0;  V = 0;
figure, F = gcf;

% warning 'Fudge factor hasn''t been removed yet'

for i = varargin
	figure(i{1}), A = gca;
	t = A.Children(1).XData;
	y1 = A.Children(1).YData;
	y2 = A.Children(2).YData;
%	if i{1}==16
	if false
		y1 = y1*3;
		y2 = y2*3;
	end
	S = S + mean([y1; y2]);
	V = V + (y1-y2).^2;
	figure(F), plot(t,y1,':k',t,y2,':k'), hold on
end

plot(t,S+sqrt(V)/2,'-k',t,S-sqrt(V)/2,'-k')

end