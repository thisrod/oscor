function wfplot()
%WFPLOT convert an XPDE plot to a waterfall chart

A = gca;  C = A.Children;
figure, S = surf(C.XData, C.YData, C.ZData);
S.MeshStyle = 'row';  S.EdgeColor = 'k';  S.FaceColor = 'none';
xlabel(A.XLabel.String), ylabel(A.YLabel.String), zlabel(A.ZLabel.String)

end