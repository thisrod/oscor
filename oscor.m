function oscor()
% 1D coherence oscillations from a coherent and bogoliubov initial state

parameters

if exist('coherent')
	coherent = configure(coherent, R, L);
	coherent.randoms = 2;
	coherent.initial = @(w,r) 1/sqrt(2*r.c.gamma) + [1 1i]*w/2;
	coherent.file = 'coh.mat';
	xsim(coherent);
end

end

function in = configure(in, R, L)
	in.dimension = 2;
	in.points(2) = R;
	in.ranges = [0 L*(R-1)/R];  % allow for extra step of circular grid
	in.linear = @(r) 1i*r.Dx.^2;
	in.da = @(a,w,r) -1i*(sqrt(2*in.c.gamma)*abs(a).^2 - 1);
	in.ensembles = load('ensembles');
end