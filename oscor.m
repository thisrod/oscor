function oscor()
% 1D coherence oscillations from a coherent and bogoliubov initial state

parameters

if exist('coherent')
	coherent = configure(coherent);
	coherent.randoms = 2;
	coherent.initial = @(w,r) (2*r.c.gamma)^(-1/4) + [1 1i]*w/2;
	coherent.compare{3} = @(t,in) 1/sqrt(2*in.c.gamma);
	coherent.compare{4} = @(t,in) 1;
	coherent.file = 'coh.mat';
	xspde(coherent);
end

end

function in = configure(in)
	in.dimension = 2;
	L = in.ranges(2);  R = in.points(2);
	in.ranges(2) = L*(R-1)/R;  % allow for extra step of circular grid
	in.linear = @(r) 1i*r.Dx.^2;
	in.da = @(a,w,r) -1i*(sqrt(2*in.c.gamma)*abs(a).^2 - 1).*a;
	in.ensembles = load('ensembles');
	
	in.olabels = { ...
		'Re \psi' ...
		'Im \psi' ...
		'n' ...
		'g^2(0)' ...
% 		'|\delta\psi|^2' ...
% 		'(Re \delta\psi)^2' ...
% 		'(Im \delta\psi)^2' ...
% 		'|\delta\psi|^4' ...
% 		'|\delta\psi|^2Re \delta\psi' ...
% 		'Re \delta\psi' ...
% 		'Im \delta\psi' ...
	};
	in.observe = { ...
		@(a,~) real(a) ...
		@(a,~) imag(a) ...
		@(a,r) abs(a).^2 - 1/(2*r.dv) ...
		@(a,r) 2*r.c.gamma*(abs(a).^4 - 2*abs(a).^2/r.dv + 1/(2*r.dv^2)) ...
% 		@(a,~) abs(a-1).^2 ...
% 		@(a,~) real(a-1).^2 ...
% 		@(a,~) imag(a-1).^2 ...
% 		@(a,~) abs(a-1).^4 ...
% 		@(a,~) abs(a-1).^2.*real(a-1) ...
% 		@(a,~) real(a-1) ...
% 		@(a,~) imag(a) ...
	};
end