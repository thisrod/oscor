function oscor()
% 1D coherence oscillations from a coherent and bogoliubov initial state

parameters

% Comparison values are for the initial state

if exist('coherent')
	coherent = configure(coherent);
	coherent.randoms = 2;
	coherent.initial = @(w,r) (2*r.c.gamma)^(-1/4) + [1 1i]*w/2;
	coherent.compare{3} = @(t,in) 1/sqrt(2*in.c.gamma);
	coherent.compare{4} = @(t,in) 1;
	coherent.file = 'coh.mat';
	xspde(coherent);
end

if exist('bogoliubov')
	bogoliubov = configure(bogoliubov);
%	bogoliubov.randoms = 2;
	bogoliubov.initial = @binit;
%	bogoliubov.compare{3} = @(t,in) 1/sqrt(2*in.c.gamma);
%	bogoliubov.compare{4} = @(t,in) 1;
	bogoliubov.file = 'bog.mat';
	xspde(bogoliubov);
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
		@(a,r) xave(abs(a).^2, r) - 1/(2*r.dv) ...
		@(a,r) 2*r.c.gamma*xave(abs(a).^4 - 2*abs(a).^2/r.dv + 1/(2*r.dv^2), r) ...
% 		@(a,~) abs(a-1).^2 ...
% 		@(a,~) real(a-1).^2 ...
% 		@(a,~) imag(a-1).^2 ...
% 		@(a,~) abs(a-1).^4 ...
% 		@(a,~) abs(a-1).^2.*real(a-1) ...
% 		@(a,~) real(a-1) ...
% 		@(a,~) imag(a) ...
	};
	in.pdimension{3} = 1;
	in.pdimension{4} = 1;
end

function a = binit(~,r)
	L = r.ranges(2);  R = r.points(2);
	L = L*R/(R-1);  % undo circular grid adjustment
	n = r.nspace;  ens = r.ensembles(1);
	a = zeros(ens, n);
	% skip k=0 because that's the condensate mode
	assert(mod(n,2)==0, ...
		'Implementation restriction: Bogoliubov grids must have an even number of points.')
	for kk = 2*pi*(1:n/2-1)/L	% sin and cos modes for each wave number
		uu = ((kk+1/kk)/sqrt(kk^2+2) + 1)/2;  vv = uu - 1;
		uu = sqrt(uu);  vv = sqrt(vv);
		z = [1 1i]*randn(2,2*ens)/2;  z = reshape(z,ens,2);
		a = a + (z(:,1)*uu-conj(z(:,1))*vv)*sin(kk*r.xc{2});
		a = a + (z(:,2)*uu-conj(z(:,2))*vv)*cos(kk*r.xc{2});
	end
	kk = 2*pi*n/2/L;		% zizag mode at the Nyquist wave number
	uu = ((kk+1/kk)/sqrt(kk^2+2) + 1)/2;  vv = uu - 1;
	uu = sqrt(uu);  vv = sqrt(vv);
	z = [1 1i]*randn(2,ens)/2;  z = reshape(z,ens,1);
	a = a + (z(:,1)*uu-conj(z(:,1))*vv)*(-1).^(1:n);
	a = sqrt(2/L)*a + (2*r.c.gamma)^(-1/4);	% order parameter
	a = reshape(a, r.d.a);
end % function init
