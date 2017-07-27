function oscor()
% 1D coherence oscillations from a coherent and bogoliubov initial state

parameters, ensembles

% Comparison values are for the initial state

if exist('coherent')
	coherent = configure(coherent);
	coherent.randoms = 2;
	coherent.c.n = 1/sqrt(2*coherent.c.gamma);  % density to normalise g2
	coherent.initial = @(w,r) (2*r.c.gamma)^(-1/4) + [1 1i]*w/2;
	coherent.file = ['coh' coherent.file];
	xsim(coherent);
end

if exist('bogoliubov')
	bogoliubov = configure(bogoliubov);
	bogoliubov.c.n = bdens(0, bogoliubov);
	bogoliubov.initial = @binit;
	bogoliubov.file = ['bog' bogoliubov.file];
	xsim(bogoliubov);
end

end

% idea: plot T on the same scale where coherent g2 is 1

function in = configure(in)
	in.dimension = 2;
	L = in.ranges(2);  R = in.points(2);
	in.ranges(2) = L*(R-1)/R;  % allow for extra step of circular grid
	in.linear = @(r) 1i*r.Dx.^2;
	in.da = @(a,w,r) -1i*(sqrt(2*in.c.gamma)*abs(a).^2 - 1).*a;
	in.file = sprintf('_%04d_%03d.mat', ...
		round(-100*log10(in.c.gamma)), ...
		round(10*(in.ranges(2))));
	
	% Think about cutting off large k for the T calculation
	obs = { ...
		'Re \psi(x)' @(a,r) xave(real(a), r); ...
		'Im \psi(x)' @(a,r) xave(imag(a), r); ...
		'n(x)' @(a,r) abs(a).^2 - 1/(2*r.dv); ...
		'n(k)' @(a,r) abs(a).^2 - 1/(2*r.dkv); ...
		'g^{(2)}(x)' @(a,r) xave(abs(a).^4 - 2*abs(a).^2/r.dv + 1/(2*r.dv^2), r) / r.c.n^2; ...
		'T(k)' @(a,r) r.kx.^2.*(abs(r.kx)<3).*r.observe{4}(a,r); ...
		'<V(x)>/L' @(a,r) r.c.gamma*r.c.n^2*r.observe{5}(a,r); ...
		'N/L (x)' @(a,r) xint(r.observe{3}(a,r), r) / L; ...
		'N/L (k)' @(a,r) xint((abs(r.kx)<3).*r.observe{4}(a,r), r.dk, r) / L; ...
		'E - <T(k)>/L' @(a,r) r.c.gamma*r.c.n^2 - xint(r.observe{6}(a,r), r.dk, r) / L; ...
	};
	in.olabels = obs(:,1)';  in.observe = obs(:,2)';

	for j = 1:length(in.olabels)
		if strfind(in.olabels{j}, '(x)')
			in.transforms{j} = [false false];
		elseif strfind(in.olabels{j}, '(k)')
			in.transforms{j} = [false true];
		else
			error('Label ''%s'' lacks x or p',  in.transforms{j})
		end
	end
	in.pdimension = num2cell(ones(size(in.observe)));
	for j = [3 4 6], in.pdimension{j} = 2; end
	
end

function o = bdens(~,in,kcut)
	% density of bosons with wave number at most kcut in a Bogoliubov state
	% kcut defaults to the Nyquist wave number of the grid
	L = in.ranges(2);  R = in.points(2);
	L = L*R/(R-1);
	if nargin > 2
		n = kcut*L/(2*pi);
	else
		n = R/2 - 1;
	end
	k = 2*pi/L*(1:n);
	vv = (k+1./k)./(2*sqrt(k.^2+2)) - 1/2;
	o = 1/sqrt(2*in.c.gamma) + 2*sum(vv)/L;
end

function a = binit(~,r)
	L = r.ranges(2);  R = r.points(2);
	L = L*R/(R-1);  % undo circular grid adjustment
	n = r.nspace;  ens = r.ensembles(1);
	a = zeros(ens, n);
	% the condensate mode k=0 is added later
	for kk = 2*pi*(1:ceil(n/2-1))/L	% sin and cos modes for each wave number
		uu = ((kk+1/kk)/sqrt(kk^2+2) + 1)/2;  vv = uu - 1;
		uu = sqrt(uu);  vv = sqrt(vv);
		z = [1 1i]*randn(2,2*ens)/2;  z = reshape(z,ens,2);
		a = a + (z(:,1)*uu-conj(z(:,1))*vv)*sin(kk*r.xc{2});
		a = a + (z(:,2)*uu-conj(z(:,2))*vv)*cos(kk*r.xc{2});
	end
	if mod(n,2) == 0	% even grids have a zizag mode
		kk = 2*pi*n/2/L;
		uu = ((kk+1/kk)/sqrt(kk^2+2) + 1)/2;  vv = uu - 1;
		uu = sqrt(uu);  vv = sqrt(vv);
		z = [1 1i]*randn(2,ens)/2;  z = reshape(z,ens,1);
		a = a + (z(:,1)*uu-conj(z(:,1))*vv)*(-1).^(1:n);
	end
	a = sqrt(2/L)*a + (2*r.c.gamma)^(-1/4);
	a = reshape(a, r.d.a);
end % function init
