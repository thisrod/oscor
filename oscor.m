function oscor()

warning 'Aliasing test version'

% 1D coherence oscillations from a coherent and bogoliubov initial state

parameters, ensembles

% Comparison values are for the initial state

if exist('coherent')
	coherent = configure(coherent);
	coherent.randoms = 2;
	coherent.c.n = 1;  % density to normalise g2
	coherent.initial = @(w,r) 1 + [1 1i]*w/2;
	coherent.file = ['coh' coherent.file];
	xsim(coherent);
end

if exist('bogoliubov')
	bogoliubov = configure(bogoliubov);
	% this gets it wrong, thesis g2 graphs were normalised
	% by reading off the N/L graphs
	bogoliubov.c.n = bdens(0, bogoliubov);
	bogoliubov.initial = @binit;
	bogoliubov.file = ['bog' bogoliubov.file];
	xsim(bogoliubov);
end

end

function in = configure(in)
	in.dimension = 2;
	in.fields = [1 2];
	L = in.ranges(2);  R = in.points(2);
	in.ranges(2) = L*(R-1)/R;  % allow for extra step of circular grid
	in.linear = @(r) (1i/2)*r.Dx.^2;
	in.da = @(a,~,r) -(1i/2)*2*r.c.gamma*abs(a(1,:)).^2.*a(1,:);
	in.define = @ftit;
	in.file = sprintf('t_%04d_%03d.mat', ...
		round(-100*log10(in.c.gamma)), ...
		round(10*(in.ranges(2))));
	
 	obs = { ...
 		'Re \psi(x)' @(a,r) xave(real(a(1,:)), r); ...
 		'Im \psi(x)' @(a,r) xave(imag(a(1,:)), r); ...
 		'n(X)' @(a,r) abs(a(1,:)).^2 - 1/(2*r.dv); ...
 		'n(K)' @(a,r) abs(a(1,:)).^2 - 1/(2*r.dkv); ...
 		'n(X) not' @(a,r) a(2,:); ...
 		'g^{(2)}(x)' @(a,r) xave(abs(a(1,:)).^4 - 2*abs(a(1,:)).^2/r.dv + 1/(2*r.dv^2), r) / r.c.n^2; ...
 		'T(K)' @(a,r) r.kx.^2.*r.observe{4}(a,r); ...
 		'V(x)' @(a,r) r.c.gamma*r.c.n^2*r.observe{6}(a,r); ...
 		'N/L (x)' @(a,r) xint(r.observe{3}(a,r), r) / L; ...
 		'N/L (k)' @(a,r) xint(r.observe{4}(a,r), r.dk, r) / L; ...
 		'<T(k)>/L' @(a,r)  xint(r.observe{7}(a,r), r.dk, r) / L; ...
 		'<V(x)>/L' @(a,r)  xint(r.observe{8}(a,r), r.dx, r) / L; ...
 		'n(K)' @(a,r) r.observe{4}(a,r).*(r.kx ~= 0); ...
		'Re n(X) not' @(a,r) real(a(2,:)); ...
		'Im n(X) not' @(a,r) imag(a(2,:)); ...
		'Re T(X) not' @(a,r) real(a(3,:)); ...
		'Im T(X) not' @(a,r) imag(a(3,:)); ...
 	};
	in.olabels = obs(:,1)';  in.observe = obs(:,2)';

 	for j = 1:length(in.olabels)
 		if strfind(in.olabels{j}, '(x)')
 			in.transforms{j} = [false false];
 			in.pdimension{j} = 1;
 		elseif strfind(in.olabels{j}, '(X)')
 			in.transforms{j} = [false false];
 			in.pdimension{j} = 2;
 		elseif strfind(in.olabels{j}, '(k)')
 			in.transforms{j} = [false true];
 			in.pdimension{j} = 1;
 		elseif strfind(in.olabels{j}, '(K)')
 			in.transforms{j} = [false true];
 			in.pdimension{j} = 2;
 		else
 			warning('Label ''%s'' lacks x or k',  in.olabels{j})
 		end
 	end
	for i = 0:3
		in.transforms{length(in.olabels) - i} = [true false];
	end

end

function o = ftit(a,~,r)
	p = zeros(r.d.aplus);  p(1,:) = a(1,:);
	a = xgraphicsfft(p, [false true], r);
	a = a(1,:);
	o = [(abs(a).^2 - 1/(2*r.dkv)).*(r.kx ~= 0); ...
		r.kx.^2.*(abs(a).^2 - 1/(2*r.dkv))];
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
	for k = 2*pi*(1:ceil(n/4-1))/L	% sin and cos modes for each wave number
		kk = r.c.healing*k;
		uu = ((kk+1/kk)/sqrt(kk^2+2) + 1)/2;  vv = uu - 1;
		uu = sqrt(uu);  vv = sqrt(vv);
		z = [1 1i]*randn(2,2*ens)/2;  z = reshape(z,ens,2);
		a = a + (z(:,1)*uu-conj(z(:,1))*vv)*sin(k*r.xc{2});
		a = a + (z(:,2)*uu-conj(z(:,2))*vv)*cos(k*r.xc{2});
	end
	if false	% even grids have a zizag mode
		k = 2*pi*n/2/L;  kk = r.c.healing*k;
		uu = ((kk+1/kk)/sqrt(kk^2+2) + 1)/2;  vv = uu - 1;
		uu = sqrt(uu);  vv = sqrt(vv);
		z = [1 1i]*randn(2,ens)/2;  z = reshape(z,ens,1);
		a = a + (z(:,1)*uu-conj(z(:,1))*vv)*(-1).^(1:n)/sqrt(2);
	end
	a = sqrt(2/L)*a + 1;
%	a = sqrt(2/L)*a + (2*r.c.gamma)^(-1/4);	% healing length units
	a = reshape(a, r.d.a);
end % function init
