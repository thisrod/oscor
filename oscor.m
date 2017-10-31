function oscor()
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
