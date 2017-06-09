function fire(letter, queue)
% FIRE  start oscor running on g2

switch queue
	case 'gstar', cpus = 10;
	otherwise error 'Unknown job queue'
end

parameters

if system(['mount | grep' ...
	' "rpolking@g2\.hpc\.swin\.edu\.au.*/Users/rpolkinghorne/Desktop/gstar" ' ...
	'> /dev/null']) ~= 0
	system 'sshfs rpolking@g2.hpc.swin.edu.au: ~/Desktop/gstar'
end

if ismember(exist(['~/Desktop/gstar/' letter]), [2 7])
	error 'Target directory already exits'
else
	mkdir(['~/Desktop/gstar/' letter])
end

copyfile('oscor.m', ['~/Desktop/gstar/' letter]);
copyfile('parameters.m', ['~/Desktop/gstar/' letter]);
xpfile = fopen(['~/Desktop/gstar/' letter '/xpsetup.m'], 'w');
fprintf(xpfile, 'addpath(genpath(''/home/rpolking/xspde_matlab''))\n');
fclose(xpfile);

ensfile = fopen(['~/Desktop/gstar/' letter '/ensembles.m'], 'w');
work = 0;
if exist('coherent')
	jobs = prod(coherent.ensembles(2:3));
	runs = ceil(jobs/cpus);
	fprintf(ensfile, 'coherent.ensembles = [%d %d %d];\n', ...
		coherent.ensembles(1), runs, cpus);
	work = work + runs*coherent.ensembles(1)*prod(coherent.points)*coherent.steps;
end
if exist('bogoliubov')
	jobs = prod(bogoliubov.ensembles(2:3));
	runs = ceil(jobs/cpus);
	fprintf(ensfile, 'bogoliubov.ensembles = [%d %d %d];\n', ...
		bogoliubov.ensembles(1), runs, cpus);
	work = work + runs*bogoliubov.ensembles(1)*prod(bogoliubov.points)*bogoliubov.steps;
end
fclose(ensfile);

runfile = fopen(['~/Desktop/gstar/' letter '/runjob'], 'w');
fprintf(runfile, '#!/bin/sh\n#PBS -q %s\n', queue);
fprintf(runfile, '#PBS -l nodes=1:ppn=%d\n#PBS -l walltime=%d:00:00\n', cpus, ceil(work/1e9));
fprintf(runfile, ['module load matlab/R2015b\ncd %s\n' ...
	'matlab -r ''parpool(%d); xpsetup; oscor; quit''\n'], ...
	['~/Desktop/gstar/' letter], cpus);
fclose(runfile);

system(sprintf('ssh rpolking@g2.hpc.swin.edu.au /opt/torque/bin/qsub %s/runjob', letter));

end