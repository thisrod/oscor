function fire(letter, queue)
% FIRE  start oscor running on g2
%
%    FIRE(COMMAND, NAME, QUEUE) copies files from the current directory to
%    ~/NAME on the Green II cluster, then runs the Matlab script or
%    function COMMAND in that directory, using the specified QUEUE.
%    
%    This requires that ssh configured to login without any interaction.
%    Your cluster login is defined in a file parameters.m, which you will
%    need to edit. You will also need to edit a few other parameters,
%    including the names of the files to copy to the cluster.
%    
%    This is intended for use with XSPDE jobs, and it has some features to
%    assist with running them in parallel. The files parameters.m and
%    ensembles.m have a special meaning.
%    
%    The script parameters.m should define a set of XSPDE input structures,
%    or at least the points, steps and ensembles fields of those
%    structures. (You can set up the structures elsewhere, then run
%    parameters to assign these fields.) It should also list the names of
%    these structures in a cell array called fire_inputs.
%    
%    One of the parameters defined is hours_per_step. Run a small job to
%    start, then set this in parameters.m, so that fire will automatically
%    request an appropriate wall time.
%    
%    The fire script generates a batch script, and also a script
%    ensembles.m. The latter assigns ensembles fields to the input
%    structures in fire_inputs, so that each input structure has the same
%    number of trajectories as requested, but with parallelism appropriate
%    to QUEUE. Your parameters.m script should run this as its last step.
%    The idea is to write parameters.m to run on your desktop machine, and
%    let the generated ensembles script reconfigure it for the cluster. The
%    job control script will request parallel Matlab resources consistent
%    with the ensembles it generates.



switch queue
	case 'gstar', cpus = 10;
	case 'sstar', cpus = 14;
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
fprintf(runfile, '#PBS -l nodes=1:ppn=%d\n', cpus);
fprintf(runfile, '#PBS -l walltime=%d:00:00\n#PBS -N %s\n', ceil(2*work/1e9), letter);
fprintf(runfile, ['module load matlab/R2015b\ncd %s\n' ...
	'echo "Working directory: %s"\n' ...
	'matlab -r ''parpool(%d); xpsetup; oscor; quit''\n'], ...
	['~/' letter], letter, cpus);
fclose(runfile);

system(sprintf('ssh rpolking@g2.hpc.swin.edu.au /opt/torque/bin/qsub %s/runjob', letter));

end