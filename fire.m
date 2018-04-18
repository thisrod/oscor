function fire(command, letter, queue)
% FIRE  start an XSPDE batch job on Green II
%
%    FIRE(COMMAND, NAME, QUEUE) copies files from the current directory
%    to ~/NAME on the Green II cluster, then runs the Matlab script or
%    function COMMAND in that directory, using the specified QUEUE.
%    
%    You must have ssh configured to login to gstar without prompting
%    for a password.
%    
%    There must be a script called parameters.m in the current directory.
%    This needs to assign the following Matlab variables:
%    
%    gstar_login = your username on Green II 
%    
%    fire_files = a cell array of filenames to copy to Green II
%    
%    fire_inputs = a cell array with the names of XSPDE input structures
%    defined in parameters.m
%    
%    hours_per_step = wall time to request per gridpoint step, as described
%    below
%    
%    This is intended for use with XSPDE jobs, and it has some features
%    to assist with running them in parallel.
%    
%    The parameters.m script should define the XSPDE input structures
%    that the job will solve, or at least the points, steps and ensembles
%    fields of those structures.  The script COMMAND can set up the rest
%    of those structures, then run parameters to assign these fields.
%    
%    FIRE generates a script ensembles.m, which the script parameters.m
%    should run as its final step.  This assigns values to the ensembles
%    fields of the structures listed in fire_inputs, which have appropriate
%    parallelism for QUEUE.  The total number of trajectories sampled
%    is kept constant, except that the generated value for ensembles(2)
%    is rounded up to an integer.  You should have a stub version of
%    ensembles.m to run on your desktop machine: an empty file will do.
%    
%    FIRE also generates a batch script called runjob, and calls qsub
%    to run it, and also The job control script will request parallel
%    Matlab resources consistent with the ensembles it generates.
%    
%    One of the parameters defined is hours_per_step. Run a small job
%    to start, then set this in parameters.m, so that fire will automatically
%    request an appropriate wall time.
%    
%    The latter assigns ensembles fields to the input structures in
%    fire_inputs, so that each input structure has the same number of
%    trajectories as requested, but with parallelism appropriate to
%    QUEUE. Your parameters.m script should run this as its last step.
%    The idea is to write parameters.m to run on your desktop machine,
%    and let the generated ensembles script reconfigure it for the
%    cluster.
%
%    Describe rounding up of wall time and ensembles
%
%    xpsetup_gstar.m gets copied to xpsetup.m so that you can have a desktop version in the same directory.

switch queue
	case 'gstar', cpus = 10;
	case 'sstar', cpus = 14;
	otherwise error 'Unknown job queue'
end

parameters
ssh_cmd = sprintf('ssh %s@g2.hpc.swin.edu.au ', gstar_login);

[status,~] = system([ssh_cmd 'ls -d ' letter]);
if status == 0
	error('Target directory %s already exits', letter)
else
	system([ssh_cmd 'mkdir -p ' letter]);
end

fire_files = [fire_files 'parameters.m'];
system(sprintf('scp -q %s %s@g2.hpc.swin.edu.au:%s', ...
		sprintf('%s ', fire_files{:}), gstar_login, letter));
if exist('xpsetup_gstar.m')
	system(sprintf('scp -q xpsetup_gstar.m %s@g2.hpc.swin.edu.au:%s/xpsetup.m', ...
		gstar_login, letter));
end

work = 0;
efilename = tempname;
ensfile = fopen(efilename, 'w');
for istruc = fire_inputs
	istruc = istruc{:};
	if exist(istruc)
		input = eval(istruc);
	end
	jobs = prod(input.ensembles(2:3));
	runs = ceil(jobs/cpus);
	fprintf(ensfile, '%s.ensembles = [%d %d %d];\n', ...
		istruc, input.ensembles(1), runs, cpus);
	work = work + runs*input.ensembles(1)*prod(input.points)*input.steps;
end
fclose(ensfile);
system(sprintf('scp -q %s %s@g2.hpc.swin.edu.au:%s/ensembles.m', ...
	efilename, gstar_login, letter));

runfilename = tempname;
runfile = fopen(runfilename, 'w');
fprintf(runfile, '#!/bin/sh\n#PBS -q %s\n', queue);
fprintf(runfile, '#PBS -l nodes=1:ppn=%d\n', cpus);
fprintf(runfile, '#PBS -l walltime=%d:00:00\n#PBS -N %s\n', ceil(2*work/1e9), letter);
fprintf(runfile, ['module load matlab/R2015b\ncd %s\n' ...
	'echo "Working directory: %s"\n' ...
	'matlab -r ''parpool(%d); xpsetup; %s; quit''\n'], ...
	['~/' letter], letter, cpus, command);
fclose(runfile);
system(sprintf('scp -q %s %s@g2.hpc.swin.edu.au:%s/runjob', ...
	runfilename, gstar_login, letter));

system(sprintf('%s "cd %s; /opt/torque/bin/qsub runjob" ', ssh_cmd, letter));

end