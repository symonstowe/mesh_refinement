function node_balance_comp()
	% Generate and save all of the meshes for the analysis 
%	if ~isfile('FMDL.mat')
%		% File does not exist. - clear all and save a blank file 
%		clear
%		save FMDL.mat
%	end
	constant_div = 5;
	mdl_fname = ['FMDL',num2str(constant_div),'.mat'];
	w = who('-file',mdl_fname);
	% Starting divisions per electrode
	% Generate the reference meshes
	ref_divs = 20;
	ref_size = 50/ref_divs/1000;
	% Netgen 
	mdl_name = 'mdl_000';
	if ismember(mdl_name,w)
		load(mdl_fname,mdl_name);
		fmdl = eval(mdl_name);
		fmdl.name = mdl_name;	
	else
		fmdl = netgen_gen()
	end
	% GMSH
	med_name = 'mdl_100';
	fmdl_ref1 = get_mesh(mdl_name);

	mesh_mult = linspace(1,2,21); % Mesh multipliers for the electrode size 
	for i = 1:length(mesh_mult) % Lets try 20 different mesh sizes (plus one reference!)
		% Netgen first 
		mdl_name = ['mdl_0' num2str(i,'%02.f')];
		fmdl = get_mesh(mdl_name);

		% GMSH
		mdl_name = ['mdl_1' num2str(i,'%02.f')];
		fmdl = get_mesh(mdl_name);


	end
	
end

function fmdl = netgen_gen(e_maxh, g_maxh)
	%opt.meshoptions.fineness = 6;
	%opt.options.meshsize = 0.05; % max elem size
	%opt.options.grading = 1; % how quickly refinement dissipates?!?!?
	%ng_write_opt(opt);
	fmdl = ng_mk_cyl_models([0.25,0.25,g_maxh],[4 0.125],[0.05,0.05,e_maxh]);
    end
	