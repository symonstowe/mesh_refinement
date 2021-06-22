function node_balance_comp()
	% Generate and save all of the meshes for the analysis 
	% Starting divisions per electrode
	constant_div = 5;
	mdl_fname = ['FMDL',num2str(constant_div),'.mat'];
	ref_mdl_fname = 'FMDL_ref.mat'; % Use a separte file for the ref meshes 
	w_r = who('-file',ref_mdl_fname);
	% Generate the reference meshes
	ref_divs = 20; % These are HUGE meshes
	ref_size = 50/ref_divs/1000;
	% Netgen 
	mdl_name = 'mdl_000';
	if ismember(mdl_name,w_r)
		load(ref_mdl_fname,mdl_name);
		fmdl = eval(mdl_name);
		fmdl.name = mdl_name;	
	else
		fmdl = netgen_gen(ref_size,ref_size);
		fmdl.name = mdl_name;	
		eval([mdl_name ' = fmdl;']);
		save(ref_mdl_fname,'-append',mdl_name);
	end
	% GMSH
	mdl_name = 'mdl_100';
	if ismember(mdl_name,w_r)
		load(ref_mdl_fname,mdl_name);
		fmdl = eval(mdl_name);
		fmdl.name = mdl_name;	
	else
		fmdl = gmsh_gen(ref_size,ref_size,0,0,1);
		fmdl.name = mdl_name;	
		eval([mdl_name ' = fmdl;']);
		save(ref_mdl_fname,'-append',mdl_name);
	end
	w = who('-file',mdl_fname);
	mesh_mult = linspace(1,2,21); % Mesh multipliers for the electrode size 
	% make the constant mesh for each upon which everything else will be based
	keyboard






	for i = 2:length(mesh_mult) % Lets try 20 different mesh sizes
		% Netgen first 
		mdl_name = ['mdl_0' num2str(i,'%02.f')];
		if ismember(mdl_name,w)
			load(ref_mdl_fname,mdl_name);
			fmdl = eval(mdl_name);
			fmdl.name = mdl_name;
		else
			fmdl = netgen_gen(mesh_mult(i));
			fmdl.name = mdl_name;	
			eval([mdl_name ' = fmdl;']);
			save(mdl_fname,'-append',ref_mdl_name);
		end
		fmdl = get_mesh(mdl_name);

		% GMSH
		mdl_name = ['mdl_1' num2str(i,'%02.f')];
		fmdl = get_mesh(mdl_name);


	end
	
end

function fmdl = netgen_gen(e_maxh, g_maxh)
	% TODO are there more options that can be added here? Is there someone that would know?
	%opt.meshoptions.fineness = 6;
	%opt.options.meshsize = 0.05; % max elem size
	%opt.options.grading = 1; % how quickly refinement dissipates?!?!?
	%ng_write_opt(opt);
	fmdl = ng_mk_cyl_models([0.25,0.25,g_maxh],[4 0.125],[0.05,0.05,e_maxh]);
	fmdl = order_elecs(fmdl);
	fmdl.stimulation = mk_stim_patterns(4,1,'{ad}','{ad}',[],1);
	fmdl.stimulation = fmdl.stimulation(2);
end

function fmdl = gmsh_gen(e_maxh, g_maxh, r1, r2, o_flag)
	% Write the gmsh file and run it to create the mesh
	% o_flag specifies optimization technique 
	fid = fopen('temp_geo.geo','w');    
	fprintf(fid,'SetFactory("OpenCASCADE"); \n');
	fprintf(fid,'tank_boundary = newv; Cylinder(tank_boundary) = {0, 0, 0,  0, 0, 0.25, 0.25, 2*Pi}; \n');
	fprintf(fid,'Box(2) = {-0.025, 0.225,0.10, 0.05, 0.05,0.05}; \n');
	fprintf(fid,'Box(3) = { 0.025,-0.225,0.10,-0.05,-0.05,0.05}; \n');
	fprintf(fid,'Box(4) = { 0.225, 0.025,0.10, 0.05,-0.05,0.05}; \n');
	fprintf(fid,'Box(5) = {-0.225,-0.025,0.10,-0.05, 0.05,0.05}; \n');
	fprintf(fid,'elecs() = BooleanDifference{ Volume{2:5}; Delete;}{Volume{1};}; \n');
	fprintf(fid,'vol() = BooleanFragments{ Volume{tank_boundary};}{Volume{elecs()};}; \n');
	fprintf(fid,'Coherence; \n');
	fprintf(fid,'Mesh.CharacteristicLengthFromPoints = 0; \n');
	fprintf(fid,'Mesh.CharacteristicLengthFromCurvature = 0; \n');
	fprintf(fid,'Mesh.CharacteristicLengthExtendFromBoundary = 0; \n');
	fprintf(fid,'Field[1] = Attractor; \n');
	fprintf(fid,'Field[1].NNodesByEdge = 1000; \n');
	fprintf(fid,'Field[1].FacesList = {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28}; \n');
	fprintf(fid,'Field[2] = Threshold; \n');
	fprintf(fid,'Field[2].IField = 1; \n');
	fprintf(fid,'Field[2].LcMin = %f; \n', e_maxh);
	fprintf(fid,'Field[2].LcMax = %f; \n', g_maxh);
	fprintf(fid,'Field[2].DistMin = %f; \n', r1);
	fprintf(fid,'Field[2].DistMax = %f; \n', r2);
	fprintf(fid,'Background Field = 2; \n');
	fprintf(fid,'Mesh 1; \n');
	fprintf(fid,'Mesh 2; \n');
	fprintf(fid,'Mesh 3; \n');
	if o_flag == 1
	fprintf(fid,'OptimizeMesh "Gmsh"; \n');
	elseif o_flag == 2
	fprintf(fid,'OptimizeMesh "Netgen"; \n');
	else
	fprintf('No 3D mesh optimization \n')
	end
	while(1)
		status = system('gmsh temp_geo.geo -save');
		if status==0; break; end
	end
	[mdl,~] = gmsh_mk_fwd_model('temp_geo.msh', [], [], [], []);
	mdl.mat_idx_to_electrode.z_contact = 0.0100;
	mdl.mat_idx_to_electrode.nodes_electrode = true;
	for i=(2:5)
	mdl = mat_idx_to_electrode(mdl, {i});
	end
	fmdl = mdl;
	fmdl = order_elecs(fmdl);
	fmdl.stimulation = mk_stim_patterns(4,1,'{ad}','{ad}',[],1);
	fmdl.stimulation = fmdl.stimulation(2);	
	end

function fmdl = order_elecs(fmdl)
	new_elec_locs = [0,0.25,0.125;
			0.25,0,0.125;
			0,-.25,0.125;
			-.25,0,0.125];
	fmdl = renumber_electrodes( fmdl , new_elec_locs, 0);
end	
	