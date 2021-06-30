function node_balance_comp()
	% Input the starting divisions per electrode edge for m-type meshes
	starting_div = 3;
	[mesh_tab,ref_meshes,m_meshes] = gen_meshes(starting_div);
	[sens_tab] = calc_all_sens(ref_meshes, m_meshes, starting_div);
	mk_sens_plots(mesh_tab,sens_tab);
	mk_latex_tables(mesh_tab,sens_tab);
end

function [mesh_tab,ref_meshes,m_meshes] = gen_meshes(constant_div)
	% Generate and save all of the meshes for the analysis 
	% input is the starting divisions per electrode
	mdl_fname = ['FMDL',num2str(constant_div),'.mat'];
	ref_mdl_fname = 'FMDL_ref.mat'; % Use a separte file for the ref meshes 
	w_r = who('-file',ref_mdl_fname);
	% Generate the reference meshes
	ref_divs = 15; % These are HUGE meshes %TODO can we use a size 25 divisions?
	ref_size = 50/ref_divs/1000;
	% Netgen 
	mdl_name = 'mdl_000';
	if ismember(mdl_name,w_r)
		load(ref_mdl_fname,mdl_name);
		fmdl = eval(mdl_name);
		fmdl.mesh_params.spec_min = ref_size*1000;
		fmdl.mesh_params.spec_max = ref_size*1000;
		eval([mdl_name ' = fmdl;']);
	else
		fmdl = netgen_gen(ref_size,ref_size);
		fmdl.name = mdl_name;	
		fmdl = calc_mesh_info(fmdl);
		eval([mdl_name ' = fmdl;']);
		save(ref_mdl_fname,'-append',mdl_name,'-v7.3');
	end
	% GMSH
	mdl_name = 'mdl_100';
	if ismember(mdl_name,w_r)
		load(ref_mdl_fname,mdl_name);
		fmdl = eval(mdl_name);
		fmdl.mesh_params.spec_min = ref_size*1000;
		fmdl.mesh_params.spec_max = ref_size*1000;
		eval([mdl_name ' = fmdl;']);
	else
		fmdl = gmsh_gen(ref_size,ref_size,0,0,1,0);
		fmdl.name = mdl_name;	
		fmdl = calc_mesh_info(fmdl);
		eval([mdl_name ' = fmdl;']);
		save(ref_mdl_fname,'-append',mdl_name,'-v7.3');
	end
	% Generate the reference parmeter table
	mesh_tab0  = struct2table(mdl_000.mesh_params);
	mesh_tab0  = [table({mdl_000.name},'VariableNames',{'Name'}), mesh_tab0] ;
	mesh_tab1  = struct2table(mdl_100.mesh_params);
	mesh_tab1  = [table({mdl_100.name},'VariableNames',{'Name'}), mesh_tab1] ;
	mesh_tab   = [mesh_tab0; mesh_tab1];
	% Generate the meshes for analysis
	w = who('-file',mdl_fname);
	mesh_mult = linspace(1,0.2,17); % Mesh multipliers for the electrode size 
	% make the constant mesh for each upon which everything else will be based
	base_size = 50/constant_div/1000;
	% Netgen
	mdl_name = 'mdl_001';
	if ismember(mdl_name,w)
		load(mdl_fname,mdl_name);
		fmdl = eval(mdl_name);
		fmdl.mesh_params.spec_min = base_size*1000;
		fmdl.mesh_params.spec_max = base_size*1000;
		eval([mdl_name ' = fmdl;']);
	else
		fmdl = netgen_gen(base_size,base_size);
		fmdl.name = mdl_name;	
		fmdl = calc_mesh_info(fmdl);
		eval([mdl_name ' = fmdl;']);
		save(mdl_fname,'-append',mdl_name);
	end
	n_ref_nodes = fmdl.mesh_params.num_nodes;
	mesh_tab_temp = struct2table(eval([mdl_name '.mesh_params']));
	mesh_tab_temp = [table({eval([mdl_name '.name'])},'VariableNames',{'Name'}), mesh_tab_temp];
	mesh_tab = [mesh_tab; mesh_tab_temp];
	% GMSH
	mdl_name = 'mdl_101';
	if ismember(mdl_name,w)
		load(mdl_fname,mdl_name);
		fmdl = eval(mdl_name);
		fmdl.mesh_params.spec_min = base_size*1000;
		fmdl.mesh_params.spec_max = base_size*1000;
		eval([mdl_name ' = fmdl;']);
	else
		fmdl = gmsh_gen(base_size,base_size,0,0,1,0);
		fmdl.name = mdl_name;	
		fmdl = calc_mesh_info(fmdl);
		eval([mdl_name ' = fmdl;']);
		save(mdl_fname,'-append',mdl_name);
	end
	g_ref_nodes = fmdl.mesh_params.num_nodes;
	mesh_tab_temp = struct2table(eval([mdl_name '.mesh_params']));
	mesh_tab_temp = [table({eval([mdl_name '.name'])},'VariableNames',{'Name'}), mesh_tab_temp];
	mesh_tab = [mesh_tab; mesh_tab_temp];
	% Run through the non-contant meshes using the reference mesh node numbers to truncate the meshing...
	for i = 2:length(mesh_mult) % Lets try 20 different mesh sizes
		% Netgen first 
		mdl_name = ['mdl_0' num2str(i,'%02.f')];
		if ismember(mdl_name,w)
			load(mdl_fname,mdl_name);
			fmdl = eval(mdl_name);
			fmdl.name = mdl_name;
		else
			num_nodes = 0;
			node_adjust = 0;
			while(abs(num_nodes-n_ref_nodes)>n_ref_nodes*0.15)
				fmdl = netgen_gen(mesh_mult(i)*base_size,(1-mesh_mult(i)+1)*base_size+(node_adjust*base_size));
				num_nodes = size(fmdl.nodes,1);
				max_size = (1-mesh_mult(i)+1)*base_size+(node_adjust*base_size);
				min_size = mesh_mult(i)*base_size;
				if num_nodes>n_ref_nodes 
					node_adjust = node_adjust+0.05;
				else
					node_adjust = node_adjust-0.1;
				end
			end
			fmdl.mesh_params.spec_min = min_size*1000;
			fmdl.mesh_params.spec_max = max_size*1000;
			fmdl.name = mdl_name;	
			fmdl = calc_mesh_info(fmdl);
			eval([mdl_name ' = fmdl;']);
			save(mdl_fname,'-append',mdl_name);
		end
		% Add to the table
		mesh_tab_temp = struct2table(eval([mdl_name '.mesh_params']));
		mesh_tab_temp = [table({eval([mdl_name '.name'])},'VariableNames',{'Name'}), mesh_tab_temp];
		mesh_tab = [mesh_tab; mesh_tab_temp];
		% GMSH
		mdl_name = ['mdl_1' num2str(i,'%02.f')];
		if ismember(mdl_name,w)
			load(mdl_fname,mdl_name);
			fmdl = eval(mdl_name);
			fmdl.name = mdl_name;
		else
			num_nodes = 0;
			node_adjust = 0;
			while(abs(num_nodes-g_ref_nodes)>g_ref_nodes*0.15)
				fmdl = gmsh_gen(mesh_mult(i)*base_size,(1-mesh_mult(i)+1)*base_size+(node_adjust*base_size),0,0.25,1,0);
				num_nodes = size(fmdl.nodes,1);
				max_size = (1-mesh_mult(i)+1)*base_size+(node_adjust*base_size);
				min_size = mesh_mult(i)*base_size;
				if num_nodes>g_ref_nodes 
					node_adjust = node_adjust+0.05;
				else
					node_adjust = node_adjust-0.1;
				end
			end
			fmdl.mesh_params.spec_min = min_size*1000;
			fmdl.mesh_params.spec_max = max_size*1000;
			fmdl.name = mdl_name;	
			fmdl = calc_mesh_info(fmdl);
			eval([mdl_name ' = fmdl;']);
			save(mdl_fname,'-append',mdl_name);
		end
		% Add to the table
		mesh_tab_temp = struct2table(eval([mdl_name '.mesh_params']));
		mesh_tab_temp = [table({eval([mdl_name '.name'])},'VariableNames',{'Name'}), mesh_tab_temp];
		mesh_tab = [mesh_tab; mesh_tab_temp];
		% GMSH 2!!! only at the edges!
		mdl_name = ['mdl_2' num2str(i,'%02.f')];
		if ismember(mdl_name,w)
			load(mdl_fname,mdl_name);
			fmdl = eval(mdl_name);
			fmdl.name = mdl_name;
		else
			num_nodes = 0;
			node_adjust = 0;
			while(abs(num_nodes-g_ref_nodes)>g_ref_nodes*0.15)
				fmdl = gmsh_gen(mesh_mult(i)*base_size,(1-mesh_mult(i)+1)*base_size+(node_adjust*base_size),0,0.25,1,1);
				num_nodes = size(fmdl.nodes,1);
				max_size = (1-mesh_mult(i)+1)*base_size+(node_adjust*base_size);
				min_size = mesh_mult(i)*base_size;
				if num_nodes>g_ref_nodes 
					node_adjust = node_adjust+0.05;
				else
					node_adjust = node_adjust-0.1;
				end
			end
			fmdl.mesh_params.spec_min = min_size*1000;
			fmdl.mesh_params.spec_max = max_size*1000;
			fmdl.name = mdl_name;	
			fmdl = calc_mesh_info(fmdl);
			eval([mdl_name ' = fmdl;']);
			save(mdl_fname,'-append',mdl_name);
		end
		% Add to the table
		mesh_tab_temp = struct2table(eval([mdl_name '.mesh_params']));
		mesh_tab_temp = [table({eval([mdl_name '.name'])},'VariableNames',{'Name'}), mesh_tab_temp];
		mesh_tab = [mesh_tab; mesh_tab_temp];
	end
	ref_meshes = load(ref_mdl_fname);
	m_meshes = load(mdl_fname);
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

function fmdl = gmsh_gen(e_maxh, g_maxh, r1, r2, o_flag, edge_flag)
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
	if edge_flag == 1
		fprintf(fid,'Field[1].EdgesList = {4,5,6,7,8,9,11,12,13,15,16,17,18,19,21,22,23,25,26,27,28,29,30,31,32,33,36,37,38,41,42,43,44,45,46,47,48,49,50,52,53,55,56}; \n');
	else
		fprintf(fid,'Field[1].FacesList = {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28}; \n');
	end
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

function fmdl = calc_mesh_info(fmdl)
	% Use fix model to get advanced metrics
	temp_fmdl = fix_model(fmdl);
	% State Units (could be helpful later)
	fmdl.mesh_params.units = 'mm or mm^3';
	% Calculate:
	% 1) Number of elements
	fmdl.mesh_params.num_elems = size(fmdl.elems,1);
	% 2) Number of nodes
	fmdl.mesh_params.num_nodes = size(fmdl.nodes,1);
	% 3) Elements per electrode (average?)
	elec_faces = ismember(fmdl.elems, fmdl.electrode(1).nodes);
    	elec_face_num = length(find(sum(elec_faces,2) == 3));
	if elec_face_num == 0
		elec_face_num = length(fmdl.electrode(1).faces);
	end  
	fmdl.mesh_params.elems_per_elec = elec_face_num;
	% 4) Max edge length
	fmdl.mesh_params.max_el = max(temp_fmdl.edge_length)*1000;
	% 5) Min Edge length
	fmdl.mesh_params.min_el = min(temp_fmdl.edge_length)*1000;
	% 6) Max element volume
	fmdl.mesh_params.max_ev = max(temp_fmdl.elem_volume)*1e+9;
	% 7) Min element volume
	fmdl.mesh_params.min_ev = min(temp_fmdl.elem_volume)*1e+9;
	% 8) Center of mass (for the mesh analysis) (As a percent!)
	[balance_pt,~,~] = node_balance(fmdl);
	fmdl.mesh_params.balance_pt = balance_pt;
end

function [balance_pt,count_nodes,I] = node_balance(fmdl)
	% input a FMDL 
	% ouput the node balance in percent of model and 
	% number of nodes in the selected area between the electroeso
	fmdl = order_elecs(fmdl);
	elec_nodes = fmdl.nodes(fmdl.electrode(1).nodes,:);
	edge_nodes = [fmdl.electrode(1).nodes,fmdl.electrode(3).nodes];
	max_z = max(elec_nodes(:,3));
	min_z = min(elec_nodes(:,3));
	max_wdth = max(elec_nodes(:,1));
	min_wdth = min(elec_nodes(:,1));
	I = find((fmdl.nodes(:,1) <= max_wdth) & (fmdl.nodes(:,1) >= min_wdth) & ...
		 (fmdl.nodes(:,3) <= max_z)    & (fmdl.nodes(:,3) >= min_z) & ...
		 (fmdl.nodes(:,2) > 0) );
	
	I = setdiff(I,edge_nodes); % Remove the nodes that are in the electrode suraface positions
	count_nodes = size(I);
	balance_pt = mean(fmdl.nodes(I,2))/0.25*100;
end

function [sens_tab] = calc_all_sens(ref_meshes, m_meshes, starting_div)
	% Ouput the sensitivity error for all regions including the next-to-electrode
	% Table format
	sens_tab=cell2table(cell(0,7),'VariableNames',{'Name','se','si','me','mi','c','t'});
	fn = fieldnames(ref_meshes);
	sens = [];
	sns_fname = ['SENS',num2str(starting_div),'.mat'];
	ref_sns_fname = 'SENS_ref.mat';
	% Reference meshes first
	for i=1:numel(fn)
		w = who('-file',ref_sns_fname);
		if ismember([fn{i}, '_sens'],w)
			load(ref_sns_fname,[fn{i}, '_sens']);	
			%sens_img = eval([fn{i}, '_sens']);
		else
			fmdl =  ref_meshes.(fn{i});
			eval([fn{i}, '_sens = calc_sens(fmdl);']);
			save('SENS_ref.mat','-append',[fn{i}, '_sens']);
		end
	end
	fn = fieldnames(m_meshes);
	for i=1:numel(fn)
		var_name = [fn{i}, '_sens'];
		d = regexp(var_name,'\d*','Match');
		if str2double(d{1}) > 100
			ref_sens = mdl_100_sens;
		else
			ref_sens = mdl_000_sens;
		end
		w = who('-file',sns_fname);
		if ismember(var_name,w)
			load(sns_fname,var_name);	
		else
			fmdl =  m_meshes.(fn{i});
			%sens_img = calc_sens(fmdl);
			eval([fn{i}, '_sens = calc_sens(fmdl);']);
			save(['SENS',num2str(starting_div),'.mat'],'-append',var_name);
		end
		sens_img = eval(var_name);
		[sens.se,sens.si,sens.me,sens.mi,sens.c,sens.t] = ... 
						compare_sens(sens_img,ref_sens);
		sens_tab_temp = struct2table(sens);
		sens.Name = fn{i};
		sens_tab = [sens_tab; struct2table(sens)];
	end
end
	
function roi =  get_roi(str)
		roi = false(512,512);
		switch(str)
		case 'se'
		    roi(462:512,193:320)=true;
		case 'si'
		    roi(410:461,193:320)=true;
		case 'me'
		    roi(192:320,1:51)=true;
		case 'mi'
		    roi(192:320,52:103)=true;
		case 'c'
		    roi(192:320,192:320)=true;
		end
end
	    
function [imgv vh] = solve_for_volts(fmdl)
	img = mk_image(fmdl,1);
	img.fwd_solve.get_all_meas = 1;
	vh = fwd_solve(img);
	imgv= rmfield(img,'elem_data');
	imgv.node_data = vh.volt;
	imgv.calc_colours.npoints = 512;
	vh = vh.meas;
end
	    
function rst = volt2raster(imgv)
	rst = calc_slices(imgv,[inf inf 0]);
end
	    
function sens = calc_sens(fmdl)
	%Given the FMDL return the sens values...
	img = mk_image(fmdl,1);
	J = calc_jacobian(img);
	J = jacobian_adjoint(img);
	S = J'./get_elem_volume(fmdl);
	tmp = img;
	tmp.elem_data = S;
	tmp.calc_colours.npoints = 512;
	Se = calc_slices(tmp,[inf inf 0.25/2]);
	ca = max(abs(Se(:)));
	
	N = 15;
	L = zeros(N,3);
	L(:,1:2) = inf;
	L(:,3) = linspace(0.05,0.2,N); % What should the seperation from the edges be?
	t = calc_slices(tmp,L);
	t = squeeze(t);
	Sa = mean(t,3);
	sens = Sa;
end
	
function [se,si,me,mi,c,t] = compare_sens(Sa,ref_Sa)
	regions = {'se', 'si', 'me', 'mi', 'c'};
	for j=0:length(regions)
		if j == 0
		roi = true(512,512);
		else
		roi = get_roi(regions{j});
		end
		ref_sa_roi = ref_Sa .* roi;
		sa_roi     = Sa .* roi;
		err = ref_sa_roi - sa_roi;
		err = sum(abs(err(~isnan(err))))/sum(sum(~isnan(err))); % What exactly do you call this?
		if j == 0
		t = err;
		else
		eval([regions{j},' = err;'])
		end
	end
end

function mk_latex_tables(mesh_tab,sens_tab)
	% Generate and format the tables
	% This is not a smart function... but the table looks nice
	fid = fopen('../docs/mesh_table.tex','w');    
	%fprintf(fid,'\\begin{table}[]\n');
	%fprintf(fid,'\\caption{\\label{tab:mesh-table}Your caption.}\n');
	fprintf(fid,'\\begin{tabular}{p{0.7cm}p{0.2cm}p{1cm}p{1cm}|p{1.4cm}p{1.4cm}p{1cm}|p{1.3cm}p{1.3cm}p{1.3cm}p{1.3cm}}\n');
	fprintf(fid,'\\multicolumn{2}{l}{Mesh ID} & ');
	fprintf(fid,'glbl. maxh [mm] & elec. maxh [mm] & ');
	fprintf(fid,'\\#~elem. & \\#~nodes & \\#~elec. elem. & ');
	fprintf(fid,'minEL\\textsuperscript{\\emph{a}} [mm] & maxEL\\textsuperscript{\\emph{b}} [mm] & ');
	fprintf(fid,'minEV\\textsuperscript{\\emph{c}} [mm\\textsuperscript{3}] & maxEV\\textsuperscript{\\emph{d}} [mm\\textsuperscript{3}]');
	fprintf(fid,'%s \\hline \n','\\');
	% Reference meshes
	for i=[linspace(1,17,9) 0]
		m_id = num2str(i,'%02.f');
		if i==0
			fprintf(fid,'\\multirow{2}{*}{REF} &');
		else
			fprintf(fid,'\\multirow{2}{*}{M-%s} &',m_id);
		end
		loc = find(strcmp(['mdl_0' m_id],mesh_tab.Name));
		fprintf(fid,'A & %0.2f & %0.2f & %d & %d & %d & %0.2f & %0.2f & %0.2f & %0.2f', ...
			mesh_tab.spec_max(loc),mesh_tab.spec_min(loc),mesh_tab.num_elems(loc), ...
			mesh_tab.num_nodes(loc),mesh_tab.elems_per_elec(loc), ... 
			mesh_tab.min_el(loc),mesh_tab.max_el(loc), ... 
			mesh_tab.min_ev(loc),mesh_tab.max_ev(loc));
		fprintf(fid,'%s \n','\\');
		loc = find(strcmp(['mdl_1' m_id],mesh_tab.Name));
		fprintf(fid,'& ');
		fprintf(fid,'B & %0.2f & %0.2f & %d & %d & %d & %0.2f & %0.2f & %0.2f & %0.2f', ...
			mesh_tab.spec_max(loc),mesh_tab.spec_min(loc),mesh_tab.num_elems(loc), ...
			mesh_tab.num_nodes(loc),mesh_tab.elems_per_elec(loc), ... 
			mesh_tab.min_el(loc),mesh_tab.max_el(loc), ... 
			mesh_tab.min_ev(loc),mesh_tab.max_ev(loc));
		if i==17
			fprintf(fid,'%s \\hline \n','\\');
		else
			fprintf(fid,'%s \n','\\');
		end
	end
	fprintf(fid,'\\multicolumn{11}{c}{\\emph{a}: minimum mesh edge length, \\emph{b}: maximum mesh edge length}');
	fprintf(fid,'%s \n','\\');
	fprintf(fid,'\\multicolumn{11}{c}{\\emph{c}: minimum mesh element volume, \\emph{d}: maximum mesh element volume}');
	fprintf(fid,'%s \n','\\');
	fprintf(fid,'\\end{tabular}\n');
	%fprintf(fid,'\\end{table}\n');
end

function mk_sens_plots_log(mesh_tab,sens_tab)
	% For every row in the sens tab find the correspoonding row in the mesh tab
	% Grab the CofM, and plot the sens error!	
	sens_names = cellstr(sens_tab.Name);
	for i=1:size(mesh_tab,1)
		var_name = mesh_tab.Name{i};
		d = regexp(var_name,'\d*','Match');
		d = str2double(d);
		if (d > 100) & (d < 200)
			loc = find(strcmp(var_name,sens_names));
			gmsh_sens(i) = sens_tab.t(loc);
			gmsh_sens_e(i) = sens_tab.se(loc)+sens_tab.me(loc);
			gmsh_sens_i(i) = sens_tab.si(loc)+sens_tab.mi(loc);
			gmsh_sens_c(i) = sens_tab.c(loc);
			gmsh_CofM(i) = mesh_tab.balance_pt(i);
		elseif (d>0) & (d<100)
			loc = find(strcmp(var_name,sens_names));
			net_sens(i) = sens_tab.t(loc);
			net_sens_e(i) = sens_tab.se(loc)+sens_tab.me(loc);
			net_sens_i(i) = sens_tab.si(loc)+sens_tab.mi(loc);
			net_sens_c(i) = sens_tab.c(loc);
			net_CofM(i) = mesh_tab.balance_pt(i);
		%elseif (d>200) % Do the advanced advanced gmsh
			loc = find(strcmp(var_name,sens_names));
			gmsh_sens_2(i) = sens_tab.t(loc);
			gmsh_CofM_2(i) = mesh_tab.balance_pt(i);
		end
	end
	%keyboard
	net_sens = net_sens(net_sens > 0);
	net_sens_c = net_sens_c(net_sens_c > 0);
	net_sens_i = net_sens_i(net_sens_i > 0);
	net_sens_e = net_sens_e(net_sens_e > 0);
	net_CofM = net_CofM(net_CofM > 0);
	gmsh_sens = gmsh_sens(gmsh_sens > 0);
	gmsh_sens_e = gmsh_sens_e(gmsh_sens_e > 0);
	gmsh_sens_i = gmsh_sens_i(gmsh_sens_i > 0);
	gmsh_sens_c = gmsh_sens_c(gmsh_sens_c > 0);
	gmsh_CofM = gmsh_CofM(gmsh_CofM > 0);

	%gmsh_sens_2 = gmsh_sens_2(gmsh_sens_2 > 0);
	%gmsh_CofM_2 = gmsh_CofM_2(gmsh_CofM_2 > 0);
	figure(1); clf;
	% Do the netgen sensitivity error plots (in blue)
	[sorted_CofM, I] = sort(net_CofM);
	sorted_sens = net_sens(I);
	n1 = plot(sorted_CofM,sorted_sens,'LineWidth',3,'Color',[  3  78 123]/255);
	hold on
	sorted_sens = net_sens_e(I);
	n2 = plot(sorted_CofM,sorted_sens,'LineWidth',2,'Color',[ 54 144 192]/255);
	sorted_sens = net_sens_i(I);
	n3 = plot(sorted_CofM,sorted_sens,'LineWidth',2,'Color',[116 169 207]/255);
	sorted_sens = net_sens_c(I);
	n4 = plot(sorted_CofM,sorted_sens,'LineWidth',2,'Color',[166 189 219]/255);

	g1 = plot(gmsh_CofM,gmsh_sens  ,'LineWidth',3,'Color',[153   0   0]/255);
	g2 = plot(gmsh_CofM,gmsh_sens_e,'LineWidth',2,'Color',[239 101  72]/255);
	g3 = plot(gmsh_CofM,gmsh_sens_i,'LineWidth',2,'Color',[252 141  89]/255);
	g4 = plot(gmsh_CofM,gmsh_sens_c,'LineWidth',2,'Color',[253 187 132]/255);
	%p2 = plot(gmsh_CofM_2,gmsh_sens_2,'LineWidth',2,'Color',[77 175 74]/255);
	hold off
	set(groot,'defaultAxesTickLabelInterpreter','latex');  
	set(groot,'defaulttextinterpreter','latex');
	set(groot,'defaultLegendInterpreter','latex');
	set(gca, 'YScale', 'log')
	set(get(gca, 'XLabel'), 'String', 'Node balance point');
	set(get(gca, 'YLabel'), 'String', 'Sensitivity Error');
	legend('Total Netgen', 'Region E Netgen','Region I Netgen', 'Region C Netgen', ...
               'Total Gmsh'  , 'Region E Gmsh'  ,'Region I Gmsh'  , 'Region C Gmsh'  , ...
           'Interpreter','latex','NumColumns',2,'FontSize',15,'Location','NorthEast');
	set(gcf,'Position',[103 584 1399 820])
	%legend('GMSH electrode refinement','Netgen electrode refinement');
	set(gca,'YGrid','on');
	print('../imgs/m-mesh_sens_error_log', '-dsvg')
end

function mk_sens_plots(mesh_tab,sens_tab)
	% For every row in the sens tab find the correspoonding row in the mesh tab
	% Grab the CofM, and plot the sens error!	
	sens_names = cellstr(sens_tab.Name);
	for i=1:size(mesh_tab,1)
		var_name = mesh_tab.Name{i};
		d = regexp(var_name,'\d*','Match');
		d = str2double(d);
		if (d > 100) & (d < 200)
			loc = find(strcmp(var_name,sens_names));
			gmsh_sens(i) = sens_tab.t(loc);
			gmsh_sens_e(i) = sens_tab.se(loc)+sens_tab.me(loc);
			gmsh_sens_i(i) = sens_tab.si(loc)+sens_tab.mi(loc);
			gmsh_sens_c(i) = sens_tab.c(loc);
			gmsh_CofM(i) = mesh_tab.balance_pt(i);
		elseif (d>0) & (d<100)
			loc = find(strcmp(var_name,sens_names));
			net_sens(i) = sens_tab.t(loc);
			net_sens_e(i) = sens_tab.se(loc)+sens_tab.me(loc);
			net_sens_i(i) = sens_tab.si(loc)+sens_tab.mi(loc);
			net_sens_c(i) = sens_tab.c(loc);
			net_CofM(i) = mesh_tab.balance_pt(i);
		%elseif (d>200) % Do the advanced advanced gmsh
			loc = find(strcmp(var_name,sens_names));
			gmsh_sens_2(i) = sens_tab.t(loc);
			gmsh_CofM_2(i) = mesh_tab.balance_pt(i);
		end
	end
	%keyboard
	net_sens = net_sens(net_sens > 0);
	net_sens_c = net_sens_c(net_sens_c > 0);
	net_sens_i = net_sens_i(net_sens_i > 0);
	net_sens_e = net_sens_e(net_sens_e > 0);
	net_CofM = net_CofM(net_CofM > 0);
	gmsh_sens = gmsh_sens(gmsh_sens > 0);
	gmsh_sens_e = gmsh_sens_e(gmsh_sens_e > 0);
	gmsh_sens_i = gmsh_sens_i(gmsh_sens_i > 0);
	gmsh_sens_c = gmsh_sens_c(gmsh_sens_c > 0);
	gmsh_CofM = gmsh_CofM(gmsh_CofM > 0);

	%gmsh_sens_2 = gmsh_sens_2(gmsh_sens_2 > 0);
	%gmsh_CofM_2 = gmsh_CofM_2(gmsh_CofM_2 > 0);
	figure(1); clf;
	tiledlayout(4,1, 'Padding', 'none', 'TileSpacing', 'compact'); 
	set(groot,'defaultAxesTickLabelInterpreter','latex');  
	set(groot,'defaulttextinterpreter','latex');
	set(groot,'defaultLegendInterpreter','latex');
	% Do the netgen sensitivity error plots (in blue)
	[sorted_CofM, I] = sort(net_CofM);
	sorted_sens = net_sens(I);
	nexttile
	n1 = plot(sorted_CofM,sorted_sens,'LineWidth',3,'Color',[  3  78 123]/255);
	hold on
	g1 = plot(gmsh_CofM,gmsh_sens  ,'LineWidth',3,'Color',[153   0   0]/255);
	grid on	
	set(get(gca, 'XLabel'), 'String', 'Node balance point as a percentage of radial distance');
	set(get(gca, 'YLabel'), 'String', 'Sensitivity Error');
	set(gca,'FontSize',15)
	legend('Total Netgen', 'Total Gmsh', ... 
           'Interpreter','latex','NumColumns',2,'FontSize',15,'Location','NorthEast');
	hold off 
	nexttile
	sorted_sens = net_sens_e(I);
	n2 = plot(sorted_CofM,sorted_sens,'LineWidth',2,'Color',[ 54 144 192]/255);
	hold on
	g2 = plot(gmsh_CofM,gmsh_sens_e,'LineWidth',2,'Color',[239 101  72]/255);
	grid on	
	set(get(gca, 'XLabel'), 'String', 'Node balance point as a percentage of radial distance');
	set(get(gca, 'YLabel'), 'String', 'Sensitivity Error');
	set(gca,'FontSize',15)
	legend('Region E Netgen', 'Region E Gmsh', ... 
           'Interpreter','latex','NumColumns',2,'FontSize',15,'Location','NorthEast');
	hold off
	nexttile
	sorted_sens = net_sens_i(I);
	n3 = plot(sorted_CofM,sorted_sens,'LineWidth',2,'Color',[116 169 207]/255);
	hold on
	g3 = plot(gmsh_CofM,gmsh_sens_i,'LineWidth',2,'Color',[252 141  89]/255);
	grid on	
	set(get(gca, 'XLabel'), 'String', 'Node balance point as a percentage of radial distance');
	set(get(gca, 'YLabel'), 'String', 'Sensitivity Error');
	set(gca,'FontSize',15)
	legend('Region I Netgen', 'Region I Gmsh', ... 
           'Interpreter','latex','NumColumns',2,'FontSize',15,'Location','NorthEast');
	hold off
	nexttile
	sorted_sens = net_sens_c(I);
	n4 = plot(sorted_CofM,sorted_sens,'LineWidth',2,'Color',[166 189 219]/255);
	hold on
	g4 = plot(gmsh_CofM,gmsh_sens_c,'LineWidth',2,'Color',[253 187 132]/255);
	grid on	
	set(get(gca, 'XLabel'), 'String', 'Node balance point as a percentage of radial distance');
	set(get(gca, 'YLabel'), 'String', 'Sensitivity Error');
	set(gca,'FontSize',15)
	legend('Region C Netgen', 'Region C Gmsh', ... 
           'Interpreter','latex','NumColumns',2,'FontSize',15,'Location','NorthEast');
	hold off
	%set(get(gca, 'XLabel'), 'String', 'Node balance point as a percentage of radial distance');
	%set(get(gca, 'YLabel'), 'String', 'Sensitivity Error');
	%legend('Total Netgen', 'Region E Netgen','Region I Netgen', 'Region C Netgen', ...
 %              'Total Gmsh'  , 'Region E Gmsh'  ,'Region I Gmsh'  , 'Region C Gmsh'  , ...
 %          'Interpreter','latex','NumColumns',2,'FontSize',15,'Location','NorthEast');
	set(gcf,'Position',[606 11 1158 1311])
	%legend('GMSH electrode refinement','Netgen electrode refinement');
	%set(gca,'YGrid','on');
	print('../imgs/m-mesh_sens_error_regions_split', '-dsvg')
end