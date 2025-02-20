	
	# Run this octave script in the directory containing the following files:

	dir_name = "./"
	
	filename = strcat(dir_name,"/heis_map_vecs.m")
	load(filename);
	filename = strcat(dir_name,"/heis_map_energies.m")
	load(filename);
	filename = strcat(dir_name,"/heis_map_sort_ind.m")
	load(filename);
	filename = strcat(dir_name,"/heis_map_proj_weights.m")
	load(filename);
	
	au2ev = 27.21165;
	au2cm = 219474.63;
	units = "cm-1";
	convert = au2cm;

	sort_ind	= heis_map_sort_ind;
	
	# Count size of blocks
	tmp1 = max(sort_ind)+1;
	count_vec = zeros(1,tmp1);
	for i = 1:length(sort_ind)
		ii = sort_ind(i)+1;
		count_vec(ii) += 1;
	end
	blocks = [];
       	for i = 1:length(count_vec)
		if count_vec(i) > 0
			blocks = [blocks, count_vec(i)];
		end
	end
	blocks

	# Should each orbital be counted as a separate site?
	sing_orb_site = 0
	if sing_orb_site
		blocks = ones(1,length(sort_ind))
	end

	proj_weights   	= heis_map_proj_weights;
	dim_heis	= sum(blocks)
	nstates		= columns(blocks)


	
	#
	#	Sort by projection onto model space
	#		However, we must also decide how many states for each site we need (todo)
       	#		Basically, I had the problem that 2 states on a single site had larger 
	#		projection than either from a different site, so it was knocked out and nothing made sense.
	#		the projection ordered approach will be necessary, but must be thought through more carefully.
	#
	#		turning off for now.
	sort_by_proj = 0
	if sort_by_proj
		[uselessVariable,permutation]=sort(proj_weights);
		permutation 	= flipud(permutation);
		#permutation 	= permutation(1:length(heis_map_energies));
		energies 	= heis_map_energies(permutation);
		vectors		= heis_map_vecs(:,permutation);
		for s = 1:nstates
			printf(" state: %5i  site: %5i   proj: %5f\n",permutation(s),sort_ind(s),proj_weights(s));
		endfor
		printf(" Sort vectors from high to low projections\n"); 
		disp(permutation(1:nstates));
		disp(proj_weights(1:nstates));
	end
	energies 	= heis_map_energies;
	vectors		= heis_map_vecs;

	
	printf("\n");
	printf("State energies\n");
	# shift energies for convenience 
	e = energies(1:nstates) +- min(energies);
	e = e +- max(e) - 1;
	for i=1:nstates
		printf("      State %4i: %12.8f\n",i,e(i));
	end
	
	
	# Sort determinant basis to order site-by-site
	printf("\n");
	printf("Projected eigenvectors\n");
	c = vectors(:,1:nstates);
	# sort determinant indices
	[uselessVariable,permutation]=sort(sort_ind);
	c = c(permutation,:);
	
	

	
	printf("\n");
	printf("Norms of projected eigenvectors (Projection #1 onto neutral determinant basis)\n");
	for s = 1:nstates
		printf("      State %4i: %12.8f\n",i,norm(c(:,s)));
	endfor

	# project onto local high-spin state
	c_new = zeros(length(blocks),nstates);
	first=1;
	bi = 1;
	for b = blocks
		last = first + b - 1;
		row = zeros(1,nstates);
		for i=first:last
			row += c(i,:);
		end
		c_new(bi,:) = row/sqrt(b);
		i = first;
		j = last;
		first = last + 1;
		bi += 1;
	end

	printf("\n");
	printf("Norms of projected eigenvectors (Projection #2 onto local high-spin states)\n");
	c = c_new;
	for s = 1:nstates
		printf("      State %4i: %12.8f\n",i,norm(c(:,s)));
	endfor
	
	printf("\n");
	printf("Overlap of projected eigenvectors\n");
	S = c'*c;
	disp(S);
	[U,l] = eig(S);
	X 	= U * l^(-1/2) * U';
	Xinv 	= U * l^(1/2) * U';
	C = c*X;
	
	printf("\n");
	printf("Orthogonalized eigenvectors\n");
	disp(C);
	
	nblocks = length(blocks)
	nstates
	if(length(blocks) != nstates)
		printf("Number of blocks not equal to number of states: die.\n")
		return;
	end

	printf("\n");
	printf("Effective Hamiltonian\n");
	Heff = C * diag(e) * C'; 
	disp(Heff);

	#
	#	Diagonalize Heff, just to recover original energies
	#
	Heff = .5 * (Heff + Heff');
	[Ceff,Eeff] = eig(Heff);
	e_eff = diag(Eeff);
	[uselessVariable,permutation]=sort(e_eff,'ascend');
	e_eff=e_eff(permutation);
	Ceff=Ceff(:,permutation);
		
	printf("\n");
	printf("State Energies:\n");
	printf("   %10s: %20s %20s (%s)\n","State","Ab Initio", "Effective Ham","units")	
	for i=1:nstates
		printf("   %10i: %20.8f %20.8f (%s)\n",i,e(i),e_eff(i),"au")	
	endfor
	printf("\n");


	printf("\n");
	printf(" Exchange coupling constants:\n");	
	printf("   	H = J * Sa * Sb \n");	
	printf("   		or\n");	
	printf("   	H = -2J' * Sa * Sb \n");	
	for i = 1:nstates
		for j = i+1:nstates
			Si = blocks(i)/2;	# total spin on site i
			Sj = blocks(j)/2;	# total spin on site j
			J = Heff(i,j)/sqrt(Si*Sj);
			printf(" J(%2i,%2i) = %20.6f  : J'(%2i,%2i) = %20.6f  (%s)\n", i,j,J*convert, i,j, J*(-.5)*convert, units);
		endfor
	endfor
	
	printf("\n");
	printf(" Norms\n");
	
	c = c*Xinv; #Unorthogonalize
	for i = 1:nstates 
		printf("  Vector %4i: %8.4f\n",i,norm(c(:,i)));
	endfor

