	
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
	printf("Norms of projected eigenvectors (Projection onto neutral determinant basis)\n");
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

	
	#
	#	Diagonalize on-site blocks of effective hamiltonian:	
	#
	U = eye(dim_heis);
	U2 = eye(dim_heis);
	first = 1;
	for b = blocks
		last = first + b - 1;
		i = first;
		j = last;
		ht = Heff(i:j,i:j);
		ct = -[[1/sqrt(2),1/sqrt(2)];[1/sqrt(2),-1/sqrt(2)]];


		[ct,et] = eig(ht);
		
		[uselessVariable,permutation]=sort(diag(et),'ascend');
		et=et(permutation,permutation);
		ct=ct(:,permutation);

		#	Because the block diagonalization eigenvectors can
		#	introduce an arbitrary sign, ensure that the max 
		#	projection onto the unit vector (high-spin state)
		#	has the same sign.
		#	
		#		This should be done properly by computing
		#		S^2 or something for identifying the HS vector,
		#		since this could look different for systems
		#		where the orbital energies are very different - 
		#		although, this model is not very appropriate for 
		#		such systems.
		#
		hs = 1/sqrt(b) * ones(b,1);	# local high-spin vector
		hs = ct*hs;
		sign_hs = 1;
		if max(hs) > -min(hs)
			sign_hs = 1;
		elseif max(hs) < -min(hs)
			sign_hs = -1;
		else
			disp("Error finding sign_hs")
			return
		endif

		# 	Fill block of transformation matrix
		U(i:j,i:j) = sign_hs*ct;

		first = last + 1;
	endfor

	printf(" Transformation matrix which rotates the Heff into the local eigenstate basis\n");
	disp(U);
	Heff = U' * Heff * U;

	printf("\n");
	printf("Heff in block diagonal (local eigenstate) basis\n");
	disp(Heff);
	
	printf("\nFor this second projection to work, each block should have only 1 nonzero diagonal element.\n");
	printf("Diagonals of Heff:\n");
	for i = 1:dim_heis
		printf("  H(%i,%i) = %20.8f \n",i,i,Heff(i,i));
	endfor

	# Just ensure symmetry numerically (though we are already symmetric)
	printf("\n");
	Heff = .5 * (Heff + Heff');
	[Ceff,Eeff] = eig(Heff);

	[uselessVariable,permutation]=sort(diag(Eeff),'ascend');
	Eeff=Eeff(permutation,permutation);
	Ceff=Ceff(:,permutation);
	
	e_eff = diag(Eeff);

	printf("\n");
	printf("State Energies:\n");
	printf("   %10s: %20s %20s %20s\n","State","Ab Initio", "Effective Ham","Relative")	
	for i=1:nstates
		printf("   %10i: %20.8f %20.8f %20.4f\n",i,e(i),e_eff(i),(e(i)-e(1))*au2ev)	
	endfor
	
	#
	# Project to local groundstate basis
	P = zeros(dim_heis,nstates);
	Q = zeros(dim_heis,nstates);
	I = eye(dim_heis,dim_heis);

	q = 1;
	for i = 1:nstates
		Q(:,i) = I(:,i);
		P(:,i) = I(:,q);
		q += blocks(i);
	endfor
	
	c = P'*Ceff*Q;
	e = Q'*Eeff*Q;
	printf("\nProjected H Eigenvectors\n");
	disp(c);
	
	printf("\n");
	printf("Overlap of projected eigenvectors\n");
	S = c'*c;
	disp(S);
	[U,l] = eig(S);
	X 	= U * l^(-1/2) * U';
	Xinv 	= U * l^(1/2) * U';
	c = c*X;
	S = c'*c;
	
	printf("\nOrthogonalized Projected H Eigenvectors\n");
	disp(c);
	printf("\n");
	
	printf("\n");
	Heff = c*e*c';
	printf("\nFinal Effective Hamiltonian in local ground state basis:\n");
	disp(Heff);

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

