au2cm = 219474.63;


# Equations to solve
#    	|1	1	1	1	1 	1| |J_12|	|E(aaaaaa)|
#  	|1	1	1	1	-1	1| |J_23|    	|E(aaaaab)|
#  	|1	1	1	-1	-1	1| |J_34|    	|E(aaaaba)|
#  -9/2	|1	1	-1	-1	1 	1| |J_45|    	|E(aaabaa)|
#  	|1	-1	-1	1	1 	1| |J_56|    	|E(aabaaa)|
#  	|-1	-1	1	1	1 	1| | h  |   =	|E(abaaaa)|
#  	|-1 	1	1	1	1 	1|    		|E(baaaaa)|
#  	|-1	1	1	1	-1	1|    		|E(baaaab)|
#  	|1	-1	1	1	1 	1|    		|E(bbaaaa)|

#
#	BS-DFT energies: PBE0
b_pbe0 = [
-10067.6534164
-10067.6536486
-10067.6539145
-10067.6539154 
-10067.6538887
-10067.6538903
-10067.6536517
-10067.6538841
-10067.6536524
];

#
#	BS-DFT energies: 50/50 
b_5050 = [
-10069.7289824
-10069.7291304
-10069.7294699
-10069.7294288
-10069.7292008
-10069.7291975
-10069.7290876
-10069.7292356
-10069.7290917
];


#
#	BS-DFT energies: B3LYP 
b_b3lyp = [
-10073.1662314
-10073.1665462
-10073.1669043
-10073.1669143
-10073.1668799
-10073.1668820
-10073.1665590
-10073.1668738
-10073.1665518
];


A = [
1	1	1	1	1 
1	1	1	1	-1
1	1	1	-1	-1
1	1	-1	-1	1 
1	-1	-1	1	1 
-1	-1	1	1	1 
-1 	1	1	1	1 
-1	1	1	1	-1
1	-1	1	1	1 
];

A = -9/2*A;

b = b_b3lyp;
b = b_5050;
b = b_pbe0;

#Make into a 2J model
do_2j = 1
if do_2j
	A = [A(:,1)+A(:,5),A(:,2)+A(:,3)+A(:,4)]
end
	
#Add a column to allow one variable to describe all the background (non-spin) energy
A = [A,ones(rows(A),1)];

printf("A\n");
disp(A)

printf("b\n");
disp(b)
printf("\n");

# 
# Only the equations (rows of our linear system) listed in "list" will be used to compute J
list = 1:rows(b); 	#select all 
list = [1,2,3,4,5,6,7]; #select only HS and 1SF


printf(" Using these equations for solving for J's:\n");
disp(list);
printf("\n");
A = A(list,:);
b = b(list);

[U,s,V] = svd(A,'econ');
rank1 = length(diag(s));

[U,s,V] = svd([A,b],'econ');
rank2 = length(diag(s));

printf(" Rank of Coefficient/Augmented matrices: %i/%i\n",rank1,rank2);
printf("\n");

x = pinv(A)*b;
printf(" H = -2*JAB*SA*SB:\n");
printf("\n");
for s=1:length(x)
	letter = "J";
	if s == length(x)
		letter = "h";
	end
	printf(" %s%-2i = %16.8f au %16.3f cm-1 \n",letter,s,x(s),x(s)*au2cm);
end
err = A*x-b;
printf("\n Norm of Error for this solution: Error = Ax-b\n");
printf(" Error: 	%12.8f au	%12.3f cm-1\n",norm(err),au2cm*norm(err));


return;


