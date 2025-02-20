au2cm = 219474.63;

A = [
	1  -1  -1  -1  -1   1
     	-1   1  -1  -1   1  -1
       	-1  -1   1   1  -1  -1
	-1  -1  -1   1   1   1
  	-1   1   1  -1  -1   1
     	1  -1   1  -1   1  -1
       	1   1  -1   1  -1  -1
	1   1   1   1   1   1
];

b = [
	-9075.3767710
	-9075.3763888
	-9075.3767716
	-9075.3767339
	-9075.3767396
	-9075.3767339
	-9075.3767396
	-9075.3770166
];

A = -2*A;

#Make into a 2J model
do_2j = 0
if do_2j
	A = [A(:,2)+A(:,5),A(:,1)+A(:,3)+A(:,4)+A(:,6)]
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


