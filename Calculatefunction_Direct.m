%% Calculate the metrics of the network(directed network)
function  [degree_all,Ave_degree00,Node_RATE00,Paths,Shortest_L00,Eff_I00,Diameter00,Giant_tu00]=Calculatefunction_Direct(A,port00,port)

num_N = size(A,1);
numall=size(port,1); 
A00=sparse(A);
degreein = sum(A)'; 
degreeout = sum(A,2); 
topoA=A;
topoA(topoA~=0)=1;
degreein_TOPO(1:num_N) = sum(topoA(1:num_N,:),1);
degreeout_TOPO(1:num_N) = sum(topoA(:,1:num_N),2);
degree_all = degreein+ degreeout;
Ave_degree00 = sum(degree_all,1)/numall; 
Number_alone = find(degree_all==0);
Node_RATE00 = size(Number_alone,1)/numall; 

[Shortest_L00,I,M01,M02,Paths,L5,NL2,D] = Shortest_path(A);
pathin = (sum(Paths)./(num_N-1))'; 
pathout = sum(Paths,2)./(num_N-1); 
Eff = 1./Paths;
idy= Eff==Inf; 
Eff(idy)=0; 
Eff_I00 = sum(sum(Eff))/(numall *(numall -1));
Diameter00 = max(max(Paths));
z2=components(A00);% Giant_tu
z3=tabulate(z2(:));
Giant_tu00=max(z3(:,2));

