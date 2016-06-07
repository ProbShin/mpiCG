A = mmread('../../Emilia/Emilia_923.mtx');
ix = load('../../Emilia/Emilia_923.metis.iperm'); 
p=16;
[m,~] = size(A);
U=ones(m,1);
ix = ix + U;
for i=1:1:p
  U(m-16+i,1)=0;
end
U(ix)=U;
b=A*U;
mmwrite('b_Emilia_order.mtx', b);



