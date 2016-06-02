A = mmread('../../Emilia/Emilia_923_order.mtx');
[m,~] = size(A);
U=ones(m,1);

b=A*U;
mmwrite('b_Emilia_order.mtx', b);



