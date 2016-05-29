A = mmread('../../StocF/StocF.mtx');
[m,~] = size(A);
U=ones(m,1);

b=A*U;
mmwrite('b_StocF.mtx', b);



