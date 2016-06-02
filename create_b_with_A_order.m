A = mmread('../../cfd2/cfd2_E16/cfd2_aug16.mtx');
ix = load('../../cfd2/cfd2_E16/cfd2_aug16.metis.iperm'); 
p=16;
[m,~] = size(A);
U=ones(m,1);
ix = ix + U;
for i=1:1:p
  U(m-16+i,1)=0;
end
U(ix)=U;
b=A*U;
mmwrite('b_cfd2_aug16_order.mtx', b);



