A = mmread('../../Geo/Geo_1438.mtx');
[m,~] = size(A);
U=ones(m,1);

b=A*U;
mmwrite('b_Geo.mtx', b);



