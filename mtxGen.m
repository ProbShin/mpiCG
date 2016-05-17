function [] = mtxgen(n)

for i = 1:n
    a(i)=i*i*0.1;
    %  b(i)=i*0.1;
    end



    A=diag(sparse(a));
    %B=diag(sparse(b));


    field ='real';
    precision = 8;


    fA=sprintf('A_%d.mtx',n)
    %comment = str2mat('matrix A');
    %[err] = mmwrite(fA, A, comment, field, precision);

    mmwrite(fA, A);

    %fB=sprintf('B_%d.mtx',n)
    %comment = str2mat('matrix B');
    %[err] = mmwrite(fB, B, comment, field, precision);



  end
