    
function x = upper_solve(U, b)
    [n,n]=size(U);
    x(n)= b(n)/U(n,n)
    for i= (n-1):-1:1
           som=0;
        for j= (i+1):n
            som=som + U(i,j)*x(j)
        end
        x(i)=(1/U(i,i))*(b(i)-som)
       end
endfunction


function x = lower_solve(L, b)
    [m,n]=size(L);
    x(1)= b(1)/L(1,1)
    for i= 1 : m
        som=0;
        for j= 1: (i-1) 
            som=som + L(i,j)*x(j)
        end
        x(i)=(1/L(i,i))*(b(i)-som)
    end
endfunction


function [U,bb]=gauss_pivot(A,b)
     [n] =length(b)
    [m1,m2]=size(A);
    bb=b;
    U=A;
    for i = 1 : n-1 
    num=i
     for k = 1 : n-1
         if U(k,i)>U(num,i)  then
             num=k
         end
         if num<>i then
             for j = i : n 
                 tmp=A(num,j)
                 A(num,j)=A(i,j)
                 A(i,j)=tmp
             end
         end
     end
     pivot = U(i,i)
     for k = (i+1):n
        fact=(U(k,i)/pivot)
      for j = i : n
          U(k,j)=U(k,j)-fact*U(i,j)
     end
     bb(k)=bb(k)-fact*bb(i)
 end
 X(n)=(bb(n)/U(n,n))

 end
endfunction






    
 
 
         
         
