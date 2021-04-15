function [r] = pagerank (M)
   [m,n]=size(M);
   I=zeros(m,n);
   s=0;
   
   for i=1 : n
       I(i,i)=1;
   end
   
   M=M-I;
 
 r=kernel(M);
 
    for i=1 : n
      s=s+r(i);
   end
   
   for i=1 : n
       r(i)=r(i)/s
   end
   
endfunction

function [P] =markov (M,n,P0)
    [m,m]=size(M);
    P=P0;
    for i=1:n
        P=M*P;
         
    end
    
 P
  
endfunction

function [P]=markov_plot(M,n,P0)
 [m,m]=size(M);
    P=P0;
    for i=1:n
        P=M*P;
         
    end
   
 
 
endfunction

function G =google_matrix(M, p)
[m,n]=size(M);
  if exists("p", "l") == 0 then
       p = 0.85; 
    end
  if(0>p | p>1) then
           p=0.85;
       end
G=((1-p)/m)*ones(m,n)+p*M;
endfunction

function [r] = pagerank (M,p)
   [m,n]=size(M);
   I=zeros(m,n);
   s=0;
   
   for i=1 : n
       I(i,i)=1;
   end
    if exists("p", "l") == 0 then
       p = 1; 
    end
  if(0>p | p>1) then
           p=1;
       end
   M=google_matrix(M, p)-I;
 
 r=kernel(M);
 
    for i=1 : n
      s=s+r(i);
   end
   
   for i=1 : n
        r(i)=r(i)/s
   end
   
endfunction

