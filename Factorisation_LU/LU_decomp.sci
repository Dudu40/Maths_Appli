
function [U,bb]=gauss_pivot(A,b)
     [n] =length(b)
    [m1,m2]=size(A);
    //on initialise les matrices U et bb aux matrices A et b en paramètre
    bb=b;
    U=A;
    for i = 1 : n-1 
    num=i
     for k = 1 : n-1
         if U(k,i)>U(num,i)  then
             num=k
         end
         //effectue un échange des lignes si la coordonnée num est infèrieur à la coordonnée k
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
         //on calcule le premier inconnue grace au pivot
        fact=(U(k,i)/pivot)
      for j = i : n
       //on determine les éléments de la matrice en remontant en se servant des valeurs U calculés aux lignes précèdentes
          U(k,j)=U(k,j)-fact*U(i,j)
     end
     //on détermine le vecteur bb en utilisant le fact obtenue dans le calcul ci dessus
     bb(k)=bb(k)-fact*bb(i)
 end

 end
endfunction

function x = upper_solve(U, b)
    [n,n]=size(U);
    x(n)= b(n)/U(n,n)
    // on récupère la valeur de l'inconnue en (n,n) de la matrice qui est déjà isolé
    for i= (n-1):-1:1
         // on effectue la somme des éléments d'une ligne fois les éléments du vecteur x déjà calculées
           som=0;
        for j= (i+1):n
            som=som + U(i,j)*x(j)
        end
       //on redescends la matrice en utilisant l'inconnue calculé à la ligne précèdente pour déterminer l'inconnue suivant
        x(i)=(1/U(i,i))*(b(i)-som)
       end
endfunction


function x = lower_solve(L, b)
    [m,n]=size(L);
    // on récupère la valeur de l'inconnue en (1,1) de la matrice qui est déjà isolé
    x(1)= b(1)/L(1,1)
    for i= 1 : m
        // on effectue la somme des éléments d'une ligne fois les éléments du vecteur x déjà calculées
        som=0;
        for j= 1: (i-1) 
            som=som + L(i,j)*x(j)
        end
       //on redescends la matrice en utilisant l'inconnue calculé à la ligne précèdente pour déterminer l'inconnue suivant
        x(i)=(1/L(i,i))*(b(i)-som)
    end
endfunction
function [LU]=my_LU(A)
    
    [n,n]=size(A);
    j=1;
    // on suppose que le pivot existe (pivot different de 0)
    bool=%t;
    // tant que le pivot est different de 0 et qu'on n'arrive pas au boout de la diagonale
    while ((j<=n-1) & (bool)) do
            // Si le pivot est nul
           if  (A(j,j)==0) then
               // la mtrice est = nan (pas de valeur)
               LU =%nan
               // on stope l'algorithme et on renvoie A=nan
               bool=%f;
              
           // si le pivot est non nul
       else
                //on effectue les opérations décrites dans la méthode générale 
                // Donc on divise les elements de la matrice par le pivot 
                //      on effectue des opérations sur les lignes sans les permuter
                for i=j+1:n
                    A(i,j)=A(i,j)/A(j,j)
                    for k=j+1:n
                        A(i,k)=A(i,k)-A(i,j)*A(j,k)
                    end
                 
                    // on renvoie la nouvelle matrice calculée après cette succesions d'opérations éléementaires
                    LU=A
                   end
            end 
        j=j+1;   
        end
endfunction


function [L,U]=extract_LU (LU)
    [n,n]=size(LU);
    
    // on initialise les  2 matrices a 0 
    L=zeros(n,n);
    U=zeros(n,n);
    
    // si la matrice LU a extraire n'existe pas alors les matrices L et U n'existent pas non plus
    if (LU==%nan) then
        L=%nan
        U=%nan
    // si LU existe 
    else
      
      // On parcour la partie supérieure de la matrice 
      
       for i=1:n
           for j=1:(i-1)
               // on remplit la partie supérieure de L par la partie supérieure de LU
               L(i,j)=LU(i,j);
           end
       end
       
       // on parcourt la diagonale
       
       for i=1: n
           // D'apres la definition, les elemnts de la diagonales sont égales a 1
           L(i,i)=1;
       end
       
       
       // on parcourt toute la matrice U (initialement nulle)
       for i=1:n
           for j=1:n
               // si l'elemnt de la matrice L est egal a 0
               if (L(i,j)==0) then
                   // alors cette partie inferieur de U prend les valeure de LU
                   U (i,j)=LU(i,j);
               end
           end
               
       end
       
       // on remplit la diagonale par les valeur de la diagonale de LU
        for i=1: n
           U(i,i)=LU(i,i);
       end
   end
           
endfunction


function [d]=my_det(A)
    [n,n]=size(A);
    // on initialise le determinant a 1
    d=1;
    
    // on calcul la decomposition LU de A
    A=my_LU(A);
    
    // On parcourt la diagonale de LU
    for i=1:n
        // on multiplie les valeurs de la diagonale
        d=d*A(i,i);
    end
    
endfunction


function x= my_linsolve(A,b)
    
    // on stocke la decomosition Lu de A dans A
    A=my_LU(A);
    
    //On extrait L et U de A 
    [L,U]=extract_LU(A);
    
    // on calcule y = L^(-1) x b
    y=lower_solve(L,b);
    // on calcule x = U^(-1) x y
    // renvoit x
    x=upper_solve(U,y);
    
    
    
endfunction






































