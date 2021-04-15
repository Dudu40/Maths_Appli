
function y = euler(f,y0,t)  // fonction euler qui determine y selon t à partir d'une fonction, de la valeur initiale de y et du point intérieur t
    y=zeros(t) //initialise le tableau des valeurs de y
    y(1)=y0
    for i = 2:length(t) 
        y(i) = y(i-1) +(t(i)-t(i-1)) * f(t(i-1),y(i-1)) //calcule y de facon réccurente
    end
endfunction

function y=rk2(f,y0,t)
    y=zeros(t)  //initialise le tableau des valeurs de y
    y(1)=y0
    for i = 2:length(t)
        y(i)=y(i-1)+(t(i)-t(i-1))*f( t(i-1)+ (t(i)-t(i-1))/2   ,y(i-1)+ ((t(i)-t(i-1))/2 )*f(t(i-1),y(i-1))) //calcule y de facon reccurente
    end
endfunction

function ydot= fs(t,s,i,r) //determine le nombre de personnes saines
    ydot = -0.9*i*s;
endfunction


function ydot=fi(t,s,i,r) //determine le nombre de personnes inféctés
    ydot =0.9*i*s - i/10;
endfunction

function ydot=fr(t,s,i,r) //determine le nombre de personnes rétablies
    ydot = i/10;
endfunction


function [s,i,r] = epidemie (a,b,c,s0,i0,r0,t)  // determine l'evolution du nombre de personnes saines,rétablies et inféctés au cours du temps
      s=zeros(t);
      i=zeros(t);
      r=zeros(t);
      s(1)=s0
      i(1)=i0
      r(1)=r0

      for j=2:length(t)
        s(j)=s(j-1)+(t(j)-t(j-1))*a(t(j-1),s(j-1),i(j-1),r(j-1))
        i(j)=i(j-1)+(t(j)-t(j-1))*b(t(j-1),s(j-1),i(j-1),r(j-1))
        r(j)=r(j-1)+(t(j)-t(j-1))*c(t(j-1),s(j-1),i(j-1),r(j-1))
      end
      
      
 
endfunction






