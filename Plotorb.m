function R=Plotorb(a,e,i,OM,om,theta0,t_span,Mu)


[N,~]=size(t_span);

R=zeros(N,3);

for j=1:N
    
t=t_span(j);
theta=time2theta(a,e,theta0,t,Mu);
[r,~]=ParOrb2RV(a,e,i,OM,om,theta,Mu);
R(j,:)=r';

end




end