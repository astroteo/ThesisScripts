function [a,em,i,OM,om,theta]=RV2OrbPar_My(r,v,mu)

v=[v(1);v(2);v(3)];
r=[r(1);r(2);r(3)];
rm=norm(r);
vm=norm(v);
a=mu/(((2*mu)/rm)-(vm^2));
h=cross(r,v);
hm=norm(h);
l=cross(v,h);
e=l/mu-r/rm;
em=norm(e);
i=acos(h(3)/hm);
kv=[0;0;1];
jv=[0;1;0];
iv=[1;0;0];

if norm(cross(kv,h))==0
n=iv;  
else
n=cross(kv,h)/norm(cross(kv,h));
end

    if dot(n,jv)>0
        OM=acos(dot(n,iv));
    else
        OM=2*pi-acos(dot(n,iv));
    end

if em==0
    
   om=0;
else
    
    if dot(e,kv)>0
        om=acos(dot(n,e)/em);
    else
        om=2*pi-acos(dot(n,e)/em);
    end
end

if em==0;

    theta=acos(dot(r/rm,iv));
    
else
    
    if dot(r,v)>0
        theta=acos(dot(r,e)/(rm*em));
    else
        theta=2*pi-acos(dot(r,e)/(rm*em));
    end
end

end