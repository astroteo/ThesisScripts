function [theta]=time2theta(a,e,theta0,dt,mu)

theta0=theta0*pi/180;
T=(2*pi*sqrt((a^3)/mu));
n=floor(dt/T);
m=ceil(dt/T);
if theta0<pi
    T_theta0=sqrt((a^3)/mu)*(2*atan(tan(theta0/2)*sqrt((1-e)/(1+e)))-e*sin(2*atan(tan(theta0/2)*sqrt((1-e)/(1+e)))));
else
    T_theta0=T+sqrt((a^3)/mu)*(2*atan(tan(theta0/2)*sqrt((1-e)/(1+e)))-e*sin(2*atan(tan(theta0/2)*sqrt((1-e)/(1+e)))));
end
if dt-n*T<T-T_theta0
    DT=T_theta0+dt-n*T;
    fE=@(x) x-e*sin(x)-DT*sqrt(mu/(a^3));
    E=fzero(fE,0);
    if DT<T/2
        theta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    else
        theta=2*pi+2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    end
else
    DT=T_theta0+dt-m*T;
    fE=@(x) x-e*sin(x)-DT*sqrt(mu/(a^3));
    E=fzero(fE,0);
    if DT<T/2
        theta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    else
        theta=2*pi+2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    end
end

end
