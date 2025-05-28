function [u,ut,gradu,D2u]=test_cases(t,z)

%global ucase;

x=z(:,1);
y=z(:,2);
u=(cos(t)*(cos(pi*x).*cos(pi*y)))';
if (nargout>1)
    ut(:,1)=-sin(t)*(cos(pi*x).*cos(pi*y));
    if (nargout>2)
        gradu(:,1)=cos(t)*(-pi*sin(pi*x).*cos(pi*y));
        gradu(:,2)=cos(t)*(-pi*cos(pi*x).*sin(pi*y));
        if (nargout>3)
            D2u(1,1,:)=cos(t)*(-pi^2*cos(pi*x).*cos(pi*y));
            D2u(1,2,:)=cos(t)*(pi^2*sin(pi*x).*sin(pi*y));
            D2u(2,1,:)=D2u(1,2,:);
            D2u(2,2,:)=D2u(1,1,:);
        end
    end
end

end
