function [u1,u2,u3]=computeDisplacementDippingShearZoneSurfaceTanhSinh( ...
    x1,x2,q1,q2,q3,L,T,W,theta,phi, ...
    epsv11p,epsv12p,epsv13p,epsv22pp,epsv23p,epsv33p,G,nu)
% function COMPUTEDISPLACEMENTDIPPINGSHEARZONESURFACETANHSINH computes the
% displacement field associated with deforming vertical shear zones
% considering the following geometry.
%
%                N (x1)
%   observation /
%     point    /
%    x1,x2 -> +----------- E (x2)
%
%                        N (x1)
%                       /
%                      /\
%                     /  \ strike     .-------+
%         source     /    \  .--------        |    E (x2)
%        q1,q2,q3 ->@-------------------------+----   
%                   |                         |   + s
%                w  :                         |  / s
%                i  |                         | / e
%                d  :                         |/ n
%                t  |                 .-------+ k
%                h  :        .--------       / c
%                   +---------              / i
%                   :       l e n g t h    / h
%                   |                     + t
%                   :                      
%                   |                 
%                   D (x3)
%
% Input:
% x1, x2             north and east coordinates of the observation point,
% q1, q2, q3         north, east and depth coordinates of the shear zone,
% L, T, W            length, thickness and width of the shear zone,
% theta (degree)     strike angle from north (from x1) of the shear zone.
% phi (degree)       dip angle from horizontal of the shear zone,
% epsvij             source strain component ij in the shear zone
%                    in the system of reference tied to the shear zone,
% G, nu              shear modulus and Poisson's ratio in the half space.
%
% Output:
% u1                 displacement component in the north direction,
% u2                 displacement component in the east direction,
% u3                 displacement component in the down direction.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - Feb 22, 2016, Singapore.
%
% SEE ALSO: unicycle

% Lame parameter
lambda=G*2*nu/(1-2*nu);

% array size
s=size(x1);

% rotate observation points to the shear-zone-centric system of coordinates
t1= (x1-q1)*cosd(theta)+(x2-q2)*sind(theta);
x2=-(x1-q1)*sind(theta)+(x2-q2)*cosd(theta);
x1=t1;

% convert strain tensor to reference coordinate system with
% R=[[sind(phi), cosd(phi)];[-cosd(phi), sind(phi)]];
epsv11=epsv11p;
epsv12= sind(phi)*epsv12p+cosd(phi)*epsv13p;
epsv13=-cosd(phi)*epsv12p+sind(phi)*epsv13p;
epsv22=+sind(phi)*(+epsv22p*sind(phi)+epsv23p*cosd(phi))+cosd(phi)*(+epsv23p*sind(phi)+epsv33p*cosd(phi));
epsv23=+sind(phi)*(-epsv22p*cosd(phi)+epsv23p*sind(phi))+cosd(phi)*(-epsv23p*cosd(phi)+epsv33p*sind(phi));
epsv33=-cosd(phi)*(-epsv22p*cosd(phi)+epsv23p*sind(phi))+sind(phi)*(-epsv23p*cosd(phi)+epsv33p*sind(phi));

% isotropic strain
epsvkk=epsv11+epsv22+epsv33;

% Green's functions
r=@(y1,y2,y3) sqrt((x1-y1).^2+(x2-y2).^2+(y3).^2);

y2=@(y2p,y3p) +y2p*sind(phi)+y3p*cosd(phi);
y3=@(y2p,y3p) -y2p*cosd(phi)+y3p*sind(phi)+q3;

G11=@(y1,y2,y3) 1/(16*pi*(1-nu))*( ...
    (3-4*nu)./r(y1,y2,y3)+1./r(y1,y2,y3)+(x1-y1).^2./r(y1,y2,y3).^3 ...
    +(3-4*nu)*(x1-y1).^2./r(y1,y2,y3).^3 ...
    +4*(1-2*nu)*(1-nu)*(r(y1,y2,y3).^2-(x1-y1).^2+r(y1,y2,y3).*y3)./(r(y1,y2,y3).*(r(y1,y2,y3)+y3).^2)...
    );
G12=@(y1,y2,y3) (x1-y1).*(x2-y2)/(16*pi*G*(1-nu)).*( ...
    1./r(y1,y2,y3).^3+(3-4*nu)./r(y1,y2,y3).^3 ...
    -4*(1-2*nu)*(1-nu)./(r(y1,y2,y3).*(r(y1,y2,y3)+y3).^2) ...
    );
G13=@(y1,y2,y3) (x1-y1)/(16*pi*G*(1-nu)).*( ...
    (-y3)./r(y1,y2,y3).^3+(3-4*nu)*(-y3)./r(y1,y2,y3).^3 ...
    +4*(1-2*nu)*(1-nu)./(r(y1,y2,y3).*(r(y1,y2,y3)+y3)) ...
    );
G21=@(y1,y2,y3) (x1-y1).*(x2-y2)/(16*pi*G*(1-nu)).*( ...
    1./r(y1,y2,y3).^3+(3-4*nu)./r(y1,y2,y3).^3 ...
    -4*(1-2*nu)*(1-nu)./(r(y1,y2,y3).*(r(y1,y2,y3)+y3).^2) ...
    );
G22=@(y1,y2,y3) 1/(16*pi*(1-nu))*( ...
    (3-4*nu)./r(y1,y2,y3)+1./r(y1,y2,y3)+(x2-y2).^2./r(y1,y2,y3).^3 ...
    +(3-4*nu)*(x2-y2).^2./r(y1,y2,y3).^3 ...
    +4*(1-2*nu)*(1-nu)*(r(y1,y2,y3).^2-(x2-y2).^2+r(y1,y2,y3).*(y3))./(r(y1,y2,y3).*(r(y1,y2,y3)+y3).^2)...
    );
G23=@(y1,y2,y3) (x2-y2)/(16*pi*G*(1-nu)).*( ...
    (-y3)./r(y1,y2,y3).^3+(3-4*nu)*(-y3)./r(y1,y2,y3).^3 ...
    +4*(1-2*nu)*(1-nu)./(r(y1,y2,y3).*(r(y1,y2,y3)+y3)) ...
    );
G31=@(y1,y2,y3) (x1-y1)/(16*pi*G*(1-nu)).*( ...
    +(4-4*nu)*(-y3)./r(y1,y2,y3).^3 ...
    -4*(1-2*nu)*(1-nu)./(r(y1,y2,y3).*(r(y1,y2,y3)+y3)) ...
    );
G32=@(y1,y2,y3) (x2-y2)/(16*pi*G*(1-nu)).*( ...
    (-y3)./r(y1,y2,y3).^3+(3-4*nu)*(-y3)./r(y1,y2,y3).^3 ...
    -4*(1-2*nu)*(1-nu)./(r(y1,y2,y3).*(r(y1,y2,y3)+y3)) ...
    );
G33=@(y1,y2,y3) 1/(16*pi*G*(1-nu))*( ...
    (3-4*nu)./r(y1,y2,y3)+(5-12*nu+8*nu^2)./r(y1,y2,y3)+(-y3).^2./r(y1,y2,y3).^3 ...
    +((3-4*nu)*(+y3).^2)./r(y1,y2,y3).^3 ...
    );


% numerical solution with tanh/sinh quadrature
h=0.01;
n=fix(1/h*2);
u1=zeros(size(x1));
u2=zeros(size(x1));
u3=zeros(size(x1));
for k=-n:n
    wk=(0.5*h*pi*cosh(k*h))./(cosh(0.5*pi*sinh(k*h))).^2;
    xk=tanh(0.5*pi*sinh(k*h));
    for j=-n:n
        wj=(0.5*h*pi*cosh(j*h))./(cosh(0.5*pi*sinh(j*h))).^2;
        xj=tanh(0.5*pi*sinh(j*h));
        u1=u1+wk*wj*IU1(xk,xj);
        u2=u2+wk*wj*IU2(xk,xj);
        u3=u3+wk*wj*IU3(xk,xj);
    end
end

% rotate displacement field to reference system of coordinates
t1=u1*cosd(theta)-u2*sind(theta);
u2=u1*sind(theta)+u2*cosd(theta);
u1=t1;


    function u = IU1(x,y)
        % function IU1 is the integrand for displacement component u1
        u=zeros(s);
        if epsv11 ~= 0 || lambda*epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv11)*T*W/4*(G11(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G11(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2)));
        end
        if epsv12 ~= 0
            u=u+           2*G*epsv12*T*W/4*(G21(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))    -G21(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
                +sind(phi)*2*G*epsv12*L*W/4*(G11((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G11((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
                +cosd(phi)*2*G*epsv12*L*T/4*(G11((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G11((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv13 ~= 0
            u=u+           2*G*epsv13*T*W/4*(G31(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))    -G31(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
                +sind(phi)*2*G*epsv13*L*T/4*(G11((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G11((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
                -cosd(phi)*2*G*epsv13*L*W/4*(G11((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G11((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2)));
        end
        if epsv22 ~= 0 || lambda*epsvkk ~= 0
            u=u+sind(phi)*(lambda*epsvkk+2*G*epsv22)*L*W/4*(G21((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G21((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               +cosd(phi)*(lambda*epsvkk+2*G*epsv22)*L*T/4*(G21((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G21((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv23 ~= 0
            u=u+sind(phi)*2*G*epsv23*L*T/4*(G21((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G21((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
               +sind(phi)*2*G*epsv23*L*W/4*(G31((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G31((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               -cosd(phi)*2*G*epsv23*L*W/4*(G21((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G21((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               +cosd(phi)*2*G*epsv23*L*T/4*(G31((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G31((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv33 ~= 0 || lambda*epsvkk ~= 0
            u=u+sind(phi)*(lambda*epsvkk+2*G*epsv33)*L*T/4*(G31((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G31((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
               -cosd(phi)*(lambda*epsvkk+2*G*epsv33)*L*W/4*(G31((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G31((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2)));
        end
    end

    function u = IU2(x,y)
        % function IU2 is the integrand for displacement component u2
        u=zeros(s);
        if epsv11 ~= 0 || lambda*epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv11)*T*W/4*(G12(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G12(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2)));
        end
        if epsv12 ~= 0
            u=u+           2*G*epsv12*T*W/4*(G22(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))    -G22(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
                +sind(phi)*2*G*epsv12*L*W/4*(G12((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G12((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
                +cosd(phi)*2*G*epsv12*L*T/4*(G12((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G12((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv13 ~= 0
            u=u+           2*G*epsv13*T*W/4*(G32(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))    -G32(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
                +cosd(phi)*2*G*epsv13*L*T/4*(G12((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G12((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
                -sind(phi)*2*G*epsv13*L*W/4*(G12((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G12((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2)));
        end
        if epsv22 ~= 0 || lambda*epsvkk ~= 0
            u=u+sind(phi)*(lambda*epsvkk+2*G*epsv22)*L*W/4*(G22((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G22((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               +cosd(phi)*(lambda*epsvkk+2*G*epsv22)*L*T/4*(G22((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G22((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv23 ~= 0
            u=u+sind(phi)*2*G*epsv23*L*T/4*(G22((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G22((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
               +sind(phi)*2*G*epsv23*L*W/4*(G32((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G32((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               -cosd(phi)*2*G*epsv23*L*W/4*(G22((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G22((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               +cosd(phi)*2*G*epsv23*L*T/4*(G32((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G32((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv33 ~= 0 || lambda*epsvkk ~= 0
            u=u+sind(phi)*(lambda*epsvkk+2*G*epsv33)*L*T/4*(G32((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G32((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
               -cosd(phi)*(lambda*epsvkk+2*G*epsv33)*L*W/4*(G32((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G32((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2)));
        end
    end

    function u = IU3(x,y)
        % function IU3 is the integrand for displacement component u3
        u=zeros(s);
        if epsv11 ~= 0 || lambda*epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv11)*T*W/4*(G13(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G13(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2)));
        end
        if epsv12 ~= 0
            u=u+           2*G*epsv12*T*W/4*(G23(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))    -G23(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
                +sind(phi)*2*G*epsv12*L*W/4*(G13((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G13((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
                +cosd(phi)*2*G*epsv12*L*T/4*(G13((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G13((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv13 ~= 0
            u=u+           2*G*epsv13*T*W/4*(G33(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))    -G33(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
                +cosd(phi)*2*G*epsv13*L*T/4*(G13((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G13((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
                -sind(phi)*2*G*epsv13*L*W/4*(G13((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G13((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2)));
        end
        if epsv22 ~= 0 || lambda*epsvkk ~= 0
            u=u+sind(phi)*(lambda*epsvkk+2*G*epsv22)*L*W/4*(G23((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G23((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               +cosd(phi)*(lambda*epsvkk+2*G*epsv22)*L*T/4*(G23((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G23((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv23 ~= 0
            u=u+sind(phi)*2*G*epsv23*L*T/4*(G23((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G23((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
               +sind(phi)*2*G*epsv23*L*W/4*(G33((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G33((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               -cosd(phi)*2*G*epsv23*L*W/4*(G23((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G23((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
               +cosd(phi)*2*G*epsv23*L*T/4*(G33((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G33((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0)));
        end
        if epsv33 ~= 0 || lambda*epsvkk ~= 0
            u=u+sind(phi)*(lambda*epsvkk+2*G*epsv33)*L*T/4*(G33((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))            -G33((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
               -cosd(phi)*(lambda*epsvkk+2*G*epsv33)*L*W/4*(G33((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G33((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2)));
        end
    end

end



% U1=@(x,y) ...
%    (lambda*epsvkk+2*G*epsv11)*T*W/4*(G11(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G11(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%                   +2*G*epsv12*T*W/4*(G21(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G21(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%                   +2*G*epsv13*T*W/4*(G31(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G31(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%    +sind(phi)*( ...
%                          2*G*epsv12*L*W/4*(G11((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G11((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv13*L*T/4*(G11((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G11((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%         +(lambda*epsvkk+2*G*epsv22)*L*W/4*(G21((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G21((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv23*L*T/4*(G21((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G21((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         +2*G*epsv23*L*W/4*(G31((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G31((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%         +(lambda*epsvkk+2*G*epsv33)*L*T/4*(G31((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G31((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ) ...
%    +cosd(phi)*( ...
%                          2*G*epsv12*L*T/4*(G11((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G11((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         -2*G*epsv13*L*W/4*(G11((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G11((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%         +(lambda*epsvkk+2*G*epsv22)*L*W/4*(G21((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G21((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         -2*G*epsv23*L*T/4*(G21((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G21((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv23*L*W/4*(G31((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G31((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%         -(lambda*epsvkk+2*G*epsv33)*L*T/4*(G31((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G31((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) );
%     
% U2=@(x,y) ...
%    (lambda*epsvkk+2*G*epsv11)*T*W/4*(G12(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G12(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%                   +2*G*epsv12*T*W/4*(G22(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G22(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%                   +2*G*epsv13*T*W/4*(G32(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G32(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%    +sind(phi)*( ...
%                          2*G*epsv12*L*W/4*(G12((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G12((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv13*L*T/4*(G12((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G12((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%         +(lambda*epsvkk+2*G*epsv22)*L*W/4*(G22((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G22((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv23*L*T/4*(G22((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G22((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         +2*G*epsv23*L*W/4*(G32((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G32((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%         +(lambda*epsvkk+2*G*epsv33)*L*T/4*(G32((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G32((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ) ...
%    +cosd(phi)*( ...
%                          2*G*epsv12*L*T/4*(G12((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G12((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         -2*G*epsv13*L*W/4*(G12((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G12((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%         +(lambda*epsvkk+2*G*epsv22)*L*W/4*(G22((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G22((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         -2*G*epsv23*L*T/4*(G22((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G22((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv23*L*W/4*(G32((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G32((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%         -(lambda*epsvkk+2*G*epsv33)*L*T/4*(G32((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G32((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) );
%     
% U3=@(x,y) ...
%    (lambda*epsvkk+2*G*epsv11)*T*W/4*(G13(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G13(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%                   +2*G*epsv12*T*W/4*(G23(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G23(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%                   +2*G*epsv13*T*W/4*(G33(L,y2(x*T/2,(1+y)*W/2),y3(x*T/2,(1+y)*W/2))-G33(0,y2(x*T/2,(1+y)*W/2),y3(x*L/2,(1+y)*W/2))) ...
%    +sind(phi)*( ...
%                          2*G*epsv12*L*W/4*(G13((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G13((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv13*L*T/4*(G13((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G13((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%         +(lambda*epsvkk+2*G*epsv22)*L*W/4*(G23((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G23((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv23*L*T/4*(G23((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G23((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         +2*G*epsv23*L*W/4*(G33((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G33((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%         +(lambda*epsvkk+2*G*epsv33)*L*T/4*(G33((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G33((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ) ...
%    +cosd(phi)*( ...
%                          2*G*epsv12*L*T/4*(G13((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G13((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         -2*G*epsv13*L*W/4*(G13((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G13((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%         +(lambda*epsvkk+2*G*epsv22)*L*W/4*(G23((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G23((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%                         -2*G*epsv23*L*T/4*(G23((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G23((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) ...
%                         +2*G*epsv23*L*W/4*(G33((1+x)*L/2,y2(y*T/2,W),y3(y*T/2,W))-G33((1+x)*L/2,y2(y*T/2,0),y3(y*T/2,0))) ...
%         -(lambda*epsvkk+2*G*epsv33)*L*T/4*(G33((1+x)*L/2,y2(T/2,(1+y)*W/2),y3(T/2,(1+y)*W/2))-G33((1+x)*L/2,y2(-T/2,(1+y)*W/2),y3(-T/2,(1+y)*W/2))) );

