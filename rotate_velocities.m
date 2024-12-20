function [vs,vd] = rotate_velocities(vE,vN,rcv)

dv = [rcv.dv(:,1), rcv.dv(:,2)]./sqrt(rcv.dv(:,1).^2 + rcv.dv(:,2).^2);

vs = vE.*rcv.sv(:,1) + vN.*rcv.sv(:,2);
vd = vE.*dv(:,1) + vN.*dv(:,2);

end