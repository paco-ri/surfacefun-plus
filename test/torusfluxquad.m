function [qnodes, qweights] = torusfluxquad(nr, nt)
%TORUSFLUXQUAD Computes quadrature for cross-section of surfacemesh.torus
%   Gauss-Legendre in r, periodic trapezoidal in theta
% 
%   Arguments:
%     nr: [int] number of r nodes
%     nt: [int] number of theta nodes
%   Returns:
%     qnodes: [double(3,nr*nt)] quadrature nodes
%     qweights: [double(1,nr*nt)] quadrature weights

[rnodes, rweights] = chebpts(nr,[0 1],1);
qnodes = zeros([3 nr*nt]);
qweights = zeros([1 nr*nt]);
for i = 1:nr
    rr = rnodes(i);
    wr = rweights(i);
    for j = 1:nt
        ij = (i-1)*nt+j;
        tt = 2*pi*(j-1)/nt;
        [g1, ~, g2] = evalTorus(tt,0);
        [dg1, ~, dg2] = duEvalTorus(tt,0);
        g1 = g1 - 4.5; 
        qnodes(:,ij) = rr*[g1; 0; g2] + [4.5; 0; 0];
        qweights(1,ij) = (2*pi/nt) * wr*rr*(g1*dg2-g2*dg1);
    end
end

end

function [x, y, z] = evalTorus(u, v)

d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.45 0];

x = zeros(size(u));
y = zeros(size(u));
z = zeros(size(u));

for i = -1:2
    for j = -1:2
        ii = i+2;
        jj = j+2;
        x = x + d(ii,jj)*cos(v).*cos((1-i)*u+j*v);
        y = y + d(ii,jj)*sin(v).*cos((1-i)*u+j*v);
        z = z + d(ii,jj)*sin((1-i)*u+j*v);
    end
end

% x = cos(v).*(4.5 + 2*cos(u));
% y = sin(v).*(4.5 + 2*cos(u));
% z = 3*sin(u);
 
end

function [x, y, z] = duEvalTorus(u, v)

d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.45 0];

x = zeros(size(u));
y = zeros(size(u));
z = zeros(size(u));

for i = -1:2
    for j = -1:2
        ii = i+2;
        jj = j+2;
        x = x - (1-i)*d(ii,jj)*cos(v).*sin((1-i)*u+j*v);
        y = y - (1-i)*d(ii,jj)*sin(v).*sin((1-i)*u+j*v);
        z = z + (1-i)*d(ii,jj)*cos((1-i)*u+j*v);
    end
end

% x = -2*cos(v).*sin(u);
% y = -2*sin(v).*sin(u);
% z = 3*cos(u);
 
end