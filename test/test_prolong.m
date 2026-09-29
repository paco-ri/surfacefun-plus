% Test prolong_leaves, which interpolates onto refine_leaves meshes.

fprintf('=== test_prolong start ===\n');

% geometry parameters
n = 6; nu = 6; nv = 3;
dom0 = surfacemesh.sharptorus(n, nu, nv);
npat0 = length(dom0.x); % number of base patches
p2q0 = [(1:npat0).', zeros(npat0,2)]; % uniform base mesh
rmax = 4; % max quadtree depth

rng(0);
% smooth test function sampled on the base mesh
f = @(x,y,z) 1 + 2*x - 3*y.^2 + x.*y.*z + z.^3;
u0 = surfacefun(f, dom0);

% Match patches between two meshes by centroid.
function idx = matchpatches(da, db)
    ca = [cellfun(@(a) mean(a(:)), da.x), cellfun(@(a) mean(a(:)), da.y), ...
          cellfun(@(a) mean(a(:)), da.z)];
    cb = [cellfun(@(a) mean(a(:)), db.x), cellfun(@(a) mean(a(:)), db.y), ...
          cellfun(@(a) mean(a(:)), db.z)];
    idx = zeros(size(ca,1),1);
    for k = 1:size(ca,1)
        [~, idx(k)] = min(sum((cb - ca(k,:)).^2, 2));
    end
end

%% A. mark everything -> must equal surfacefun/refine(u0,1) exactly
[dA, ~, qA] = surfacemesh.refine_leaves(dom0, p2q0, 1:npat0, rmax);
uA = prolong_leaves(u0, p2q0, dA, qA);
uref = refine(u0, 1);
dref = uref.domain;
% the two meshes order patches differently
idx = matchpatches(dA, dref);
eA = 0; enode = 0;
for k = 1:length(dA.x)
    eA = max(eA, max(abs(uA.vals{k} - uref.vals{idx(k)}),[],'all'));
    enode = max(enode, max(abs(dA.x{k} - dref.x{idx(k)}),[],'all'));
end
fprintf('TESTA vs surfacefun/refine: node err %.3e, value err %.3e\n', enode, eA);

%% B. one direct two-level prolong == two chained one-level prolongs
% refine half the base patches, then a third of the resulting leaves
marked1 = randperm(npat0, max(1,floor(npat0/2)));
[d1, ~, q1] = surfacemesh.refine_leaves(dom0, p2q0, marked1, rmax);
marked2 = randperm(size(q1,1), max(1,floor(size(q1,1)/3)));
[d2, ~, q2] = surfacemesh.refine_leaves(dom0, q1, marked2, rmax);

u_chain  = prolong_leaves(prolong_leaves(u0, p2q0, d1, q1), q1, d2, q2);
u_direct = prolong_leaves(u0, p2q0, d2, q2);
eB = 0;
for k = 1:length(d2.x)
    eB = max(eB, max(abs(u_chain.vals{k} - u_direct.vals{k}),[],'all'));
end
fprintf('TESTB chained == direct: %d patches, err %.3e\n', length(d2.x), eB);
fprintf('TESTB levels: %s\n', mat2str(unique(q2(:,2)).'));

%% C. exactness on data that really is polynomial in reference coords
%   Build vals as a degree-(n-1) polynomial in the parent (u,v) directly.
xc = chebpts(n,[-1 1]);
[UU, VV] = meshgrid(xc, xc);          % VV = rows, UU = cols
pol = @(U,V) 3 - U + 2*V.^2 + U.^3.*V - 0.5*V.^4;
vals0 = cell(npat0,1);
for k = 1:npat0, vals0{k} = pol(UU, VV); end
w0 = surfacefun(vals0, dom0);
[dC, ~, qC] = surfacemesh.refine_leaves(dom0, p2q0, marked1, rmax);
wC = prolong_leaves(w0, p2q0, dC, qC);
eC = 0;
for k = 1:length(dC.x)
    l = qC(k,2); m = uint64(qC(k,3));
    % exact values: pol evaluated at this leaf's nodes in parent coordinates
    if l == 0
        want = pol(UU, VV);
    else
        [mx,my] = quadforest.deinterleave(m,l);
        h = 2/2^l;
        su = chebpts(n,[-1+h*double(mx), -1+h*(double(mx)+1)]);
        sv = chebpts(n,[-1+h*double(my), -1+h*(double(my)+1)]);
        [SU,SV] = meshgrid(su,sv);
        want = pol(SU,SV);
    end
    eC = max(eC, max(abs(wC.vals{k} - want),[],'all'));
end
fprintf('TESTC polynomial exactness: max err %.3e\n', eC);

fprintf('DONE_SENTINEL\n');
