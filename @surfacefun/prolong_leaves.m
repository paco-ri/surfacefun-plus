function u = prolong_leaves(u, p2q_coarse, dom_fine, p2q_fine)
%PROLONG_LEAVES   Transfer a SURFACEFUN onto a refinement of its mesh.
%   V = PROLONG_LEAVES(U, P2Q_COARSE, DOM_FINE, P2Q_FINE), using polynomial 
%   interpolation, evaluates the SURFACEFUN U, defined on the mesh described 
%   by P2Q_COARSE, at the nodes of the finer mesh DOM_FINE described by 
%   P2Q_FINE, and returns the result as a SURFACEFUN on DOM_FINE.
%
%   Both leaf lists must refer to the same base (level-0) mesh, and every
%   leaf of P2Q_FINE must be a descendant of, or equal to, some leaf of
%   P2Q_COARSE -- which is what SURFACEMESH/REFINE_LEAVES guarantees.
%

arguments (Input)
    u
    p2q_coarse (:,3) double
    dom_fine
    p2q_fine   (:,3) double
end

vals = u.vals;
if ( size(p2q_coarse, 1) ~= numel(vals) )
    error('SURFACEFUN:prolong_leaves:size', ...
        'P2Q_COARSE has %d rows but U has %d patches.', ...
        size(p2q_coarse, 1), numel(vals));
end
if ( size(p2q_fine, 1) ~= length(dom_fine.x) )
    error('SURFACEFUN:prolong_leaves:size', ...
        'P2Q_FINE has %d rows but DOM_FINE has %d patches.', ...
        size(p2q_fine, 1), length(dom_fine.x));
end

n = size(vals{1}, 1);
xc = chebpts(n, [-1 1]);

% Hash the coarse leaves on (root, level, morton) for O(1) ancestor lookup.
key = @(t, l, m) sprintf('%d_%d_%d', t, l, m);
lookup = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:size(p2q_coarse, 1)
    lookup(key(p2q_coarse(i,1), p2q_coarse(i,2), p2q_coarse(i,3))) = i;
end

npat = size(p2q_fine, 1);
newvals = cell(npat, 1);
for j = 1:npat
    t = p2q_fine(j, 1);
    l = p2q_fine(j, 2);
    m = uint64(p2q_fine(j, 3));

    % Walk up from this leaf until we land on a coarse leaf.
    ic = 0;
    for lc = l:-1:0
        mc = bitshift(m, -2*(l - lc));
        k = key(t, lc, double(mc));
        if ( isKey(lookup, k) )
            ic = lookup(k);
            break
        end
    end
    if ( ic == 0 )
        error('SURFACEFUN:prolong_leaves:ancestor', ...
            ['Fine patch %d (root %d, level %d, morton %d) has no ancestor ' ...
             'among the coarse leaves.'], j, t, l, p2q_fine(j,3));
    end

    if ( lc == l )
        % Same cell on both meshes; nothing to interpolate.
        newvals{j} = vals{ic};
        continue
    end

    % Parameter square of the coarse leaf and of the fine leaf, both in the
    % [-1,1]^2 coordinates of their common level-0 ancestor.
    [cx, cy] = quadforest.deinterleave(bitshift(m, -2*(l - lc)), lc);
    [fx, fy] = quadforest.deinterleave(m, l);
    hc = 2/2^lc;
    hf = 2/2^l;
    au = -1 + hc*double(cx);  av = -1 + hc*double(cy);
    cu = -1 + hf*double(fx);  cv = -1 + hf*double(fy);

    % Map the fine square into the coarse leaf's own [-1,1]^2 coordinates.
    su = [-1 + 2*(cu - au)/hc, -1 + 2*(cu + hf - au)/hc];
    sv = [-1 + 2*(cv - av)/hc, -1 + 2*(cv + hf - av)/hc];

    Bu = barymat(chebpts(n, su), xc);
    Bv = barymat(chebpts(n, sv), xc);
    newvals{j} = Bv * vals{ic} * Bu.';
end

u = surfacefun(newvals, dom_fine);

end
