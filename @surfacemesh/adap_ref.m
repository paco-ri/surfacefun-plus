function [dom, qf, p2q] = adap_ref(dom, amr_tol, rmax, mode, marked)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    dom
    amr_tol
    rmax
    mode
    marked
end

arguments (Output)
    dom
    qf
    p2q
end

% mode = 2 -> refine w.r.t. second fundamental form
n = length(dom.x{1}); % polynomial order

% number of patches
npat = length(dom.x);
npat0 = npat;

% record points and first fundamental form
x = dom.x;
y = dom.y;
z = dom.z;
E = dom.E;
F = dom.F;
G = dom.G;

% allocate mapping from patch number to Q index
% q. index: (tree #, tree level, Morton encoding)
p2q = nan([npat0*4^rmax 3]);
% initially, no quadtree on any patch, because all patches are unrefined
p2q(1:npat0, 1) = 1:npat0;
p2q(1:npat0, 2:3) = zeros([npat0 2]);

if mode == 2
    % store second fundamental form
    L = cell(size(x));
    M = cell(size(x));
    N = cell(size(x));
    for i = 1:npat0
        [L{i}, M{i}, N{i}] = patch_ff2(n, x{i}, y{i}, z{i});
    end
end

r = 1; % refinement level
dont_have_to_refine = 0;
while r <= rmax && ~dont_have_to_refine
    dont_have_to_refine = 1;
    marked_new = zeros(1,4*length(marked)); % marked patches at next iter
    % initialize to maximum possible length, will cull zeroes later
    i = 1; % index for marked_new
    
    for p = marked

        if mode == 2
            [x_fin,y_fin,z_fin,...
            E_fin,F_fin,G_fin,...
            L_fin,M_fin,N_fin,...
            E_err,F_err,G_err,...
            L_err,M_err,N_err] ...
            = refine_patch_ff2(...
            r,n,x{p},y{p},z{p},E{p},F{p},G{p},L{p},M{p},N{p});
        else
            [x_fin,y_fin,z_fin,E_fin,F_fin,G_fin,E_err,F_err,G_err] ...
            = refine_patch_ff1(r,n,x{p},y{p},z{p},E{p},F{p},G{p});
        end

        % mark fine patches for refinement if error is large
        add_fine_patches = 0;
        for j = 1:4

            % allowed error on patch is amr_tol
            % absolute error is enough (think about flat patch)
            if mode == 2
                ff_err = patchL2maxnorm({E_err{j} F_err{j} G_err{j} L_err{j} M_err{j} N_err{j}}, E_fin{j}.*G_fin{j} - F_fin{j}.^2);
            else
                ff_err = patchL2maxnorm({E_err{j} F_err{j} G_err{j}}, E_fin{j}.*G_fin{j} - F_fin{j}.^2);
            end

            if ff_err > amr_tol % refine
                if j == 1
                    marked_new(i) = p; % first small patch inherits number of parent
                else
                    marked_new(i) = npat + j - 1; % rest get new numbers
                end
                i = i + 1;
                dont_have_to_refine = 0;
                add_fine_patches = 1;
            end

        end

    if add_fine_patches
        % put new patches in x, y, z, E, F, G, L, M, N
        x{p} = x_fin{1};
        y{p} = y_fin{1};
        z{p} = z_fin{1};
        E{p} = E_fin{1};
        F{p} = F_fin{1};
        G{p} = G_fin{1};
        if mode == 2
            L{p} = L_fin{1};
            M{p} = M_fin{1};
            N{p} = N_fin{1};
        end
        p2q(p, 2) = p2q(p, 2) + 1; % patch p is one quadtree level higher
        p2q(p, 3) = 4 * p2q(p, 3); % update Morton encoding accordingly
        tree_label = p2q(p, 1);
        new_morton = p2q(p, 3);

        for k = 1:3
            x{npat+k} = x_fin{k+1};
            y{npat+k} = y_fin{k+1};
            z{npat+k} = z_fin{k+1};
            E{npat+k} = E_fin{k+1};
            F{npat+k} = F_fin{k+1};
            G{npat+k} = G_fin{k+1};
            if mode == 2
                L{npat+k} = L_fin{k+1};
                M{npat+k} = M_fin{k+1};
                N{npat+k} = N_fin{k+1};
            end
            p2q(npat+k, 1) = tree_label;
            p2q(npat+k, 2) = r; % quadtree level = refinement level
            p2q(npat+k, 3) = new_morton + k;
        end

        npat = npat + 3; % three new patches
    end    
    
    end

    marked = marked_new(marked_new~=0);
    r = r + 1;

end

% initialize quadforest
Q = cell(npat0, 1); 
for i = 1:npat0
    Q{i} = cell(rmax,1);
    for j = 1:rmax
        Q{i}{j} = uint64([]);
    end
end

% use p2q to populate Q
p2q = rmmissing(p2q); % remove NaNs
for i = 1:npat
    level = p2q(i, 2);
    if level > 0
        Q{p2q(i, 1)}{level} = [Q{p2q(i, 1)}{level} p2q(i, 3)];
    end
end

% for i in 1:npat0, if p2q(i, 2) is not zero, then there is a tree at Q{i}
tree_roots = find(p2q(:, 2) > 0);
qf = quadforest(Q, rmax, dom.connectivity.elem2elem, tree_roots);
% TODO get new patches required for balancing
% idea: at each level, collect list of patches to be refined and call
% refine routines below
% example: refine 1 --> 1 76 77 78; refine 76 --> 76 79 80 81
addl_p2q = nan([npat0*4^rmax 3]);
i = 1;
% columns of p2q are tree root index, quadtree level, Morton code
for t = 1:npat0
    for l = 1:r-1
        if ~isempty(qf.addl_patches{t}{l})
            p2q(t, 2) = p2q(t, 2) + 1; % patch is now one level deeper
            root = p2q(t, 1);
            [x_ref, y_ref, z_ref] = refine_patch(n, x{root}, y{root}, z{root});
            x{root} = x_ref{1}; y{root} = y_ref{1}; z{root} = z_ref{1};
            for p = qf.addl_patches{t}{l}
                addl_p2q(i, :) = [t l p];
                i = i + 1;
            end
        end
    end
end
addl_p2q = rmmissing(addl_p2q);
npat1 = size(p2q, 1);
p2q = unique([p2q; addl_p2q], 'rows', 'stable');
npat = size(p2q, 1);
% for each row [tree root, level, Morton], create a patch at tree root at
% level with Morton code using refine()
tree_root = nan; level = nan;
% x_par = x{1}; y_par = y{1}; z_par = z{1}; 
x_chi = x{1}; y_chi = y{1}; z_chi = z{1};
x_chi_sav = x_chi; y_chi_sav = y_chi; z_chi_sav = z_chi;

% computation of new patches generated during balancing
for i = (npat1 + 1):npat
    % fprintf('Processing patch %d of %d\n', i, npat)
    if p2q(i, 1) == tree_root
        % use previously computed parent
        if p2q(i, 2) ~= level
            level = p2q(i, 2);
            morton = p2q(i, 3);
            parent_morton = bitshift(morton, -2);
            x_par = x_chi_sav{mod(parent_morton,4)+1};
            y_par = y_chi_sav{mod(parent_morton,4)+1};
            z_par = z_chi_sav{mod(parent_morton,4)+1};
            [x_chi, y_chi, z_chi] = refine_patch(n, x_par, y_par, z_par);
        end
    else
        tree_root = p2q(i, 1);
        level = p2q(i, 2);
        morton = p2q(i, 3);
        % find new row of p2q that is [tree_root, level-1, parent_morton]
        x_par = dom.x{tree_root}; y_par = dom.y{tree_root}; z_par = dom.z{tree_root};
        for l = 1:(level - 1)
            [x_chi, y_chi, z_chi] = refine_patch(n, x_par, y_par, z_par);
            anc_morton = bitshift(morton, -2 * (r - level + l));
            x_par = x_chi{mod(anc_morton,4)+1}; 
            y_par = y_chi{mod(anc_morton,4)+1}; 
            z_par = z_chi{mod(anc_morton,4)+1};
        end
        [x_chi, y_chi, z_chi] = refine_patch(n, x_par, y_par, z_par);
        x_chi_sav = x_chi; y_chi_sav = y_chi; z_chi_sav = z_chi;
    end
    x{i} = x_chi{mod(p2q(i, 3),4)+1};
    y{i} = y_chi{mod(p2q(i, 3),4)+1};
    z{i} = z_chi{mod(p2q(i, 3),4)+1};
end

% plot(dom)
% hold on
% set(gcf, 'Color', 'w')
% ctr_ix = idivide(n, int32(2)) + 1;
% for i = 1:npat
%     text(1.1 * x{i}(ctr_ix, ctr_ix), 1.1 * y{i}(ctr_ix, ctr_ix), 1.1 * z{i}(ctr_ix, ctr_ix), int2str(i), 'Color', 'red');
% end

% get split info
split = qf.get_split(p2q);
% writematrix(dom.connectivity.elem2elem, "orig_conn.txt")
% writecell(split, "split.txt")

dom = surfacemesh(x, y, z, split);

end

function [x_fin,y_fin,z_fin] = refine_patch(n,x_patch,y_patch,z_patch)
%REFINE_PATCH refine a patch

% tools for refining a patch
x  = chebpts(n, [-1 1]);
xL = chebpts(n, [-1 0]);
xR = chebpts(n, [ 0 1]);
BL = barymat(xL, x);
BR = barymat(xR, x);

% get nodes of fine patches
x_fin = cell(4,1);
y_fin = cell(4,1);
z_fin = cell(4,1);
x_fin{1} = BL * x_patch * BL.';
y_fin{1} = BL * y_patch * BL.';
z_fin{1} = BL * z_patch * BL.';
x_fin{2} = BL * x_patch * BR.';
y_fin{2} = BL * y_patch * BR.';
z_fin{2} = BL * z_patch * BR.';
x_fin{3} = BR * x_patch * BL.';
y_fin{3} = BR * y_patch * BL.';
z_fin{3} = BR * z_patch * BL.';
x_fin{4} = BR * x_patch * BR.';
y_fin{4} = BR * y_patch * BR.';
z_fin{4} = BR * z_patch * BR.';

end

function I = patchL2maxnorm(fvals_list, J)
%PATCHL2MAXNORM   Largest patch L^2 norm over a cell array of n x n values.
I = 0;
for k = 1:numel(fvals_list)
    I = max(I, patchL2norm(fvals_list{k}, J));
end
end
