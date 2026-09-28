function [dom, qf, p2q, hist] = adap_ref(dom, amr_tol, rmax, mode, marked, opts)
%ADAP_REF   Refine a surface mesh until its fundamental forms are resolved.
%   [DOM, QF, P2Q] = ADAP_REF(DOM0, AMR_TOL, RMAX) repeatedly marks the leaves
%   of the base mesh DOM0 whose FF_INDICATOR exceeds AMR_TOL, and splits them
%   with REFINE_LEAVES, for up to RMAX passes. RMAX is also the maximum
%   quadtree depth. P2Q is the leaf list of DOM, as described in
%   REFINE_LEAVES, and QF the balanced quadforest.
%
%   ADAP_REF(DOM0, AMR_TOL, RMAX, MODE) selects the indicator: MODE = 1
%   (default) uses the first fundamental form, MODE = 2 both.
%
%   ADAP_REF(DOM0, AMR_TOL, RMAX, MODE, MARKED) refines only within the base
%   patches listed in MARKED (default: all). Leaves added elsewhere to keep
%   the 2:1 balance are not refined further.
%
%   ADAP_REF(..., marking='dorfler', theta=THETA) marks, on each pass, the
%   smallest set of leaves carrying a fraction THETA of the sum of ETA.^2,
%   in place of the AMR_TOL threshold. The default is marking='threshold'.
%
%   [DOM, QF, P2Q, HIST] = ADAP_REF(...) also returns a struct array with one
%   entry per mesh, starting from DOM0, with fields
%
%     p2q     the leaf list; REFINE_LEAVES(DOM0, P2Q, [], RMAX) rebuilds it
%     eta     FF_INDICATOR of each leaf
%     marked  the leaves marked for the next pass, including any already at
%             depth RMAX (these are not split)
%
%   See also SURFACEMESH/FF_INDICATOR, SURFACEMESH/REFINE_LEAVES.

arguments (Input)
    dom
    amr_tol (1,1) double
    rmax    (1,1) double {mustBeInteger, mustBePositive}
    mode    (1,1) double {mustBeMember(mode, [1 2])} = 1
    marked  = 1:length(dom.x)
    opts.marking {mustBeMember(opts.marking, {'threshold', 'dorfler'})} = 'threshold'
    opts.theta (1,1) double = 0.5
end

arguments (Output)
    dom
    qf
    p2q
    hist
end

% Every mesh is rebuilt from DOM0 and the leaf list, so start from one leaf
% per base patch.
dom0 = dom;
npat0 = length(dom0.x);
p2q = [(1:npat0).', zeros(npat0, 2)];
[dom, qf, p2q] = surfacemesh.refine_leaves(dom0, p2q, [], rmax);
inroot = false(npat0, 1);
inroot(marked) = true;

hist = struct('p2q', {}, 'eta', {}, 'marked', {});
% Up to RMAX refinement passes, plus a last one that only records ETA on
% the final mesh.
for pass = 0:rmax
    eta = surfacemesh.ff_indicator(dom, p2q, mode);
    % Zero ETA outside MARKED so those leaves are never marked and don't
    % count toward the Dorfler total.
    e = eta .* inroot(p2q(:, 1));
    switch opts.marking
        case 'threshold'
            mk = find(e > amr_tol);
        case 'dorfler'
            [se, ord] = sort(e, 'descend');
            c = cumsum(se.^2);
            mk = ord(1:find(c >= opts.theta*c(end), 1));
            % If every E is zero, FIND still returns one leaf.
            mk = mk(e(mk) > 0);
    end
    hist(end+1) = struct('p2q', p2q, 'eta', eta, 'marked', mk); %#ok<AGROW>

    % Leaves at depth RMAX stay marked in HIST but can't be split.
    mk = mk(p2q(mk, 2) < rmax);
    if ( pass == rmax || isempty(mk) )
        break
    end
    [dom, qf, p2q] = surfacemesh.refine_leaves(dom0, p2q, mk, rmax);
end

end
