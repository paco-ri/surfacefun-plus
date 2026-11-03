function [dom, qf, p2q, hist] = adap_ref(dom, amr_tol, rmax, mode, marked)
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
    mk = find(eta > amr_tol & inroot(p2q(:, 1)));
    hist(end+1) = struct('p2q', p2q, 'eta', eta, 'marked', mk); %#ok<AGROW>

    % Leaves at depth RMAX stay marked in HIST but can't be split.
    mk = mk(p2q(mk, 2) < rmax);
    if ( pass == rmax || isempty(mk) )
        break
    end
    [dom, qf, p2q] = surfacemesh.refine_leaves(dom0, p2q, mk, rmax);
end

end
