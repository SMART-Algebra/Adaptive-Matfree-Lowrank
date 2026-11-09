classdef MatFreeOp
% MatFreeOp  Matrix-free linear operator for the randQB_* routines.
%
%   A = MatFreeOp(m, n, Afun)          % forward products only (randQB_AF_*)
%   A = MatFreeOp(m, n, Afun, Atfun)   % forward and adjoint (all routines)
%
%   Afun(X)  must return A*X   for an n-by-p block X (p = the block size b),
%   Atfun(Y) must return A'*Y  for an m-by-p block Y (A' = conjugate transpose).
%
%   The routines use only the following operations on A, and the class
%   provides exactly these:
%       [m, n] = size(A),  size(A, dim)
%       A * X        ->  Afun(X)
%       X * A        ->  (A' * X')'  =  (Atfun(X'))'      (the routines call Qi' * A)
%       A'           ->  the adjoint operator (the two handles swapped)
%
%   Products are applied block by block, so a handle that accepts only one
%   column must loop over the columns of X (Test_MatFreeOp.m shows a wrapper).
%   Anything that needs the entries of A, such as norm(A) or A - Q*B, is not
%   available for an operator, so the routines' compute_exact option must
%   stay false.  Note that randQB_MF_Spec estimates ||A||_2 from its first
%   block, so a relative tolerance still works.
%
%   Example: an operator available only through sparse solves.
%       dK = decomposition(K);                              % factorize once
%       A  = MatFreeOp(n, n, @(X) dK \ (M*X), @(Y) M' * (dK' \ Y));
%       [Q, B]    = randQB_MF_Fro(A, 16, 1e-6);
%       [U, S, V] = QB_to_SVD(Q, B);
%   The handles may capture whatever they need (factorizations, meshes,
%   solver objects).  A PDE Jacobian J is passed the same way: Afun runs the
%   linearized forward solves and Atfun the adjoint solves.

    properties (SetAccess = immutable)
        m        % number of rows
        n        % number of columns
        Afun     % X -> A*X
        Atfun    % Y -> A'*Y, or [] when no adjoint is available
    end

    methods
        function A = MatFreeOp(m, n, Afun, Atfun)
            if nargin < 4
                Atfun = [];
            end
            A.m    = m;
            A.n    = n;
            A.Afun = Afun;
            A.Atfun = Atfun;
        end

        function varargout = size(A, dim)
            sz = [A.m, A.n];
            if nargin > 1
                if dim <= 2
                    varargout{1} = sz(dim);
                else
                    varargout{1} = 1;
                end
            elseif nargout <= 1
                varargout{1} = sz;
            else
                varargout = num2cell([sz, ones(1, nargout - 2)]);
            end
        end

        function C = mtimes(X, Y)
            if isa(X, 'MatFreeOp') && isa(Y, 'MatFreeOp')
                error('MatFreeOp:mtimes', ...
                      'The product of two MatFreeOp objects is not supported.');
            elseif isa(X, 'MatFreeOp')
                C = X.Afun(Y);              % A * Y
            else
                C = (Y' * X')';             % X * A = (A' * X')'
            end
        end

        function At = ctranspose(A)
            if isempty(A.Atfun)
                error('MatFreeOp:noAdjoint', ...
                      ['No adjoint handle was given; use randQB_AF_Fro or ', ...
                       'randQB_AF_Spec, which need only A*X.']);
            end
            At = MatFreeOp(A.n, A.m, A.Atfun, A.Afun);
        end
    end
end
