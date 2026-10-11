function [c, x, loop, h_norm] = topAC99(nelx, nely, volfrac, penal, w, beta)
%% 0. PHYSICAL INTERFACE PARAMETERS
ell = w / (sqrt(2)*log(9));
eps = ell * sqrt(beta);
%% 1. MATERIAL & FINITE ELEMENT SETUP
E0 = 1.0; Emin = 1e-9; nu = 0.3;
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12; A12' A11] + nu*[B11 B12; B12' B11]);
nodenrs = reshape(1:(nelx+1)*(nely+1), nely+1, nelx+1);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely, 1);
edofMat = repmat(edofVec,1,8) + repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1);
iK = reshape(kron(edofMat, ones(8,1))', 64*nelx*nely, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nelx*nely, 1);
U = zeros(2*(nely+1)*(nelx+1), 1);
F = sparse(2, 1, -1, 2*(nely+1)*(nelx+1), 1);
fixeddofs = [1:2:2*(nely+1), 2*(nelx+1)*(nely+1)];
freedofs = setdiff(1:2*(nely+1)*(nelx+1), fixeddofs);
%% 2. CONFORMING LAPLACIAN & MASS OPERATORS
Dx = spdiags([-ones(nelx+1,1) 2*ones(nelx+1,1) -ones(nelx+1,1)], -1:1, nelx+1, nelx+1);
Dx(1,1) = 1; Dx(end,end) = 1;
Dy = spdiags([-ones(nely+1,1) 2*ones(nely+1,1) -ones(nely+1,1)], -1:1, nely+1, nely+1);
Dy(1,1) = 1; Dy(end,end) = 1;
mx = [0.5; ones(nelx-1,1); 0.5]; my = [0.5; ones(nely-1,1); 0.5];
M = my * mx'; Minv = 1 ./ M; Omega = nelx * nely;
Klap = kron(Dx, spdiags(my, 0, nely+1, nely+1)) + kron(spdiags(mx, 0, nelx+1, nelx+1), Dy);
%% 3. INITIALIZATION
phi = repmat(volfrac, nely+1, nelx+1); p_old1 = phi; p_old2 = phi;
low = phi - 0.5; upp = phi + 0.5;
lambda_h = zeros(nely+1, nelx+1); rho_h = 0.1; c_old = 1e6; move = 0.10;
%% 4. AUGMENTED LAGRANGIAN MMA OPTIMIZATION LOOP
for loop = 1:500
    x = conv2(phi, [1 1; 1 1]/4, 'valid');
    % FE-Analysis & Sensitivity
    sK = reshape(KE(:)*(Emin + (E0 - Emin)*x(:)'.^penal), 64*nelx*nely, 1);
    K = sparse(iK, jK, sK); K = (K + K')/2;
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat), 2), nely, nelx);
    c = sum(sum((Emin + (E0 - Emin)*x.^penal).*ce));
    dc = -penal * (E0 - Emin) * x.^(penal-1) .* ce;
    grad_c = conv2(dc, [1 1; 1 1]/4, 'full');
    % Conserved Allen-Cahn Residual & Zero-Adjoint Gradient
    dW = beta * (2*phi - 6*phi.^2 + 4*phi.^3);
    d2W = beta * (2 - 12*phi + 12*phi.^2);
    mu = sum(M(:).*dW(:)) / Omega;
    h = eps^2 * reshape(Klap*phi(:), nely+1, nelx+1) + M.*(dW - mu);
    h_norm = sqrt(sum(h(:).^2 ./ M(:)) / Omega);
    z = lambda_h + rho_h * (h .* Minv);
    z_bar = sum(M(:).*z(:)) / Omega;
    grad_h = eps^2 * reshape(Klap*z(:), nely+1, nelx+1) + M.*d2W.*(z - z_bar);
    % Total Sensitivity & Moving Asymptotes Update
    w_h = min(1, 0.05 + 0.06 * loop);
    df0 = grad_c + w_h * grad_h;
    df1 = M / Omega;
    if loop > 30 && abs(c - c_old)/c < 1e-3 && h_norm <= 0.28*beta, move = max(0.005, move*0.95); end
    if loop > 2
        gam = ones(size(phi));
        gam((phi - p_old1).*(p_old1 - p_old2) > 0) = 1.2;
        gam((phi - p_old1).*(p_old1 - p_old2) < 0) = 0.7;
        low = max(phi - 10, min(phi - 0.01, phi - gam.*(p_old1 - low)));
        upp = min(phi + 10, max(phi + 0.01, phi + gam.*(upp - p_old1)));
    end
    p_old2 = p_old1; p_old1 = phi;
    % Closed-Form 1D-Dual MMA Subproblem Solver
    p0 = (max(0, df0) + 0.001*abs(df0) + 1e-6).*(upp - phi).^2;
    q0 = (max(0, -df0) + 0.001*abs(df0) + 1e-6).*(phi - low).^2;
    p1 = df1.*(upp - phi).^2;
    r1 = sum(df1(:).*phi(:)) - volfrac - sum(p1(:)./(upp(:) - phi(:)));
    a_lim = max(0, max(phi - move, low + 0.1*(phi - low)));
    b_lim = min(1, min(phi + move, upp - 0.1*(upp - phi)));
    phi_trial = max(a_lim, min(b_lim, (sqrt(p0).*low + sqrt(q0).*upp)./(sqrt(p0) + sqrt(q0))));
    if r1 + sum(p1(:)./(upp(:) - phi_trial(:))) > 0
        l1 = 0; l2 = 1e6; phi_feas = a_lim;
        for it = 1:25
            lmid = 0.5 * (l1 + l2);
            pj = p0 + lmid * p1;
            phi_trial = max(a_lim, min(b_lim, (sqrt(pj).*low + sqrt(q0).*upp)./(sqrt(pj) + sqrt(q0))));
            if r1 + sum(p1(:)./(upp(:) - phi_trial(:))) > 0
                l1 = lmid; else, l2 = lmid; phi_feas = phi_trial;
            end
        end; phi_trial = phi_feas;
    end
    ch = max(abs(phi_trial(:) - phi(:)));
    % Display Iteration Status & Topology
    fprintf('It.:%3d  Comp.:%7.3f  ||h||_*:%7.2e  ch.:%6.4f\n', loop, c, h_norm, ch);
    colormap(gray); imagesc(1 - x); caxis([0 1]);
    axis equal; axis tight; axis off; drawnow;
    if loop >= 30 && ch < 0.01, break; end
    phi = phi_trial;
    % ALM Multiplier & Penalty Schedule
    lambda_h = lambda_h + rho_h * (h .* Minv);
    lambda_h = max(-10, min(10, lambda_h));
    lambda_h = lambda_h - sum(M(:).*lambda_h(:)) / Omega;
    if mod(loop, 10) == 0 && h_norm > 1e-3, rho_h = min(10.0, rho_h * 1.1); end
    c_old = c;
end
end

%% ===========================================================================
%% topAC99: Allen-Cahn Topology Optimization (ACTO) via ALM
%% ===========================================================================
% Author: Run Du (rdu@swjtu.edu.cn)
% Repository: https://github.com/smolab/topAC99

% Copyright (c) 2026 Run Du (Southwest Jiaotong University, China)
%
% This program is free software: you can redistribute it and/or 
% modify it under the terms of the GNU Affero General Public License (AGPLv3)
% as published by the Free Software Foundation, either version 3 of the License, 
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Affero General Public License for more details.
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%% ===========================================================================
% Usage:
%   [c, x, loop, h_norm] = topAC99(120, 80, 0.30, 3.0, 2.1, 0.065);
%
% Inputs:
%   nelx, nely : Number of elements in horizontal and vertical directions
%   volfrac    : Prescribed volume fraction limit (e.g., 0.30 or 0.40)
%   penal      : SIMP stiffness penalization exponent (baseline: 3.0)
%   w          : Diffuse interface width in element units (baseline: 2.1)
%   beta       : Double-well potential barrier height (baseline: 0.065)
%
% Outputs:
%   c          : Optimized structural compliance (F^T * U)
%   x          : Physical element density matrix (nely x nelx)
%   loop       : Number of iterations executed until convergence
%   h_norm     : Mass-weighted continuous dual RMS residual norm ||h||_*
%% ===========================================================================
