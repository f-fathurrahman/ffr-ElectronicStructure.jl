function [Etotvec, X, vtot, rho] = chebyscf(mol,varargin)
%
%  Usage: [Etotvec, X, vtot, rho] = chebyscf(mol);
%
%         [Etotvec, X, vtot, rho] = chebyscf(mol,'maxscfiter', 15); 
%
%         [Etotvec, X, vtot, rho] = chebyscf(mol,'scftol', 1e-6, 'maxscfiter',15); 
%
%         [Etotvec, X, vtot, rho] = chebyscf(mol,'scftol', 1e-6, 'degree',5, ...
%                                            'X0',X0); 
%         [Etotvec, X, vtot, rho] = chebyscf(mol,'verbose','off');
%
%  Purpose:
%     Use Self Consistent Field (SCF) iteration to find the ground state 
%     minimum total energy and the corresponding wave functions.
%     A Chebyshev polynomial filter is used in each SCF iteration
%     to construct an approximate invariant subspace associated
%     with the nocc smallest eigenvalues of the KS Hamiltonian.
%
%  Input:
%     mol                                a Molecule object
%     (option_name,option)value) pairs   option_name is a predefined string,
%                                        option_value is a user provided value
%     The following option_names and option_values are allowed:
%
%   option_name        option_value type            purpose
%   -----------        -----------------            -------
%   'maxscfiter'          integer        maximum # of inner SCF iterations allowed 
%                                        (default 15)
%
%   'degree'              integer        the degree of the Chebyshev polynomial 
%                                        filter (default 10)
%
%   'scftol'             float            residual tolerance for SCF
%   'X0'                 Wavefun object   initial guess of the wavefunctions
%   'verbose'            'on' or 'off'    turn on/off print statements.
%
%  Output:
%   Etotvec     an array of total energy values computed at each DCM iteration
%   X           final approximate wavefunctions (Wavefun object)
%   vtot        final total (local) potential (3D array)
%   rho         final total charge density (3D array)
%


%
% check the mol input argument is valid.
%

% set defaults
maxscfiter = 10;
scftol     = 1e-8;
maxcgiter  = 10;
X0         = [];
X          = [];
vtot       = [];
rho        = [];
pmix       = 1; % turn Pulay mixing on by default
kmix       = 1; % turn Kerker mixing on by default
Etotvec    = [];
degree     = 10;
verbose    = 1;
tstart     = cputime;
%
% parse input options and change defaults if necessary
%
options = varargin;
while (length(options) >= 2)
optname = options{1};
value   = options{2};
options = options(3:end);
switch optname
    case 'maxscfiter'
        maxscfiter = value;  
    case 'scftol'
        scftol = value;
    case 'maxcgiter'
        maxcgiter = value; 
    case 'X0'
        X0 = value; 
    case 'degree'
        degree = value; 
    case 'pulaymix'
        if (strcmp(value,'On') || ...
	    strcmp(value,'ON') || ...
	    strcmp(value,'on') || ...
	    strcmp(value,'oN') )
	    pmix = 1;
        else
	    pmix = 0;
	end;
    case 'kerkermix'
	if ( strcmp(value,'On') || ...
	     strcmp(value,'ON') || ...
	     strcmp(value,'on') || ...
	     strcmp(value,'oN') )
	    kmix = 1;
        else
	    kmix = 0;
	end;
    case 'verbose'
	if ( strcmp(value,'Off') || ...
             strcmp(value,'OFF') || ...
	     strcmp(value,'off') || ...
	     strcmp(value,'oFF') )
            verbose = 0;
        else
	    verbose = 1;
	end;
    end; 
end;
%
% we do not distinguish spins (nspin = 1), i.e.
% the number of occupied state = total number of valence electrons/2
%
nspin = get(mol,'nspin');
nocc  = getnel(mol)/2*nspin;
if ( isempty(X0) )
   %
   % generate at least two more wavefunctions
   % these extra wavefunctions are needed for computing
   % bounds for an interval in which the Chebyshev
   % polynomial has a damping effect.
   %
   X0 = genX0(mol,nocc+2);
end;
%
% check the dimension of X0
%
if ( get(X0,'ncols') < nocc )
   fprintf('Error: The number of columns in X0 is less than nocc\n');
   fprintf('Error: size(X0,2) = %d, nocc = %d\n', get(X0,'ncols'), nocc);
   return;
end;
%
n1 = get(mol,'n1');
n2 = get(mol,'n2');
n3 = get(mol,'n3');
n123 = n1*n2*n3;
m1 = get(X0,'n1');
m2 = get(X0,'n2');
m3 = get(X0,'n3');
if (m1 ~= n1 | m2 ~= n2 | m3 ~= n3)
   error('scf: dimension of the molecule does not match that of the wavefunction')
end;
%
iterscf = 1;
%
% construct the initial Hamiltonian
%
H    = Ham(mol);
vion = get(H,'vion');
vext = get(H,'vext');
vin  = get(H,'vtot');
%
% calculate Ewald and Ealphat (one time calculation)
%
Ewald     = getewald(mol);
Ealphat   = getealphat(mol);
%fprintf('Ewald   = %11.4e\n', Ewald);
%fprintf('Ealphat = %11.4e\n', Ealphat);
%
% run 2*nocc steps of Lanczos iteration to obtain
% the lower and upper bounds of the interval in
% which the eigenvalues are filtered out.
%
v0 = genX0(mol,1);
[T,Q,f]=lanczos(H,v0,nocc*2);
d = sort(real(eig(T)));
lb = d(nocc+2);
ub = d(2*nocc);
%
fprintf('Beging ChebySCF calculation...\n');
V = [];
R = [];
while ( iterscf <= maxscfiter )
    %
    % apply Chebyshev filter polynomial to the wavefunctions contained in X0
    %
    X = chebyfilt(H, X0, degree, lb, ub);
    %
    % orthonormalize the new basis vectors
    %
    [Q,R0]=qr(X,0);
    %
    % compute Rayleigh quotient
    %
    T = Q'*(H*Q); T = (T+T')/2;
    [U,D] = eig(T);
    X = Q*U;
    [ev,id] = sort(real(diag(D)));
    X = X(:,id);
    %
    % save potentials for charge mixing
    %
    V = [V reshape(vin, n123, 1)];
    %
    % compute charge density
    % 
    rho = getcharge(mol,X(:,1:nocc),nspin);
    %
    % Kinetic energy and some additional energy terms 
    %
    Ekin = (2/nspin)*sum(ev(1:nocc));
    %
    % ionic and external potential energy was included in Ekin
    % along with incorrect Ecoul and Exc. Need to correct them
    % later;
    %
    Ecor = getEcor(mol, rho, vin, vion, vext);
    %
    % Compute Hartree and exchange correlation energy and potential
    % using the new charge density; update the total potential
    %
    [vhart,vxc,uxc2,rho]=getvhxc(mol,rho);
    vout = getvtot(mol, vion, vext, vhart, vxc);
    %
    % Calculate the potential energy based on the new potential
    %
    Ecoul = getEcoul(mol,abs(rho),vhart);
    Exc   = getExc(mol,abs(rho),uxc2);
    %
    %fprintf('Ekin  = %22.12e\n', Ekin);
    %fprintf('Ecor  = %22.12e\n', Ecor);
    %fprintf('Ecoul = %22.12e\n\n', Ecoul);
    %fprintf('Exc   = %22.12e\n\n', Exc);
    %
    Etot = Ewald + Ealphat + Ekin + Ecor + Ecoul + Exc;
    Etotvec = [Etotvec Etot];
    %
    % check self-consistency in potential
    %
    r = reshape(vout-vin, n123, 1);
    verr = norm(r); 
    fprintf('\nSCF iter %2d:\n', iterscf);
    if (verbose == 1) 
       fprintf('norm(vout-vin) = %10.3e\n', verr);
    end
    %
    fprintf('Total energy   = %20.13e\n', Etot);
    if (verr < scftol)
       % converged
       break;
    end;
    %
    % perform charge mixing to update the potential
    % 
    vupdated = 0;
    R = [R r];
    if (iterscf > 1 && pmix == 1) 
       vin = pulaymix(vin,vout,V,R);
       vupdated = 1;
    end
    if (kmix == 1)
       vin = kerkmix(mol,vin,vout);
       vupdated = 1;
    end
    if (~vupdated)
       vin = vout; 
    end
    %
    % update the total potential and hence the Hamiltonian
    %
    H = set(H,'vtot',vin);
    %
    % check residual error
    %
    HX = H*X;
    G = X'*HX; G = (G+G')/2;
    RES = HX-X*G;
    if (verbose==1)
       for j = 1:nocc
          fprintf(' resnrm = %11.3e\n', norm(RES(:,j)));
       end
       fprintf('-------------\n\n');
    end
    %
    X0 = X;
    %
    % update lb and ub
    %
    lb = ev(nocc+2);
    [T,W,f]=lanczos(H,v0,10);
    ub = norm(T) + norm(f);
    iterscf = iterscf + 1;
    %pause;
end;
vtot = vin;
timetot = cputime - tstart;
fprintf('Total time used = %11.3e\n', timetot);
fprintf('norm(HX-XD)     = %11.3e\n', norm(RES,'fro'));
rho = fftshift(real(rho));
vtot = fftshift(real(vtot));
