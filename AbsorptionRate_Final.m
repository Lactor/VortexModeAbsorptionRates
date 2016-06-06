function [Gamma] = AbsorptionRate_Final...
                           (ni,li,mi, nf,lf,mf, displacement, s, R,  lvortex)

   % By Francisco Machado (fmachado@mit.edu)
   
   % Calculation of the interaction matrix element between state ni,li,mi, and
   % nf,lf,mf in hydrogen, through a plasmon mode with confinement
   % factor s.
   %
   % displacement [dx, dy,dz] - displacement in nm of the atom with respect to the
   % vortex mode. dz = 0 implies the atom is at the surface of the
   % SP(h)P supporting material
   displacement = displacement*10^-9; %Transforming into meters
   %
   % lvortex - the OAM value of the vortex mode of interest, use -1
   % for the plane wave SP(h)P mode.
   %
   % R is the upper integration limit in the radial integration. A
   % value of 300 was found to be enough for the convergence of the value.
   
   
   %% Defining x = r/a0 the dimensionless distance
   
   %% Defines fundamental constants
   % Compute k to match the energy of the transition
   Ryd = 1.0973731569 * 10^7; % (m^(-1))
   hbar = 1.05*10^(-34); %(J.s)
   c = 3*10^8; % (m/s)
   a0 = 5.26 * 10^(-11); % (m)
   e = 1.6*10^(-19); % (C)
   m = 9.11 * 10^(-31); % (kg)
   alpha = (137.036)^(-1); 
   epsilonbar = (1+3.9)/2; % Relative average between vaccum/air                       
   % and SiO2 (source: Analysis and Design of Analog Integrated
   % Circuits 2009 p. 156)

   %% Transition matching momentum and energy
   k = 2*pi * Ryd * abs( 1/nf^2 - 1/ni^2);
   omega = c*Ryd*abs(1/ni^2 - 1/nf^2);

   %% Defines derivative Step
   delta = 10^(-5);
   
   
   %% Defines Plane Wave vector potential
   % In the electrostatic approximation q=K we have
   % that the plasmon propagating in x direction is given by:
   %    
   % A = e^{-ksz + iksr} [\hat{x} - i\hat{z}]/\sqrt{2}
   % 
   VectorPotX = @(x, theta,phi) (1/sqrt(2)) * exp(-(k*s*a0).*x.*cos(theta) -1i.*(s*k*a0).*x.*sin(theta).*cos(phi)) ...
       .* (1i.*cos(theta) - sin(theta).*cos(phi));
   VectorPotT = @(x, theta,phi) (1/sqrt(2)) * exp(-(k*s*a0).*x.*cos(theta) -1i.*(s*k*a0).*x.*sin(theta).*cos(phi)) ...
       .* (-1i.*sin(theta) - cos(theta).*cos(phi));
   VectorPotP = @(x, theta,phi) (1/sqrt(2)) * exp(-(k*s*a0).*x.*cos(theta) -1i.*(s*k*a0).*x.*sin(theta).*cos(phi)) ...
       .*  sin(phi);

   
   %% OAM carrying vortex mode vector potential defined in
   %% functions at the end of the file
   
   %% Define Electron Wavefunction in Hydrogren
   N = @(n,l,m) (-1)^m * 1/n^2* 2^(l+1) * sqrt( (factorial(n-l-1)*(2*l+1)*factorial(l-m))/...
                    (factorial(n+l) * 4*pi*factorial(l+m)));
   Psi = @(n,l,m, x,theta,phi) N(n,l,m) * exp(-x./n) .* (x./n).^l .* polyval(LaguerreGen(n-l-1, 2.*l+1), 2.*x./n)...
       .*exp(1i.*m.*phi) .* legendre_newP(l,m,cos(theta));
    
   
   %% Define Grad of Wavefunction
   gradX = @(n,m,l,x,theta,phi) (Psi(n,m,l,x+delta, theta,phi) -Psi(n,m,l,x,theta,phi))./(delta);
   gradT = @(n,m,l,x,theta,phi) 1./x .* (Psi(n,m,l,x, theta+delta,phi) -Psi(n,m,l,x,theta,phi))./(delta);
   gradP = @(n,m,l,x,theta,phi) 1./(x.*sin(theta)).*(Psi(n,m,l,x, ...
                                                     theta,phi+delta) -Psi(n,m,l,x,theta,phi))./(delta); 
   
   %% Dimension of wavefunction is canceled by the volume
   %% integration, so all quantities are kept dimensionlesss

    
   
   %% Confirmation of wavefunction definition
   % for n=1:8
   %     for l=0:(n-1)
   %         for ml=-l:l
   %             %n,l,m
   %             integrand = @(x, theta, phi)  x.^2 .* sin(theta) .* conj(Psi(n,l,ml, x, theta, phi) ) .* Psi(n,l,ml, x, theta, phi);
   %             integral3(integrand, 0, 300, 0, pi, 0, 2*pi)
   %         end
   %     end
   % end

   
   %% Computed Transition for Plane Wave SPP
   if (lvortex == -1)

       %% Computes the matrix element through a simple 3
       %% dimensional integral. The integrand is computed as the
       %% inner product of the vector potential and the gradient
        integrand = @(x,theta,phi) - x.^2.*sin(theta).*conj(Psi(nf,lf,mf, x, theta, phi)).* (1i ) .* ...
            (VectorPotX(x,theta,phi) .* gradX(ni,li,mi, x,theta,phi) + ...
             VectorPotT(x,theta,phi) .* gradT(ni,li,mi, x,theta,phi) + ...
             VectorPotP(x,theta,phi) .* gradP(ni,li,mi, x,theta,phi));
        
        MatElement = integral3(integrand, 0, R, 0, pi, 0, 2*pi, 'AbsTol',1e-15, 'RelTol', 1e-3 );
        
   %% Computes Transition for Vortex Mode for centered Atom
   elseif lvortex > -1 && (sqrt(displacement(1)^2 + displacement(2)^2)) ...
           < 10^-15

       % Due to the rotational symmetry of the problem the
       % azimuthal integration can be done analytical and yield 2pi
       % when lvortex = mf - mi, and 0 otherwise. The integral is
       % then computed using the value of phi = 0, although any
       % would yield the same answer, after absolute value of the
       % matrix is taken.
       
       if mi + lvortex - mf == 0
           PotX = @(x,theta) VortexlSpX(s, k, lvortex, x, theta, 0);
           PotT = @(x,theta) VortexlSpT(s, k, lvortex, x, theta, 0);
           PotP = @(x,theta) VortexlSpP(s, k, lvortex, x, theta, 0);
           
           integrand = @(x,theta) x.^2.*sin(theta).*conj(Psi(nf,lf,mf, x, theta, 0)).* (1i ) .* ...
               (PotX(x,theta) .* gradX(ni,li,mi, x,theta,0) + ...
                PotT(x,theta) .* gradT(ni,li,mi, x,theta,0) + ...
                PotP(x,theta) .* gradP(ni,li,mi, x,theta,0)); 
           MatElement = 2*pi * quad2d(integrand, 0, R, 0, pi, ...
                                      'AbsTol',1e-15, 'RelTol', ...
                                      1e-3, 'MaxFunEvals', 50000);
       else
           MatElement = 0;
       end
       
   %% Computes Transition for Not Centered Atom (by doing 3d integral)
   elseif lvortex > -1

       %% Computes the vector potential of the displaced vortex
       %% mode in spherical coordinates centered at the vortex
       %% center
       
       PotXPrime = @(x,theta, phi) VortexLDispX(s, k, lvortex,displacement(1), ...
                                                displacement(2), x, theta, phi);
       PotTPrime = @(x,theta, phi) VortexLDispT(s, k, lvortex,displacement(1), ...
                                                displacement(2), x, theta, phi);
       PotPPrime = @(x,theta, phi) VortexLDispP(s, k, lvortex,displacement(1), ...
                                                displacement(2), x, theta, phi);
       
       %% In order to compute the inner product between the
       %% gradient and the vector potential it is necessary to
       %% relate the two sets of basis vectors through a change of
       %% basis matrix

       
       % Coordinates with respect to the electronic system
       X = @(x,theta,phi) a0*x.*sin(theta).*cos(phi);
       Y = @(x,theta,phi) a0*x.*sin(theta).*sin(phi);
       Z = @(x,theta,phi) a0*x.*cos(theta);
       
       r = @(x,theta,phi) a0*x;
       rho = @(x,theta,phi) a0*x.*sin(theta);

       % Coordinates with respect to the vortex center
       Xprime = @(x,theta,phi) X(x,theta,phi) -displacement(1);
       Yprime = @(x,theta,phi) Y(x,theta,phi) -displacement(2);       
       Zprime = @(x,theta,phi) Z(x,theta,phi);
       
       rhoprime = @(x,theta,phi) sqrt(Xprime(x,theta,phi).^2 + ...
                                      Yprime(x,theta,phi).^2);
       Rprime = @(x,theta,phi) sqrt(rhoprime(x,theta,phi).^2 + ...
                                    Zprime(x,theta,phi).^2);

       % C - Cosine ; S - Sine
       % TPrime - ThetaPrime  ; PPrime - PhiPrime
       % Computes the trigonometric function of angles in the
       % vortex coordinate system
       CTprime = @(x,theta,phi) Zprime(x,theta,phi)./Rprime(x, ...
                                                         theta,phi);
       STprime = @(x,theta,phi) rhoprime(x,theta,phi)./ Rprime(x, ...
                                                         theta,phi);
       CPprime = @(x,theta,phi) Xprime(x,theta,phi) ./ rhoprime(x, ...
                                                         theta,phi);
       SPprime = @(x,theta,phi) Yprime(x,theta,phi) ./ rhoprime(x, ...
                                                         theta,phi);
       

       % R - \hat{r} ; RPrime - \hat{rprime} and analogous to
       % T - Theta and P - Phi

       % Computes the inner product between the basis vectors of
       % the two coordinate system forming a 3x3 change of basis matrix
       
       RRprime = @(x,theta,phi) ...
                 (sin(theta).*cos(phi).*STprime(x,theta,phi).*CPprime(x,theta,phi) + ...
                  sin(theta).*sin(phi).*STprime(x,theta,phi).*SPprime(x,theta,phi) + ...
                  cos(theta).*CTprime(x,theta,phi));
       
       RTprime =  @(x,theta,phi) ...
           (sin(theta).*cos(phi).*CPprime(x,theta,phi).*CTprime(x,theta,phi) + ...
            sin(theta).*sin(phi).*SPprime(x,theta,phi).*CTprime(x,theta,phi) + ...
            - cos(theta).*STprime(x,theta,phi));
       
       RPprime =  @(x,theta,phi) ... 
           (-sin(theta).*cos(phi).*SPprime(x,theta,phi) + ...
            sin(theta).*sin(phi).*CPprime(x,theta,phi));
       
       TRprime =  @(x,theta,phi) ... 
           (cos(phi).*cos(theta).*STprime(x,theta,phi).*CPprime(x,theta,phi) + ...
            sin(phi).*cos(theta).*STprime(x,theta,phi).*SPprime(x,theta,phi) + ...
            - sin(theta).*CTprime(x,theta,phi));
       
       TTprime =  @(x,theta,phi) ... 
           (cos(phi).*cos(theta).*CPprime(x,theta,phi).* ...
            CTprime(x,theta,phi) + ...
            sin(phi).*cos(theta).*SPprime(x,theta,phi).* ...
            CTprime(x,theta,phi)+ ...
            sin(theta).*STprime(x,theta,phi));
       
       TPprime =  @(x,theta,phi) ... 
           (-cos(phi).*cos(theta).*SPprime(x,theta,phi)+...
            sin(phi).*cos(theta).*CPprime(x,theta,phi));
       
       
       PRprime =  @(x,theta,phi) ... 
           (-sin(phi).*STprime(x,theta,phi).*CPprime(x,theta,phi) +...
            cos(phi).*STprime(x,theta,phi).*SPprime(x,theta,phi));
       
       PTprime =  @(x,theta,phi) ... 
           (-sin(phi).*CPprime(x,theta,phi).* ...
            CTprime(x,theta,phi) + ...
            cos(phi).*SPprime(x,theta,phi).* ...
            CTprime(x,theta,phi));
       
       PPprime =  @(x,theta,phi) ... 
           (sin(phi).*SPprime(x,theta,phi) +...
            cos(phi).*CPprime(x,theta,phi));
       
       
       
       %% Defines the integrand of our 3 dimensional integral as
       %% the inner product of the two vectors, computed with the
       %% change of basis matrix.
       
       integrand = @(x,theta,phi) x.^2.*sin(theta).*conj(Psi(nf,lf,mf, x, theta, phi)).* (1i ) .* ...
            (gradX(ni,li,mi, x,theta,phi) .* ...
                (PotXPrime(x,theta,phi).*RRprime(x,theta,phi) +...
                 PotTPrime(x,theta,phi).*RTprime(x,theta,phi) +...
                 PotPPrime(x,theta,phi).*RPprime(x,theta,phi))+ ...
              gradT(ni,li,mi, x,theta,phi) .* ...
                 (PotXPrime(x,theta,phi).*TRprime(x,theta,phi) +...
                  PotTPrime(x,theta,phi).*TTprime(x,theta,phi) +...
                  PotPPrime(x,theta,phi).*TPprime(x,theta,phi))+ ...
            gradP(ni,li,mi, x,theta,phi) .* ...
                (PotXPrime(x,theta,phi).*PRprime(x,theta,phi) +...
                 PotTPrime(x,theta,phi).*PTprime(x,theta,phi) +...
                 PotPPrime(x,theta,phi).*PPprime(x,theta,phi)));

       % Computes the integral.
       MatElement = integral3(integrand, 0, R, 0, pi, 0, 2*pi, ...
                            'AbsTol',1e-8, 'RelTol', 1e-5);
   else
       disp('ERROR The case is Incorrect')
   end
   
  
   %% Transition Rate computation
   % Simple formula arises from the field normalization.
   Gamma = (2*pi/epsilonbar) * s^3 * omega * alpha^3 * abs(MatElement)^2 ...
       * exp(-2*k*s*displacement(3));
   
end


function [Ar] = VortexlSpX(s, k, l,  x, theta, phi)
%% Compute the radial component of the vortex mode vector potential
% r - radial component
% theta - angle with z axis 
% phi - azimuthal angle
    
    % Bohr Constant 
    a0 = 5.26 * 10^(-11);

    % Recover the length scale of the problem
    r = x*a0;
    K = k*s;
    Ar = 1/(2*pi)* 1/sqrt(2) * (1i)^(1-l) * pi * exp(-K*r.*cos(theta)) .* exp(1i*l .*phi) ...
         .*(2*besselj(-l, K*r.*sin(theta)) .*cos(theta) + ...
            (besselj(1-l, K*r.*sin(theta)) - besselj(-1-l, K*r.*sin(theta))).*sin(theta));

end 


function [Ar] = VortexLDispX(s, k, l,x0,y0,  x, theta, phi)

%% Compute the radial component of the displaced vortex mode
%% vector potential
    
% r - radial component on the atom basis
% theta - angle with z axis on the atom basis
% phi - azimuthal angle on the atom basis

% x0 - displacement in m
% y0 - displacement in m
% z0 - displacement in m is treated at the end

    % Bohr Constant 
    a0 = 5.26 * 10^(-11);
    
    % Recover the length scale of the problem
    K = k*s;
    r = x*a0;

    x = r.*sin(theta).*cos(phi) - x0;
    y = r.*sin(theta).*sin(phi) - y0;
    z = r.*cos(theta);

    rho = sqrt(x.^2 + y.^2);
    S = rho./sqrt(z.^2 + rho.^2);
    C = z./sqrt(z.^2 + rho.^2);

    Ar =  1/(2*pi)* 1i^(1-l) * pi/ sqrt(2) * ((x+1i*y)./rho).^l .* exp(-K*z) ....
          .*(2*besselj(-l, K*rho) .* C + (besselj(1-l, K*rho) - ...
                                          besselj(-1-l, K*rho)) .*S );
end 


function [At] = VortexlSpT(s, k, l, x, theta, phi)
%% Compute the \hat{theta} component of the vortex mode vector potential
    
% See VortexlSpX
    K = s*k;
    a0 = 5.26 * 10^(-11);
    r = x*a0;
    At =- 1/(2*pi)* 1/sqrt(2) * (1i)^(1-l) * pi * exp(-K*r.*cos(theta)) .*  exp(1i*l .*phi) .*...
        ( 2*besselj(-l, K*r.*sin(theta)) .*sin(theta) + ...
          (besselj(-l-1, K*r.*sin(theta)) - besselj(1-l, K*r.*sin(theta))).*cos(theta));

end 

function [At] = VortexLDispT(s, k, l,x0,y0, x, theta, phi)
%% Compute the \hat{theta} component of the displaced vortex mode vector potential

% For information see VortexLDispX
    K = s*k;
    a0 = 5.26 * 10^(-11);
    r = x*a0;

    x = r.*sin(theta).*cos(phi) - x0;
    y = r.*sin(theta).*sin(phi) - y0;
    z = r.*cos(theta);

    rho = sqrt(x.^2 + y.^2);
    S = rho./sqrt(z.^2 + rho.^2);
    C = z./sqrt(z.^2 + rho.^2);

    At =   - 1/(2*pi)*1i^(1-l) * pi/ sqrt(2) * ((x+1i*y)./rho).^l .* exp(-K*z) ....
           .*(2*besselj(-l, K*rho) .* S + ( besselj(-l-1, K*rho) - ...
                                            besselj(1-l, K*rho)) .*C );
end 

function [Ap] = VortexlSpP(s, k,l,  x, theta, phi)
%% Compute the azimuthal component of the vortex mode vector potential
% For information see VortexSp
    
    K = k*s;
    a0 = 5.26 * 10^(-11);
    r= x*a0;
    Ap= 1/(2*pi)* sqrt(2) * (1i)^(-l) * pi * l * exp(-K*r.*cos(theta)) .*  exp(1i*l*phi) ...
        .* (besselj(-l, K*r.*sin(theta)).* csc(theta))./(K*r);

end 

function [Ap] = VortexLDispP(s, k,l,x0,y0,  x, theta, phi)
%% Compute the azimuthal component of the displaced vortex mode vector potential
% For information see VortexLDispX
    
    K = s*k;
    a0 = 5.26 * 10^(-11);
    r = x*a0;

    x = r.*sin(theta).*cos(phi) - x0;
    y = r.*sin(theta).*sin(phi) - y0;
    z = r.*cos(theta);

    rho = sqrt(x.^2 + y.^2);
    S = rho./sqrt(z.^2 + rho.^2);
    C = z./sqrt(z.^2 + rho.^2);

    Ap= 1/(2*pi)* sqrt(2) * (1i)^(-l) * pi * l * ((x+1i*y)./rho).^l .* exp(-K*z) .* ...
        (besselj(-l, K*rho))./(K*rho);
end 
