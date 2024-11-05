1;

search_algorithms;

function [xmin, fmin, nbiter, iters, CONVCRIT]=minimize(x0, f, gr, varargin)
  % Default parameters
  tol = 0.01;
  iterlimit = 400;
  dkeps = 1e-10;
  flagx = false;
  flag_alpha = 1;
  flag_beta = 0;
  flag_dfp = false;
  flag_bfgs = false;

  iter.x = [];
  iter.f = [];
  iter.alpha = [];
  iter.beta = [];
  iter.dk = [];
  iter.gr = [];
  iter.s = [];
  iter.S = [];

  DIM = rows(x0);

  % Parse number of output arguments
  if nargout==4
    flagx = true;
  endif

  % Parse input arguments
  k = 1;
  if nargin>3

    % First optional input
    if rem(nargin,2)==0
      method=varargin(1);
      k+=1;
    endif

    % Argument pairs to be parsed
    while k<length(varargin)
      s = cell2mat(varargin(k));
      k+=1;
      a = cell2mat(varargin(k));
      k+=1;
      switch (s)
        case 'alphamethod'
          switch a;
            case 'aramijo'
              flag_alpha = 1;
            case 'armijo'
              flag_alpha = 1;
            case 'parabolic'
              flag_alpha = 2;
          endswitch
        case 'betamethod'
          switch a;
            case 'none'
              flag_beta = 0;
            case 'fletcher'
              flag_beta = 1;
            case 'ribière'
              flag_beta = 2;
          endswitch
        case 'tol'
          tol = a;
        case 'iterlimit'
          iterlimit = a;
        case 'dkeps'
          dkeps = a;
        case 'newtonmethod'
          switch a
            case 'dfp'
              flag_dfp = true;
            case 'bfgs'
              flag_bfgs = true;
          endswitch
      endswitch
    endwhile
  endif

  % First iteration

  xk = x0;
  grk = gr(xk);
  ngrk = norm(grk);
  if ngrk < dkeps
    CONVCRIT = "Norm of the gradient too small";
  endif

  dk = -grk;
  Sk = eye(DIM);
  dk = Sk*dk;

  if flag_alpha==1
    alphak = aramijo_alpha_search(xk, dk, grk, f);
  else
    alphak = parabolic_alpha_search(xk, dk, grk, f);
  endif

  sk = alphak*dk;
  xk1 = xk + sk;
  fk1 = f(xk1);
  fk = f(xk);

  iter.x = [xk];
  iter.f = [fk];
  iter.alpha = [alphak];
  iter.beta = [0];
  iter.dk = [dk];
  iter.gr = [grk];
  iter.s = [sk];
  iter.S = [Sk];
  iters = iter;

  nbiter = 1;
  precrel = norm(xk1 - xk)/norm(xk);

  % Iteration loop
  while precrel > tol && nbiter < iterlimit
    grk1 = gr(xk1);

    Sk1 = Sk;
    betak1 = 0;
    dk1 = -grk1;
    if flag_beta
      if rem(nbiter, DIM)
        if flag_beta==1
          betak1 = beta_keanureeves_search(grk, grk1);
        elseif flag_beta==2
          betak1 = beta_biere_search(grk, grk1);
        end
        dk1 = -grk1 + betak1*dk;
      endif
    elseif flag_dfp
      gama_k = grk1 - grk;
      delta_k = xk1 - xk;
      vk = delta_k - Sk*gama_k;
      ak = 1/(vk'*gama_k);
      Ck = ak*vk*vk';

      Sk1 = Sk + Ck;
      dk1 = Sk1*dk1;
    elseif flag_bfgs
      gama_k = grk1 - grk;
      delta_k = xk1 - xk;
      Ck = (1 + gama_k' * Sk * gama_k / (delta_k' * gama_k));
      Ck -= (delta_k * gama_k' * Sk + Sk * gama_k * delta_k') / (delta_k' * gama_k);

      Sk1 = Sk + Ck;
      dk1 = Sk1*dk1;
    endif

    ndk1 = norm(dk1);
    if ndk1 < dkeps
      CONVCRIT = "Step direction tolerance achieved";
      break;
    endif
    dk1 /= ndk1;

    if flag_alpha==1
      alphak1 = aramijo_alpha_search(xk1, dk1, grk1, f);
    elseif flag_alpha==2
      alphak1 = parabolic_alpha_search(xk1, dk1, grk1, f);
    endif

    sk1 = alphak1*dk1;
    xk2 = xk1 + sk1;
    precrel = norm(xk2 - xk1)/norm(xk1);
    if precrel <= tol
      CONVCRIT = "Convergence tolerance achieved";
    endif

    iter.x = xk1;
    iter.f = fk1;
    iter.alpha = alphak1;
    iter.beta = betak1;
    iter.dk = dk1;
    iter.gr = grk1;
    iter.s = sk1;
    iter.S = Sk1;
    iters(end+1) = iter;

    xk = xk1;
    xk1 = xk2;
    dk = dk1;
    fk1 = f(xk1);
    fk = f(xk);
    grk = grk1;
    nbiter += 1;
    if nbiter >= iterlimit
      CONVCRIT = "Maximum number of iterations achieved";
    endif
  endwhile

  xmin = xk2;
  fmin = f(xk2);

  iter.x = xk2;
  iter.f = fmin;
  iter.alpha = [];
  iter.beta = [];
  iter.dk = [];
  iter.gr = [];
  iter.s = [];
  iter.S = [];
  iters(end+1) = iter;
end
