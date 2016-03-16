# MPC-R
# Model Predictive Control (MPC) toolbox in R
# Examples, solvers and support functions
# Run e.g. MPCexample1()
# Requires packages: expm, Matrix, lpSolve, quadprog
# Rasmus Halvgaard, rhal@dtu.dk, 2016-3-10

# Example of Economic MPC and a first order system
MPCexample1 <- function() {
  Np <- 50 # [samples] Prediction horizon
  Ns <- 50 # [samples] Simulation horizon
  
  # Example simulation scenario data
  d <- rep(0,max(Np,Ns+Np)) # Disturbance
  p <- rep(0.9,max(Np,Ns+Np)) # Price
  p[30:35] <- 0.2
  ymin <- rep(1,Ns+Np) # minimum output constraint
  ymin[15:20] <- 2
  u <- rep(0,Ns+1) # previous control action u[0]
  y <- rep(-0.5,Ns+1) # intial output
  
  # Predictive model used in MPC
  tau <- 10; K <- 5
  H <- ImpulseResponse(A=-1/tau,B=K/tau,E=0,C=1,Np,Ts=1)
  dH <- DifferenceMatrix(Np)
  
  # Closed-loop MPC
  for (k in 1:Ns) {
    # Update predictions and constraints
    dp <- d[k:(k+Np-1)] # disturbance prediction
    pp <- p[k:(k+Np-1)] # price prediction
    yminp <- ymin[k:(k+Np-1)] # min output constraint
    c <- SetConstraints1(Np,H,dH,dp,x0=y[k],um1=u[k],yminp)
    
    mpc <- BuildMPClp(c,Np,pp)
    opt <- SolveMPClp(mpc)
    
    if (opt$status) {
      warning("Infeasible")
    } else {
      u[k+1] <- opt$u[1] # only use first control action
      y[k+1] <- opt$y[2] # save output
    }
  }
  
  if (Ns==1) { # Save open-loop profile
    u <- opt$u
    y <- opt$y
  } else {
    u <- tail(u,-1) # discard u[0]
  }
  
  # Results
  # The Economic MPC minimizes costs and maintains operating constraints.
  # (Time)  
  # (0:5)   The initial condition is lower than the output constraint.
  #         So the controller brings it back to minimum.
  #         The slack variable makes sure that the problem is solvable
  #         even though constraints are violated.
  # (15:20) The minimum output constraint is increased and controller
  #         reacts ahead to be above the constraint.
  # (20:25) The control input is u=0, but the output slowly decreases
  #         accoring to the first order dynamics.
  # (30:35) The input price is cheap and the controller increases u and
  #         consequently also y that still stays within constraints.
  plot(y,ylab="(u,y,p)",xlab="Time",type="s",cex.lab=1.5,cex.axis=1.5)
  lines(u,type="s",col="blue")
  lines(p,col="red",type="s")
  lines(c(NA,ymin),lty=2,type="s")
  lines(rep(c$umin,Ns),lty=2,col="blue",type="s")
  lines(rep(c$ymax,Ns),lty=2,type="s")
  legend(0,2.5,c("y","u","p","umin","ymin"),lty=c(1,1,1,2,2),lwd=2.5,
         col=c("black","blue","red","blue","black"),cex=1.2)
}
SetConstraints1 <- function(Np,H,dH,d,x0,um1,ymin=rep(1,Np),du=2) {
  # Setup MPC constraints
  umin <- rep(0,Np)
  umax <- rep(1,Np)
  dumin <- -umax/du
  dumax <- umax/du
  ymin <- ymin
  ymax <- rep(4,Np)
  dI0 <- rep(0,Np); dI0[1] <- um1
  
  return(constraints = list(um1=um1,x0=x0,y0=x0,umin=umin,umax=umax,dumin=dumin,dumax=dumax,
                            ymin=ymin,ymax=ymax,dH=dH,dI0=dI0,H=H,d=d))
}

# Example of Economic MPC and a storage unit
MPCexample2 <- function() {
  Np <- 12 # [samples] Prediction horizon
  Ns <- 24 # [samples] Simulation horizon
  
  # Simulation data
  d <- rep(0.2,max(Np,Ns+Np)) # Disturbance
  p <- rep(0.9,max(Np,Ns+Np)) # Price
  p[12:15] <- 0.2
  
  # Predictive model used in MPC
  H <- ImpulseResponse(A=0,B=-1,E=1,C=1,Np,Ts=1)
  dH <- DifferenceMatrix(Np)
  
  # Closed-loop MPC
  y <- rep(1,Ns+1)
  u <- rep(0,Ns+1)
  for (k in 1:Ns) {
    constraints <- SetConstraints2(Np,H,dH,d=d[k:(k+Np-1)],x0=y[k],um1=u[k],du=2)
    mpc <- BuildMPClp(constraints,Np,p[k:(k+Np-1)])
    opt <- SolveMPClp(mpc)
    u[k+1] <- opt$u[1] # only use first control action
    y[k+1] <- pmax(0,opt$y[2]) # output
    if (opt$status) { warning("Infeasible") }
  }
  
  u <- tail(u,-1) # discard um1
  if (Ns==1) { # openloop when Ns=1
    u <- opt$u
    y <- opt$y
  } 
  
  # Plot results
  plot(y,ylim=c(0,4.1),
       ylab="Normalized signals",xlab="Time [h]",type="l")
  lines(u,type="s",col="blue")
  lines(p,col="red",type="s")
  lines(d,col="brown",type="s")
  lines(rep(constraints$ymin,Ns),lty=2)
  lines(rep(constraints$ymax,Ns),lty=2)
  lines(rep(constraints$umin,Ns),lty=2,col="red")
}
SetConstraints2 <- function(Np,H,dH,d,x0,um1,du=5,ymax=4,ymin=1,umin=0,umax=1) {
  # Setup MPC constraints
  dumin <- rep(-umax/du,Np)
  dumax <- rep(umax/du,Np)
  umin <- rep(umin,Np)
  umax <- rep(umax,Np)  
  ymin <- rep(ymin,Np)
  ymax <- rep(ymax,Np)
  dI0 <- rep(0,Np); dI0[1] <- um1
  
  return(constraints = list(um1=um1,x0=x0,y0=x0,umin=umin,umax=umax,dumin=dumin,dumax=dumax,
                            ymin=ymin,ymax=ymax,dH=dH,dI0=dI0,H=H,d=d))
}

BuildMPClp <- function(con,Np,p) {
  ## Setup MPC problem on standard form
  # minimize f'*x
  # subject to
  #   A*x <= b
  #   x >= 0
  # =
  # minimize( p'*u + r'*s  )
  #         subject to
  #             y = Hu*u + Hx*x0 + Hd*d
  #             ymin - s <= y <= ymax + s
  #             umin <= u <= umax
  #             dumin <= du <= dumax
  #             s>=0, u>=0
  
  # Objective function including slack variable penalty
  f <- as.matrix(c(p,rep(1e6,Np))) # (u,s)
  
  # Constraint matrices
  I <- diag(Np); O <- matrix(0,Np,Np)
  
  A <- rbind(cbind(con$H$u,-I), cbind(-con$H$u,-I),
             cbind(I,O), cbind(-I,O),
             cbind(con$dH,O), cbind(-con$dH,O))
  
  y0d <- con$H$x %*% con$x0 + con$H$d %*% con$d
  b <- as.matrix(c( con$ymax-y0d,
                    -con$ymin+y0d,
                    con$umax,-con$umin,
                    con$dumax+con$dI0,-con$dumin-con$dI0))
  
  if (length(f) != ncol(A) ) {
    stop("MPC objective function (f) length not correct")
  }
  
  return(list(f=f,A=A,b=b,H=con$H,Np=Np,d=con$d,
              x0=con$x0,y0=con$y0,
              umin=con$umin,umax=con$umax,ymin=con$ymin,ymax=con$ymax))
}
SolveMPClp <- function(mpc) {
  # Solve MPC LP problem
  x <- lp("min", mpc$f, mpc$A, rep("<=", nrow(mpc$A)), mpc$b)
  #x <- lp("min", mpc$f, mpc$A, rep("<=", nrow(mpc$A)), mpc$b,binary.vec = 1:mpc$Np)
  
  # Extract optimal solution
  u <- head(x$solution,mpc$Np)
  s <- tail(x$solution,mpc$Np)
  
  # Predict output
  y <- mpc$H$u %*% u + mpc$H$d %*% mpc$d + mpc$H$x*mpc$x0
  y <- c(mpc$y0,y)[1:mpc$Np]
  
  # Check slack variables
  if (sum(s) > 1e-3) { warning("Output constraints violated.") }
  
  return(list(status=x$status,u=u,y=y,s=s))
}

BuildMPCqp <- function(con,Np,r) {
  ## Setup MPC problem on standard form
  # minimize 1/2 x^T HH x + f'*x
  # subject to
  #   A*x <= b
  # 
  # minimize( 1/2 u^T HH u + p'*u )  ||r-y||^2 
  #         subject to
  #             y = Hu*u + Hx*x0 + Hd*d
  #             umin <= u <= umax
  #             dumin <= du <= dumax
  
  # Objective function matrices including slack variable penalty
  Q <- diag(rep(1,Np)) # Weight in L2-objective 
  HH <- t(con$H$u) %*% Q %*% con$H$u # QP Hessian
  
  y0d <- con$H$x%*%con$x0 + con$H$d%*%con$d
  g <- r - y0d
  f <- -t(con$H$u) %*% Q %*% g
  
  # Constraint matrices
  I <- diag(Np)#; O <- matrix(0,Np,Np)
  
  A <- rbind(con$dH, -con$dH,
             I, -I)
  
  b <- as.matrix(c(  con$dumax+con$dI0,
                    -con$dumin-con$dI0,
                     con$umax,
                    -con$umin))
  
  if (length(f) != ncol(A)) {
    stop("MPC objective function (f) length not correct")
  }
  
  return(list(HH=HH,f=f,A=A,b=b,H=con$H,Np=Np,d=con$d,
              x0=con$x0,y0=con$y0,
              umin=con$umin,umax=con$umax))
}
SolveMPCqp <- function(mpc) {
  # Solve MPC QP problem and extract solution
  x <- quadprog(mpc$HH,mpc$f,mpc$A,mpc$b)
  x$status <- is.na(x$solution[1])
  
  u <- head(x$solution,mpc$Np)
  
  y <- mpc$H$u%*%u + mpc$H$d%*%mpc$d + mpc$H$x*mpc$x0 # output prediction
  #y <- c(mpc$y0,y)[1:mpc$Np]

  return(list(status=x$status,u=u,y=y))
}

BuildMPCqps <- function(con,Np,p) {
  ## Setup MPC problem on standard form
  # minimize 1/2 x^T HH x + f'*x
  # subject to
  #   A*x <= b
  # 
  # minimize( 1/2 u^T HH u + p'*u + q'*s  )  ||r-y||^2
  #         subject to
  #             y = Hu*u + Hx*x0 + Hd*d
  #             ymin - s <= y <= ymax + s
  #             umin <= u <= umax
  #             dumin <= du <= dumax
  #             s >= 0
  
  # Objective function matrices including slack variable penalty
  Q <- diag(c(rep(1e-6,Np),rep(1e6,Np))) 
  HH <- t(con$Hu) %*% Q %*% con$Hu # QP Hessian
  
  y0d <- con$H$x%*%con$x0 + con$H$d%*%con$d
  b <- con$r - y0d
  p <- t(con$Hu)
  
  f <- as.matrix(c(g,rep(1e6,Np))) # (p,q) (u,s)
  
  # Constraint matrices
  I <- diag(Np); O <- matrix(0,Np,Np)
  
  A <- rbind(cbind(con$H$u,-I), cbind(-con$H$u,-I),
             cbind(con$dH,O), cbind(-con$dH,O),
             cbind(I,O), cbind(-I,O),
             cbind(O,-I))
  
  
  b <- as.matrix(c( con$ymax-y0d,
                    -con$ymin+y0d,
                    con$dumax+con$dI0,
                    -con$dumin-con$dI0,
                    con$umax,
                    -con$umin,
                    rep(0,Np)))
  
  
  
  if (length(f) != ncol(A)) {
    stop("MPC objective function (f) length not correct")
  }
  
  return(list(f=f,A=A,b=b,H=con$H,Np=Np,d=con$d,
              x0=con$x0,y0=con$y0,
              umin=con$umin,umax=con$umax,ymin=con$ymin,ymax=con$ymax))
}
SolveMPCqps <- function(mpc) {
  # Solve MPC QP problem and extract solution
  print("Solving MPC")
  x <- quadprog(mpc$HH,mpc$f,mpc$A,mpc$b)
  x$status <- is.na(x$solution[1])
  
  u <- head(x$solution,mpc$Np) 
  s <- tail(x$solution,mpc$Np)
  
  y <- mpc$H$u%*%q + mpc$H$d%*%mpc$d + mpc$H$x*mpc$x0 # output prediction
  #y <- c(mpc$y0,y)[1:mpc$Np]
  
  if (sum(s) > 1e-3) { warning("Output constraints violated.") }
  
  return(list(status=x$status,u=u,y=y,s=s))
}

## Generic toolbox functions
# Continous to discrete time models
library(expm) # Matrix exponential for discretization
c2d <- function(A,B,E,Ts,nx=nrow(B),nu=ncol(B),nd=ncol(E)) {
  # Discretize a continuous time state space model:
  # dx/dt = A x + B u + E d, y = C x, (sampling period Ts)
  M  <- rbind(cbind(A,B,E), matrix(0,nu+nd,nx+nu+nd))
  eM <- expm(M*Ts)  
  Ad <- eM[1:nx,1:nx]
  Bd <- eM[1:nx,(nx+1):(nx+nu)]
  Ed <- eM[1:nx,(nx+nu+1):(nx+nu+nd)]
  return(list(A=Ad,B=Bd,E=Ed))
}
ImpulseResponse <- function(A,B,E,C,N,Ts=1,con=T) {
  # Compute impulse response model:
  # Y = Hu*U + Hd*D + Hx*x0
  # from continuous or discrete time state space model
  # dx/dt = A x + B u + E d, y = C x
  # con = T/F: specifies if (A,B,C,E) is continuous or discrete time state space model
  # N: Prediction horizon
  # Ts: Sampling period
  A <- as.matrix(A)
  B <- as.matrix(B)
  E <- as.matrix(E)
  C <- as.matrix(C)
  
  nu <- ncol(B)
  nx <- nrow(B)
  ny <- nrow(C)
  nd <- ncol(E)
  H0 <- matrix(0,N*ny,nx)
  Hu <- matrix(0,N*ny,N*nu)
  Hd <- matrix(0,N*ny,N*nd)
  
  if (con) { # Discretize
    d <- c2d(A,B,E,Ts,nx,nu,nd)
  } else {
    d <- list(A=A,B=B,E=E)
  }
  # Compute first column with all impulse response coefficients
  M <- C # temporary variable
  k1 <- 1
  k2 <- ny
  for (k in 1:N) {
    Hu[k1:k2,1:nu] <- M %*% d$B
    Hd[k1:k2,1:nd] <- M %*% d$E
    M <- M %*% d$A
    H0[k1:k2,1:nx] <- M
    k1 <- k1+ny
    k2 <- k2+ny
  }
  
  # Copy coefficients and fill out remaining columns
  k1row <- ny+1
  k2row <- N*ny
  k1col <- nu+1
  k2col <- nu+nu
  kk <- N*ny-ny
  for (k in 2:N) {
    Hu[k1row:k2row,k1col:k2col] <- Hu[1:kk,1:nu]
    k1row <- k1row+ny
    k1col <- k1col+nu
    k2col <- k2col+nu
    kk <- kk-ny
  }
  
  # Repeat for Hd
  k1row <- ny+1
  k2row <- N*ny
  k1col <- nd+1
  k2col <- nd+nd
  kk <- N*ny-ny
  for (k in 2:N) {
    Hd[k1row:k2row,k1col:k2col] <- Hd[1:kk,1:nd]
    k1row <- k1row+ny
    k1col <- k1col+nd
    k2col <- k2col+nd
    kk <- kk-ny
  }  
  
  return(list(u=Hu,d=Hd,x=H0))
}

# Solvers
library(lpSolve)
library(quadprog)
linprog <- function(f,A,b,G,h) {
  # Solves the LP:
  # minimize f'*x
  # subject to
  #   A*x <= b
  #   G*x = h (optional)
  #   x >= 0
    
  if (nargs() == 5) {
    # Solve with equality constraints (Gx = h)
    x <- lp("min", f, rbind(A,G), c(rep("<=", nrow(A)),rep("=", nrow(G))), c(b,h))    
  } else {
    x <- lp("min", f, A, rep("<=", nrow(A)), b)
  }
  # else if binary.vec / int.vec defined (indices of binary / integer variables )
  
  if (!x$status) { 
    print("Optimization converged :)")
  } else { 
    print("Infeasible optimization problem :(")
  }
  return(x)
}
quadprog <- function(H,f,A,b,G,h) {
  # Solves the QP:
  # minimize 1/2 x'Hx + f'x
  # subject to Ax <= b
  #            Gx = h (optional)
  
  if (nargs() > 4) {
    # Solve with equality constraints (Gx = h)
    x <- solve.QP(Dmat=H,dvec=-f,Amat=-t(rbind(G,A)),bvec=-rbind(h,b),meq=length(h))
  } else {
    x <- solve.QP(Dmat=H,dvec=-f,Amat=-t(A),bvec=-b,meq=0)
  }
  
  return(x)
}

# Support functions
DifferenceMatrix <- function(N=3,nu=1) {
  # Computes difference matrix for rate of change constraints
  Iu <- diag(nu)
  Id <- BelowDiagonal(N) %x% Iu
  dH <- (diag(N+1) %x% Iu) - Id
  return(dH[1:(N*nu),1:(N*nu)])
}
BelowDiagonal <- function(x) {
  if (length(x) == 1) {
    x <- rep(1,x)
  }
  return( rbind( rep(0,length(x)+1), cbind(diag(x),rep(0,length(x))) ) )
}
shuffle <- function(X) {
  # Converts MIMO inputs to single column vector
  return(matrix(t(X),ncol=1))
}
unshuffle <- function(x,n) {
  # Converts single column vector to MIMO inputs
  return(matrix(x,ncol=n,byrow=T))
}
testshuffle <- function() {
  x <- 1:3
  y <- 4:6
  z <- 7:9
  xyz <- c(x,y,z)
  xyz.shuffled <- shuffle(matrix(c(x,y,z), ncol=3))
  xyz.unshuffled <- unshuffle(xyz.shuffled,3)
  tst <- (x == xyz.unshuffled[,1]) & 
    (y == xyz.unshuffled[,2]) & (z == xyz.unshuffled[,3])
  print(tst)
  return(data.frame(xyz,xyz.shuffled,xyz.unshuffled))
}
