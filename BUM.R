# MAP Estimation in the Beta Unfolding Model with censored boundary handling
# Based on: Noel, Y. (2014). A beta unfolding model for continuous bounded responses, Psychometrika, 79(4), 647-674.
# This version adds:
#   1. Gaussian priors on item parameters (delta, lambda, tau) for MAP estimation,
#      regularizing extreme item locations that lie far from the person distribution.
#   2. Censored likelihood for boundary responses (0 and 1): instead of evaluating
#      the Beta density at the boundary (which is infinite or zero), boundary observations
#      contribute P(X <= eps) or P(X >= 1-eps) to the likelihood, treating them as
#      censored events (Tobit-like approach).
require(R6)
require(tensor)

BUM = R6Class("BUM",
  public = list(

    # Data
    N = NULL,
    p = NULL,
    x = NULL,

    # Model parameters
    K = 30, 
    theta = NULL, 
    delta = NULL,
    lambda = NULL, 
    tau = NULL,
    theta.k = NULL, 
    pi.k = NULL, 
    history = NULL,

    # MAP prior hyperparameters (Gaussian priors on item parameters)
    # Set sd to Inf (or very large) to effectively disable a prior
    prior.delta  = list(mean = 0, sd = 3),
    prior.lambda = list(mean = 0, sd = 10),
    prior.tau    = list(mean = -1, sd = 10),

    # Tobit (censored likelihood) approach for boundary responses (0 and 1).
    # When TRUE, observations <= eps are left-censored and >= 1-eps are right-censored.
    # When FALSE (default), boundary values are clipped to [small, 1-small] as in the original BUM.
    tobit = FALSE,
    eps = 0.01,

    # Loglikelihoods
    marginalLogLik = NULL,
    jointLogLik = NULL,

    # Display parameters
    mfrow = c(3, 3), 

    #----------------------------------- Constructor ----------------------------------
    initialize = function (N = 200, p = 20, data = NULL, inverted=NULL, ...) {

      if(!is.null(data)) {
        self$x = as.matrix(data)
        self$N = nrow(data)
        self$p = ncol(data)
      }
      
      else private$Simulate(N,p,...)
      
      # Data within [0;1]
      stopifnot( (min(self$x)>=0) && (max(self$x)<=1) )
      
      # Manage U-shaped ICC if any expected
      if(!is.null(inverted)) {
        self$x[,inverted] = 1 - self$x[,inverted]
        private$inverted = inverted
      }

      private$dsec = array(0,dim=c(self$p,3,3),list(colnames(self$x),c("lambda","delta","tau"),c("lambda","delta","tau")))

    },
    #----------------------------------- Public methods -------------------------------
    Estimate = function (unimodal.constraint = FALSE, fixed = NULL, jml.niter = 2, start.method="CA", prior="gaussian",
                         prior.delta = NULL, prior.lambda = NULL, prior.tau = NULL, tobit = FALSE, bayesian = TRUE, display=TRUE,...)  {

      private$fixed = fixed
      private$unimodal.constraint = unimodal.constraint
      private$prior = prior
      self$tobit = tobit

      # When bayesian=FALSE, use diffuse priors (effectively MLE)
      if(!bayesian) {
        self$prior.delta  = list(mean = 0, sd = 1e6)
        self$prior.lambda = list(mean = 0, sd = 1e6)
        self$prior.tau    = list(mean = 0, sd = 1e6)
      }

      # Update prior hyperparameters if explicitly provided (overrides bayesian flag)
      if(!is.null(prior.delta))  self$prior.delta  = prior.delta
      if(!is.null(prior.lambda)) self$prior.lambda = prior.lambda
      if(!is.null(prior.tau))    self$prior.tau    = prior.tau

      if( (jml.niter>0) && (!unimodal.constraint)) private$JML(fixed=fixed,niter=jml.niter,start.method=start.method,display=display,...)
      private$EM(unimodal.constraint=unimodal.constraint,fixed=fixed,prior=prior,display=display,...)

      self$marginalLogLik = private$logLik
      self$jointLogLik = private$jointLogLikelihood()
    },

    Summary = function (delta.order = TRUE, persons = FALSE, print.it = TRUE) {

      # Standard errors of item parameter estimates
      ase = private$ASE()
          
      # Reorder items by location
      k = 1:length(self$delta)
      if(delta.order) k = sort(self$delta,index.return=TRUE)$ix
      
      # Outfit et infit
      f = self$Fitted()
      eij = self$x - f$mu
      item.outfit = colMeans((eij**2)/f$sigma2)[k]
      item.infit  = colSums(eij**2)[k]/colSums(f$sigma2)[k]

      fit1 = cbind(DELTA=self$delta[k])
      rownames(fit1) = colnames(self$x)[k]

      if("delta" %in% private$fixed)  fit1 = cbind(fit1,A.S.E=0,LAMBDA=self$lambda[k])
      else                            fit1 = cbind(fit1,A.S.E=ase[k,"delta"],LAMBDA=self$lambda[k])

      if("lambda" %in% private$fixed) fit1 = cbind(fit1,A.S.E=0,TAU=self$tau[k])
      else                            fit1 = cbind(fit1,A.S.E=ase[k,"lambda"],TAU=self$tau[k])

      if("tau" %in% private$fixed)    fit1 = cbind(fit1,A.S.E=0)
      else                            fit1 = cbind(fit1,A.S.E=ase[k,"tau"])
      
      fit1 = cbind(fit1,OUTFIT=item.outfit)
      fit1 = cbind(fit1,INFIT=item.infit)
      
      fit1 = as.data.frame(fit1)
      bad = which( (fit1[,"INFIT"] < 0.7) | (fit1[,"INFIT"] > 1.3) | (fit1[,"OUTFIT"] < 0.5) | (fit1[,"OUTFIT"] > 1.5) )
      fit1$Misfit = ""
      fit1$Misfit[bad] = "***"

      # Person parameters and person fit statistics
      fit2 = NULL
      if(persons) {
        postSD = private$thetaPosteriorSD()
        fit2 = cbind(THETA=self$theta,postSD = postSD)
        if("theta" %in% private$fixed) fit2[,"postSD"] = 0
        
        if(print.it) {
          cat("\nPERSON PARAMETERS AND FIT STATISTICS\n\n")
          print(round(fit2,3))
        }
      }

      if(print.it) {
        cat("\nITEM PARAMETERS AND STATISTICS\n\n")
        # print(round(fit1,3))
        print(fit1)
        edf = self$effectiveDF()
        cat("\nEffective df =",round(edf,2),"(nominal =",private$npar,")\n")
        cat("AIC (effective) =",round(self$Aic(effective=TRUE),2),"\n")
        cat("AIC (nominal)   =",round(self$Aic(effective=FALSE),2),"\n")
      }

      invisible(list(items=fit1,persons=fit2))
    },

    # -------- Fonction principale --------
    toLatex = function(df, caption="Titre de la table.", label="tab:tableau", digits=2, symmetric=FALSE) {

      stopifnot(is.data.frame(df))

      # Charger les noms de lignes si présents
      rn <- rownames(df)
      if (!is.null(rn)) {
        df <- cbind(Row = rn, df)
      }

      # Arrondir les valeurs numériques
      df_fmt <- df
      for (j in seq_along(df_fmt)) {
        if (is.numeric(df_fmt[[j]])) {
          df_fmt[[j]] <- formatC(df_fmt[[j]], format = "f", digits = digits)
        }
      }

      # For symmetric matrices (e.g. correlations), blank out the upper triangle
      # (must be done after rounding, to avoid coercing numeric columns to character too early)
      if (symmetric) {
        for (i in seq_len(nrow(df_fmt))) {
          for (j in seq_along(df_fmt)[-1]) {  # skip the Row column
            if (j - 1 > i) df_fmt[i, j] <- ""
          }
        }
      }
      
      # Colonnes : l pour la première, c pour les autres
      colspec <- paste0("l", paste(rep("c", ncol(df_fmt) - 1), collapse = ""))
      
      # En-têtes
      headers <- paste(names(df_fmt), collapse = " & ")
      
      # Lignes
      body <- apply(df_fmt, 1, function(x) paste(x, collapse = " & "))
      body_text <- paste(body, collapse = " \\\\\n")
      
      # Génération du code LaTeX
      out <- paste0(
      "\\begin{table}[tbp]
      \\centering
      \\caption{", caption, "}
      \\label{", label, "}
      \\begin{tabular}{", colspec, "}
      \\toprule
      ", headers, " \\\\
      \\midrule
      ", body_text, " \\\\
      \\bottomrule
      \\end{tabular}
      \\end{table}"
      )
      
      return(out)
    },

    Aic = function (joint=TRUE, effective=TRUE) {

      if(effective) np = self$effectiveDF()
      else          np = private$npar

      if(joint) AIC = -2*self$jointLogLik + 2*np
      else      AIC = -2*private$logLik[length(private$logLik)] + 2*np

      AIC
    },

    # Effective degrees of freedom under MAP estimation.
    # For each item, computes p_eff_j = trace(H_post^{-1} H_lik) where
    #   H_post = Hessian of negative log-posterior (data + prior)
    #   H_lik  = H_post - P  (data only, removing the prior precision)
    #   P      = diag(1/sd_lambda^2, 1/sd_delta^2, 1/sd_tau^2)
    # Returns the total effective number of parameters (items + persons/latent dist).
    effectiveDF = function (verbose = FALSE) {

      P.prior = diag(c(1/self$prior.lambda$sd^2,
                        1/self$prior.delta$sd^2,
                        1/self$prior.tau$sd^2))

      # Choose the right objective/gradient depending on estimation method
      if(private$method == "JML") {
        objfn  = private$itemLogLikFunction
        grfn   = private$itemGradient
      } else {
        objfn  = private$itemEMFunction
        grfn   = private$itemEMGradient
      }

      item.edf = numeric(self$p)
      for(j in 1:self$p) {
        params = c(lambda=unname(self$lambda[j]), delta=unname(self$delta[j]), tau=unname(self$tau[j]))

        # Numerical Hessian of the penalized (MAP) objective at convergence
        H.post = optimHess(params, objfn, grfn, item=j)

        # Effective df: 3 - trace(H_post^{-1} P)
        # = number of free params minus the shrinkage from priors
        edf = tryCatch({
          3 - sum(diag(solve(H.post) %*% P.prior))
        }, error = function(e) 3)  # fallback to nominal if Hessian singular

        item.edf[j] = max(0, min(3, edf))  # clamp to [0, 3]
      }

      if(verbose) {
        cat("Effective df per item:\n")
        edf.df = data.frame(Item=colnames(self$x), edf=round(item.edf, 3))
        print(edf.df, row.names=FALSE)
        cat("Total item edf:", round(sum(item.edf), 2), "out of nominal", 3*self$p, "\n")
      }

      # Total: item effective df + person/latent distribution parameters
      total = sum(item.edf)

      if(private$method == "JML") {
        total = total + self$N - 1  # N person params - 1 location constraint
        if("theta" %in% private$fixed) total = total - (self$N - 1)
      } else {
        if(private$prior != "gaussian") total = total + (self$K - 1) - 1
        # If Gaussian prior: K-1 class probs are fixed, only item params contribute
      }

      # Corrections for fixed item parameters (they have edf=0 already from the
      # Hessian, but if they were truly fixed, optimHess wasn't called with the
      # right bounds. Override to be safe.)
      # Not needed: if a parameter is fixed, H_post for that dimension → ∞,
      # so solve(H_post) → 0, and edf → 3 - trace → correct.
      # But the bounds in optim make lb=ub, so optimHess may fail. Already handled by tryCatch.

      total
    },

    Predict = function (tt = NULL) {

      if(is.null(tt)) tt = seq(min(self$theta),max(self$theta),length=100)
      N = length(tt)

      T = outer(rep(1,N),self$tau,"*")
      D = outer(tt,self$delta,"-")
      L = outer(rep(1,N),self$lambda,"*")

      m = exp(L)                          # Acceptation
      n = exp(D+T)+exp(-D+T)              # Refusal

      mu = m / (m+n)                      # Expectation
      colnames(mu) = colnames(self$x)
      md = (m-1) / (m+n-2)                # Mode
      colnames(md) = colnames(self$x)
      md[md>1] = 1
      md[md<0] = 0
      med = qbeta(.50,m,n)                # Median
      colnames(med) = colnames(self$x)
      q5 = qbeta(.05,m,n)                 # 5th percentile
      colnames(q5) = colnames(self$x)
      q95 = qbeta(.95,m,n)                # 95th percentile
      colnames(q95) = colnames(self$x)
      sigma2 = mu * (1-mu) / (m + n + 1)  # Variance
      colnames(mu) = colnames(self$x)

      list(theta=tt,m=m,n=n,mu=mu,mode=md,median=med,q5=q5,q95=q95,sigma2=sigma2)
    },

    Fitted = function() {

      T = outer(rep(1,self$N),self$tau,"*")
      D = outer(self$theta,self$delta,"-")
      L = outer(rep(1,self$N),self$lambda,"*")

      m = exp(L)
      n = exp(D+T)+exp(-D+T)

      mu = m / (m+n)
      sigma2 = mu * (1-mu) / (m + n + 1)

      list(mu=mu,sigma2=sigma2)
    },

    Information = function(tt = NULL) {

      if(is.null(tt)) tt = seq(min(self$theta),max(self$theta),length=100)
      N = length(tt)
      
      T = outer(rep(1,N),self$tau,"*")
      D = outer(tt,self$delta,"-")
      L = outer(rep(1,N),self$lambda,"*")

      m = exp(L)
      n = exp(D+T)+exp(-D+T)

      -4*((exp(T)*sinh(D))**2) * (trigamma(m+n)-trigamma(n))

    },  

    Plot = function (type = "erf", items = NULL, delta.order = TRUE, plot.data = TRUE, plot.exp=TRUE, lwd=3,
        plot.modal = FALSE, plot.median = FALSE, plot.ci = FALSE, conf.level = 0.9,
        plot.locations = FALSE, plot.smooth = FALSE, gam.k = 4, nlevels = 10, true.par = NULL,
        param1 = "delta", param2 = "lambda", theta.angle = -30, phi.angle = 25,
        f = 2/3, xlim, ylim, main = NULL, mfrow = NULL, mar=c(4,4.5,3,2), pch=19, cex=.4, col.point=gray(0.3,0.4), col.exp="#4285f4", col.modal="#0f9d58", col.median="#db4437", col.ci="lightgrey", col.smooth="#f4b400",col.locations="#333333",...) {


      if(type=="erf") {

        model = self$Predict()
        if(is.null(items)) items = 1:self$p
        if(!is.numeric(items)) items = match(items,colnames(self$x)) # convert column names into indices

        k = items
        if(delta.order) { k = sort(self$delta[items],index.return=TRUE)$ix ; items = items[k] }
        titles = if(is.null(main)) colnames(self$x)[items] else main

        if(!is.null(mfrow)) par(mfrow=mfrow)

        panel = 0

        for(i in items) {
          panel = panel + 1
          if(!is.null(mar)) par(mar=mar)

          x = self$x[,i]
          expected = model$mu[,i]
          modal = model$mode[,i]
          medn = model$median[,i]
          q5 = model$q5[,i]
          q95 = model$q95[,i]

          # Manage inverted items
          if(colnames(self$x)[i] %in% private$inverted) {
            x = 1-x
            expected = 1-expected
            modal = 1-modal
            medn = 1-medn
            q5 = 1-q5
            q95 = 1-q95
          }

          plot(self$theta,x,type="n",xlab=expression(hat(theta)[i]),ylab=expression(x[ij]),main=titles[panel],ylim=c(0,1),...)
          if(plot.ci) {
            x.ci = c(model$theta,rev(model$theta),model$theta[1])
            y.ci = c(q5,rev(q95),q5[1])
            polygon(x.ci,y.ci,col=gray(0.9),border=NA)
          }
          if(plot.data)      points(self$theta,x,pch=pch,cex=cex,col=col.point,...)
          if(plot.exp)       lines(model$theta,expected,lwd=lwd,col=col.exp,...)
          if(plot.modal)     lines(model$theta,modal,lwd=lwd,col=col.modal,...)
          if(plot.median)    lines(model$theta,medn,lwd=lwd,col=col.median,...)
          if(plot.locations) lines(rbind(c(self$delta[i],0),c(self$delta[i],1)),col=col.locations,...)
          if(plot.smooth) {
            # lines(lowess(self$theta,x,f=f),col=col.smooth,lty=2,lwd=2)
            require(mgcv)
            tt = self$theta
            m.gam = gam(x~s(tt,k=gam.k,fx=TRUE),family="betar")
            tt2 = seq(min(self$theta),max(self$theta),len=100)
            gam.fit = predict(m.gam,newdata=data.frame(tt=tt2),type="response")
            lines(tt2,gam.fit,col=col.smooth,lty=2,lwd=2)
          }
        }
        
        # par(opar)
      }

      # Draw item response functions from a simulated model
      else if(type=="true.model") {

        if(is.null(private$true.par)) return()
        theta =  seq(min(private$true.par$theta),max(private$true.par$theta),len=200)
        model =  private$Predict2(tt=theta,conf.level=conf.level)
        
        delta =  private$true.par$delta
        lambda = private$true.par$lambda
        tau =    private$true.par$tau
        
        if(is.null(items)) items = 1:self$p
        
        if(!is.null(mfrow)) par(mfrow=mfrow)

        k = items
        if(delta.order) { k = sort(delta[items],index.return=TRUE)$ix ; items = items[k] }
        titles = if(is.null(main)) colnames(self$x)[items] else main
        true.theta = private$true.par$theta
        panel = 0

        for(i in items) {
          panel = panel + 1
          if(!is.null(mar)) par(mar=mar,cex.main=1.5,cex.lab=1.5)
          x = self$x[,i]
          expected = model$mu[,i]
          modal = model$mode[,i]
          medn = model$median[,i]
          q5 = model$q5[,i]
          q95 = model$q95[,i]

          # Manage inverted items
          if(colnames(self$x)[i] %in% private$inverted) {
            x = 1-x
            expected = 1-expected
            modal = 1-modal
            medn = 1-medn
            q5 = 1-q5
            q95 = 1-q95
          }

          plot(true.theta,x,type="n",xlab=expression(hat(theta)[i]),ylab=expression(x[ij]),main=titles[panel],ylim=c(0,1),...)
          if(plot.ci) {
            x.ci = c(theta,rev(theta),theta[1])
            y.ci = c(q5,rev(q95),q5[1])
            polygon(x.ci,y.ci,col=col.ci,border=NA)
          }
          if(plot.data)      points(true.theta,x,pch=pch,cex=cex,col=col.point,...)
          if(plot.exp)       lines(theta,expected,lwd=lwd,col=col.exp,...)
          if(plot.modal)     lines(theta,modal,lwd=2,col=col.modal,...)
          if(plot.median)    lines(theta,medn,lwd=2,col=col.median,...)
          if(plot.locations) lines(rbind(c(delta[i],0),c(delta[i],1)),col=col.locations,...)
        }
        
        # par(opar)
      
      }

      else if(type=="theta.dist") {
      
        opar = par(mfrow=c(1,1))
        hist(self$theta,freq=FALSE,xlab=expression(hat(theta)),ylab="Density",main="Person parameter distribution")
        par(opar)
      }
      
      else if(type=="latent.dist") {
      
        opar = par(mfrow=c(1,1))
        plot(self$theta.k,self$pi.k,xlab=expression(theta[k]),ylab=expression(pi[k]),main="Latent distribution",type="h",lwd=2)
        par(opar)
      }
      
      else if (type=="convergence") {
        par(mfrow=c(1,1))
        plot(private$target,xlab="Iterations",ylab="loglikelihood",main="Gain",type="l",lwd=2)
      }
      
      else if (type=="simulation") {

        if(is.null(private$true.par)) return

        opar = par(mfrow=c(3,2))
        
        plot(private$true.par$theta,self$theta,xlab=expression(theta[i]),ylab=expression(hat(theta[i])),main="Attitudes")
        abline(0,1,col="red")

        plot(private$true.par$delta,self$delta,xlab=expression(delta[j]),ylab=expression(hat(delta[j])),main="Item locations",type="n")
        abline(0,1,col="red")
        text(private$true.par$delta,self$delta,paste(1:self$p),cex=.6)
        
        plot(private$true.par$lambda,self$lambda,xlab=expression(lambda[j]),ylab=expression(hat(lambda[j])),main="Acceptances",type="n")
        abline(0,1,col="red")
        text(private$true.par$lambda,self$lambda,paste(1:self$p),cex=.6)
        
        plot(private$true.par$tau,self$tau,xlab=expression(tau[j]),ylab=expression(hat(tau[j])),
                                  main="Dispersions",type="n")
        abline(0,1,col="red")
        text(private$true.par$tau,self$tau,paste(1:self$p),cex=.6)
        
        if(!is.null(private$logLik)) plot(private$logLik,xlab="Iterations",ylab="Gain",main="Log Likelihood",type="l",lwd=2)
        
        par(opar)

      }

      else if (type=="3D") {
        require(rgl)
        model = self$Predict()
        pc = princomp(self$x)
        pred = predict(pc,model$mu)
        rgl.bg(color="white")
        rgl.clear()
        rgl.points(pc$scores[,1:3],color="blue",size=4)
        rgl.lines(pred[,1:3],size=3,color="red")
        rgl.bbox()
      }
      
      else if(type=="contour") {

        tt = seq(min(self$theta),max(self$theta),len=100)
        model = self$Predict(tt)
        if(is.null(items)) items = 1:self$p
        
        opar = par(mfrow=mfrow,mar=c(4,4.5,3,2),cex.main=1.5,cex.lab=1.5)
        k = items
        if(delta.order) { k = sort(self$delta[items],index.return=TRUE)$ix ; items = items[k] }
        titles = if(is.null(main)) colnames(self$x)[items] else main

        y = seq(0,1,len=50)
        panel = 0

        for(i in items) {
          panel = panel + 1

          # TODO: Manage item inversion
          x = self$x[,i]
          expected = model$mu[,i]

          zz = NULL
          for(score in (1:length(tt))) {
            zz = rbind(zz,dbeta(y,model$m[score,i],model$n[score,i]))
          }
          plot(self$theta,x,main=titles[panel],xlab=expression(hat(theta)[i]),ylab=expression(x[ij]),col="grey50",pch=19,cex=.5)
          lines(model$theta,expected,lwd=3,col="grey30")
          contour(tt,y,zz,add=TRUE,nlevels=nlevels)
        }
        par(opar)
      }
      
      else if(type=="filled.contour") {

        # TODO: Manage item inversion

        tt = seq(min(self$theta),max(self$theta),len=100)
        model = self$Predict(tt)
        if(is.null(items)) items = 1
        
        y = seq(0,1,len=50)
        i = items[1]   # One plot only: Filled.contour does not allow for multiple panels
        zz = NULL
        for(score in (1:length(tt))) { 
          zz = rbind(zz,dbeta(y,model$m[score,i],model$n[score,i]))
        }
        
        zz[!is.finite(zz)] = -1
        zz[zz == -1] = max(zz)

        rgb.palette = colorRampPalette(c("white", "blue", "black"),space = "rgb")
        pal = rgb.palette(10)

        filled.contour(tt,y,zz,nlevels=nlevels,color.palette=rgb.palette,main=colnames(self$x)[i],xlab=expression(hat(theta)[i]),ylab=expression(x[ij]),
                      plot.axes = { 
                        points(self$theta,x,col="grey50",pch=19,cex=.5) ;
                        lines(model$theta,expected,lwd=3,col="grey30") ; 
                        if(plot.smooth)    lines(lowess(self$theta,x,f=f),col="darkgoldenrod2",lty=2,lwd=2)
                        axis(1) ; axis (2) } )

      }
      
      else if(type=="surface2") {

        tt = seq(min(self$theta),max(self$theta),len=30)
        model = self$Predict(tt)
        if(is.null(items)) items = 1:self$p
        
        opar = par(mfrow=mfrow,mar=c(1,1,2.5,1),cex.main=1.5,cex.lab=1.5)
        k = items
        if(delta.order) { k = sort(self$delta[items],index.return=TRUE)$ix ; items = items[k] }
        
        y = seq(0,1,len=32)[2:31]
        
        for(i in items) {
          zz = NULL
          for(score in (1:length(tt))) { 
            zz = rbind(zz,dbeta(y,model$m[score,i],model$n[score,i]))
          }
          rgb.palette = colorRampPalette(c("white","white","blue"),space = "rgb")

          fcol = matrix("white", nr = length(tt)-1, nc = ncol(zz)-1) 
          zi = (zz[-1,-1] + zz[-1,-ncol(zz)] + zz[-length(tt),-1] + zz[-length(tt),-ncol(zz)])/4
          fcol = rgb.palette(20)[cut(zi, quantile(zi, seq(0,1, len = 21)),include.lowest = TRUE)]
          zz[!is.finite(zz)] = NA
          persp(tt,y,zz,col=fcol,xlab="Attitude",ylab="Response",zlab="Density",theta=-50,phi=55,ticktype="detailed",r=2.5,main=colnames(self$x)[i],...)
        }
        par(opar)
      }
      
      else if(type=="surface") {

        tt = seq(min(self$theta),max(self$theta),len=30)
        model = self$Predict(tt)
        if(is.null(items)) items = 1:self$p
        
        opar = par(mfrow=mfrow,mar=c(1,1,2.5,1),cex.main=1.5,cex.lab=1.5)
        k = items
        if(delta.order) { k = sort(self$delta[items],index.return=TRUE)$ix ; items = items[k] }
        
        y = seq(0,1,len=32)[2:31]
        
        for(i in items) {
          zz = NULL
          for(score in (1:length(tt))) { 
            zz = rbind(zz,dbeta(y,model$m[score,i],model$n[score,i]))
          }
          zz[!is.finite(zz)] = NA
          persp(tt,y,zz,xlab="Attitude",ylab="Response",zlab="Density",theta=-50,phi=55,ticktype="detailed",shade=.75,r=2.5,main=colnames(self$x)[i],...)
        }
        par(opar)
      }
      
      else if(type=="mosaique") {
        require(MASS)
        tt = seq(min(self$theta),max(self$theta),len=100)
        model = self$Predict(tt)
        if(is.null(items)) items = 1:self$p
        
        opar = par(mfrow=mfrow,mar=c(4,4.5,3,2),cex.main=1.5,cex.lab=1.5)
        k = items
        if(delta.order) { k = sort(self$delta[items],index.return=TRUE)$ix ; items = items[k] }
        
        y = seq(0,1,len=50)
        
        for(i in items) {
          zz = NULL
          for(score in (1:length(tt))) { 
            zz = rbind(zz,dbeta(y,model$m[score,i],model$n[score,i]))
          }
          f1 = kde2d(self$theta,self$x[,i])
          image(f1,col=rev(gray.colors(30)),main=colnames(self$x)[i])
          contour(tt,y,zz,add=TRUE,col="red")
        }
        par(opar)
      }
      
      else if(type=="item.information") {
        if(is.null(items)) items = 1:self$p
        
        opar = par(mfrow=mfrow,mar=c(4,4.5,3,2),cex.main=1.5,cex.lab=1.5)
        k = items
        if(delta.order) { k = sort(self$delta[items],index.return=TRUE)$ix ; items = items[k] }
        
        tt = seq(min(self$theta),max(self$theta),length=100)
        I = self$Information(tt)

        for(i in items) {
          plot(tt,I[,i],main=colnames(self$x)[i],xlab=expression(theta),ylab="Information",ylim=c(min(I),max(I)),type="l",lwd=2)
        }
        
        par(opar)
      }
      
      else if(type=="test.information") {
        
        tt = seq(min(self$theta),max(self$theta),length=100)
        I = rowSums(self$Information(tt))

        plot(tt,I,main="Test information function",ylim=c(0,max(I)),xlab=expression(theta),ylab="Fisher information",type="l",lwd=2)
      }
      
      else if(type=="exponents") {

        tt = seq(min(self$theta),max(self$theta),len=100)
        model = self$Predict(tt)
        if(is.null(items)) items = 1:self$p
        
        opar = par(mfrow=mfrow,mar=c(4,4.5,3,2),cex.main=1.5,cex.lab=1.5)
        k = items
        if(delta.order) { k = sort(self$delta[items],index.return=TRUE)$ix ; items = items[k] }
        
        for(i in items) {
          plot(tt,model$n[,i],main=colnames(self$x)[i],xlab=expression(hat(theta)[i]),ylab="Exponents",type="l",lwd=2,ylim=c(0,10),col="red")
          lines(tt,model$m[,i],col="blue",lwd=2)
          abline(h=1,lty=2)
        }
        par(opar)
      }
      
      else if(type=="logl.surface") {
      
        opar = par(mfrow=mfrow,mar=c(2,2.5,2,2))
        rgb.palette = colorRampPalette(c("white","white","black"),space = "rgb")
        nticks = 21

        for(i in items) {
          
          model = list(delta=self$delta[i],lambda=self$lambda[i],tau=self$tau[i],x=self$x[,i])
          sim = model
          logl1 = matrix(0,nticks,nticks)
          
          if(missing(xlim)) param1.range = model[[param1]] + seq(-2,2,len=nticks)
          else param1.range = seq(xlim[1],xlim[2],len=nticks)
          if(missing(ylim)) param2.range = model[[param2]] + seq(-2,2,len=nticks)
          else param2.range = seq(ylim[1],ylim[2],len=nticks)
          dimnames(logl1)=list(param1.range,param2.range)
          
          # Compute surface
          for(i2 in 1:nticks) {
            for(j2 in 1:nticks) {
              sim[[param1]] = param1.range[i2]
              sim[[param2]] = param2.range[j2]
              logl1[i2,j2] = private$itemJointLogLikelihoodFunction(sim)
            }
          }
          
          logl1 = logl1/1000 # Rescale to avoid text superposition
          fcol = matrix("white", nr = nticks-1, nc = nticks-1) 
          zi = (logl1[-1,-1] + logl1[-1,-nticks] + logl1[-nticks,-1] + logl1[-nticks,-nticks])/4
          fcol = rgb.palette(20)[cut(zi,quantile(zi,seq(0,1, len = 21),na.rm=TRUE),include.lowest = TRUE)]
          
          # Plot loglikelihood surface in the (param1,param2) subspace
          p = persp(param1.range,param2.range,logl1,col=fcol,xlab=param1,ylab=param2,zlab="Loglikelihood (x 1000)",theta=theta.angle,phi=phi.angle,ticktype="detailed",r=2.5,main=eval(substitute(expression(paste("Item ",i,": ",delta==a," ",lambda==b," ",tau==c)), list(i=i,a=round(model$delta,1),b=round(model$lambda,1),c=round(model$tau,1)))))
          
          # Plot contour lines at the bottom
          ct = contourLines(param1.range,param2.range,logl1,levels=seq(min(logl1),0,len=50))
          for(i2 in 1:length(ct)) lines(trans3d(ct[[i2]]$x,ct[[i2]]$y,min(logl1),pmat=p))

          # Plot point of true parameter values, if available
          if(!is.null(private$true.par)) points(trans3d(private$true.par[[param1]][i],private$true.par[[param2]][i],min(logl1),pmat=p),pch=4,col="blue")
          
          # Plot point of estimated parameter values
          points(trans3d(model[[param1]],model[[param2]],max(logl1),pmat=p),pch=4,col="red")
          lines(trans3d(c(model[[param1]],model[[param1]]),c(model[[param2]],model[[param2]]),c(min(logl1),max(logl1)),pmat=p),lty=2,col="red")
          lines(trans3d(c(model[[param1]],model[[param1]]),c(model[[param2]],param2.range[1]),c(min(logl1),min(logl1)),pmat=p),lty=2,col="red")
          lines(trans3d(c(model[[param1]],param1.range[1]),c(model[[param2]],model[[param2]]),c(min(logl1),min(logl1)),pmat=p),lty=2,col="red")
          
          # plot estimation trajectory in the parameter space
          if(!is.null(self$history))
          lines(trans3d(self$history[,paste(param1,i,sep=".")],self$history[,paste(param2,i,sep=".")],rep(min(logl1),nrow(self$history)),pmat=p),col="darkgreen")
        }
        par(opar)
      }

    }
  ),
  private = list(
    #----------------- Private methods

    # Log-likelihood contributions for a single item.
    # When tobit=TRUE: censored likelihood for boundary observations.
    # When tobit=FALSE: standard dbeta (data already clipped to [small, 1-small]).
    # Returns a list with:
    #   $ll  : vector of length N with per-observation log-likelihood contributions
    #   $left, $right, $interior : logical vectors marking censoring status
    censoredLogLik = function(x, m, n) {
      N = length(x)
      ll = numeric(N)

      if(self$tobit) {
        left  = (x <= self$eps)
        right = (x >= 1 - self$eps)
        interior = !(left | right)

        if(any(interior)) ll[interior] = dbeta(x[interior], m[interior], n[interior], log=TRUE)
        if(any(left))     ll[left]     = pbeta(self$eps,     m[left],     n[left],     log.p=TRUE)
        if(any(right))    ll[right]    = pbeta(1 - self$eps, m[right],    n[right],    lower.tail=FALSE, log.p=TRUE)
      } else {
        left  = rep(FALSE, N)
        right = rep(FALSE, N)
        interior = rep(TRUE, N)
        ll = dbeta(x, m, n, log=TRUE)
      }

      list(ll=ll, left=left, right=right, interior=interior)
    },

    # Gradient weights for censored observations.
    # For interior obs: d/dparam log f(x|m,n) uses the standard dbeta gradient.
    # For left-censored:  d/dparam log F(eps|m,n) = f(eps|m,n) / F(eps|m,n) * d/dparam F(eps|m,n)/f(eps|m,n)
    #   but more conveniently, we need d log P(X<=eps) / d(n) = d/dn log pbeta(eps,m,n).
    # We compute the "hazard ratio" h = dbeta(boundary, m, n) / pbeta(boundary, m, n)
    # and then chain-rule through the Beta density gradients w.r.t. m, n.
    #
    # d/dm log pbeta(eps,m,n) = [d/dm log Beta_density] * dbeta(eps,m,n) / F(eps,m,n)
    #   ... but this is NOT exact because pbeta depends on m,n through the integral.
    #
    # Correct approach: use numerical differentiation of log pbeta w.r.t. (m,n).
    # However, this is expensive. Instead we use the identity:
    #   d/dm log P(X<=eps) = E[log(X) | X<=eps] - (digamma(m) - digamma(m+n))  ... NO, this is not right.
    #
    # Actually, the simplest correct approach for the censored Beta:
    #   log L_censored = log pbeta(eps, m, n)  [for left-censored]
    # We need d/d(param) log pbeta(eps, m, n) where param affects m and n.
    #
    # Using the chain rule and the fact that d/da I_x(a,b) = I_x(a,b) * [log(x) - psi(a) + psi(a+b)]
    # ... this is NOT the correct derivative of the regularized incomplete beta.
    #
    # The correct derivatives are (see Boik & Robison-Cox, 1998):
    #   d/da log I_x(a,b) = [d/da I_x(a,b)] / I_x(a,b)
    # where d/da I_x(a,b) can be computed via numerical differentiation.
    #
    # For computational tractability, we use a small finite difference.
    # This is computed in the gradient functions themselves.

    itemLogLikFunction = function (params, item) {

      DA = (self$theta - params['delta'])

      m = exp(params['lambda'])
      n = exp(DA+params['tau']) + exp(-DA+params['tau'])

      # Censored log-likelihood
      cl = private$censoredLogLik(self$x[,item], rep(m, self$N), n)

      # Negative log-posterior = negative log-likelihood + negative log-prior
      nll = -sum(cl$ll)
      nll = nll - dnorm(params['delta'],  self$prior.delta$mean,  self$prior.delta$sd,  log=TRUE)
      nll = nll - dnorm(params['lambda'], self$prior.lambda$mean, self$prior.lambda$sd, log=TRUE)
      nll = nll - dnorm(params['tau'],    self$prior.tau$mean,    self$prior.tau$sd,    log=TRUE)
      nll

    },

    itemGradient = function (params, item) {

      D = self$theta - params['delta']
      DA = D

      m = exp(params['lambda'])
      n = exp(DA+params['tau']) + exp(-DA+params['tau'])
      x = self$x[,item]

      if(self$tobit) {
        left  = (x <= self$eps)
        right = (x >= 1 - self$eps)
      } else {
        left  = rep(FALSE, self$N)
        right = rep(FALSE, self$N)
      }
      interior = !(left | right)

      # --- Interior observations: standard Beta gradient ---
      g.lambda = 0; g.delta = 0; g.tau = 0

      if(any(interior)) {
        g.lambda = g.lambda - m*sum( (digamma(m+n[interior])-digamma(m)+log(x[interior])) )
        g.delta  = g.delta  + 2*exp(params['tau'])*sum( sinh(DA[interior])*(digamma(m+n[interior])-digamma(n[interior])+log(1-x[interior])) )
        g.tau    = g.tau    - sum( n[interior]*(digamma(m+n[interior])-digamma(n[interior])+log(1-x[interior])) )
      }

      # --- Censored observations: finite-difference gradient of log pbeta ---
      if(self$tobit) {
        h = 1e-5  # step size for finite differences

        if(any(left)) {
          lp0 = pbeta(self$eps, m, n[left], log.p=TRUE)
          m1 = exp(params['lambda'] + h)
          lp1 = pbeta(self$eps, m1, n[left], log.p=TRUE)
          g.lambda = g.lambda - sum((lp1 - lp0) / h)

          n1 = exp(DA[left]+h+params['tau']) + exp(-(DA[left]+h)+params['tau'])
          lp1 = pbeta(self$eps, m, n1, log.p=TRUE)
          g.delta = g.delta + sum((lp1 - lp0) / h)

          n1 = exp(DA[left]+params['tau']+h) + exp(-DA[left]+params['tau']+h)
          lp1 = pbeta(self$eps, m, n1, log.p=TRUE)
          g.tau = g.tau - sum((lp1 - lp0) / h)
        }

        if(any(right)) {
          lp0 = pbeta(1-self$eps, m, n[right], lower.tail=FALSE, log.p=TRUE)
          m1 = exp(params['lambda'] + h)
          lp1 = pbeta(1-self$eps, m1, n[right], lower.tail=FALSE, log.p=TRUE)
          g.lambda = g.lambda - sum((lp1 - lp0) / h)

          n1 = exp(DA[right]+h+params['tau']) + exp(-(DA[right]+h)+params['tau'])
          lp1 = pbeta(1-self$eps, m, n1, lower.tail=FALSE, log.p=TRUE)
          g.delta = g.delta + sum((lp1 - lp0) / h)

          n1 = exp(DA[right]+params['tau']+h) + exp(-DA[right]+params['tau']+h)
          lp1 = pbeta(1-self$eps, m, n1, lower.tail=FALSE, log.p=TRUE)
          g.tau = g.tau - sum((lp1 - lp0) / h)
        }
      }

      # Add gradient of negative log-prior: d/dx [-log N(x|mu,sd)] = (x - mu) / sd^2
      g.lambda = g.lambda + (params['lambda'] - self$prior.lambda$mean) / self$prior.lambda$sd^2
      g.delta  = g.delta  + (params['delta']  - self$prior.delta$mean)  / self$prior.delta$sd^2
      g.tau    = g.tau    + (params['tau']    - self$prior.tau$mean)    / self$prior.tau$sd^2

      c(lambda = g.lambda, delta = g.delta, tau = g.tau)

    },

    # Estimation by Joint Maximum Likelihood (just for a preliminary configuration)
    JML = function (niter = 200, lower = c(-5, -10, -5), upper = c(5, 10, 5), plot.type = "erf",
        fixed = NULL, unimodal.constraint = FALSE,
        save.plot = FALSE, save.history = TRUE, save.path = getwd(), display = TRUE,  
        verbose = TRUE, mfrow = c(3, 3), start.method = "CA", ...) {

      # To constrain the response density, you have to fit the model uncontrained first, to detect BRD items
      if(unimodal.constraint && !start.config) {
        cat("Error: You have to fit the model uncontrained first, to detect BRD items.\n")
        return()
      }
      
      private$start.time = proc.time()
      private$method = "JML"
      private$logLik = NULL
      self$mfrow = mfrow
      private$fixed=fixed

      # Clip boundary values unless tobit is enabled
      if(!self$tobit) {
        for(j in 1:self$p) {
          self$x[self$x[,j] <    private$small,j] = private$small
          self$x[self$x[,j] > (1-private$small),j] = 1 - private$small
        }
      }

      # Starting configuration
      if(!private$start.config)  private$CA(verbose=verbose)

      paramsTable = cbind(lambda=self$lambda,delta=self$delta,tau=self$tau)
      control = list(maxit=10)

      # Save history
      if(save.history) {
        self$history = rbind(self$history,matrix(c(0,0,self$delta,self$lambda,self$tau), 1, 3*self$p +2))
        colnames(self$history) = c("Algo","Iter",paste("delta",1:self$p,sep="."),paste("lambda",1:self$p,sep="."),paste("tau",1:self$p,sep="."))
      }

      # Boundary constraints
      lb = outer(rep(1,self$p),lower,"*")
      colnames(lb) = c("lambda","delta","tau")
      ub = outer(rep(1,self$p),upper,"*")
      colnames(ub) = c("lambda","delta","tau")

      # Fix parameters to their current value, if required
      if("lambda" %in% fixed) lb[,"lambda"] <- ub[,"lambda"] <- self$lambda
      if("delta" %in% fixed)  lb[,"delta"]  <- ub[,"delta"]  <- self$delta
      if("tau" %in% fixed)    lb[,"tau"]    <- ub[,"tau"]    <- self$tau

      # Main loop
      for(iter in 1:niter) {
        
        # Plot current solution
        if(display) {

          # Displays the relationship to true values if provided
          if(!is.null(private$true.par)) self$Plot(type="simulation")

          # Otherwise displays as a function of estimates
          else {
            self$Plot(type=plot.type,mfrow=mfrow,...)
            if(save.plot && (plot.type=="3D")) {
              rgl.snapshot(paste(save.path,"/bum",formatC(t,digits=4,flag="0"),".png",sep=""))
            }
          }
        }

        # Item update
        ll = 0
        for(item in 1:self$p) {

          if(verbose) cat(sprintf('\rJML cycle %d/%d — item %d/%d   ', iter, niter, item, self$p))
          params = paramsTable[item,]
          names(params) = c("lambda","delta","tau")

          newpar = optim(params,private$itemLogLikFunction,private$itemGradient,control=control,method="L-BFGS-B",lower=lb[item,],upper=ub[item,],item=item)
          
          # Store results
          paramsTable[item,] = newpar$par
          ll = ll + newpar$value
        }
        
        self$lambda = paramsTable[,'lambda']
        self$delta  = paramsTable[,'delta']
        self$tau    = paramsTable[,'tau']

        # Person update
        if(verbose) {
          cat("\n")
          cat(sprintf('\rJML cycle %d/%d — persons       ', iter, niter))
          cat("\n")
        }
        private$JML.theta()

        # Center final solution so that theta distrib has 0 mean
        if(!("delta" %in% fixed)) {
        
          m = mean(self$theta)
          self$theta = self$theta - m
          self$delta = self$delta - m
          
        }
        
        private$target = c(private$target,-ll)
        it = length(private$target)
        if( (iter>1) && (((private$target[it]-private$target[it-1])/private$target[it]) < private$reltol) ) break

        private$logLik = c(private$logLik,private$target[it])   
      }

      private$end.time = proc.time()
      delay = (private$end.time - private$start.time)[3]
      minutes = floor(delay/60)
      seconds = round(delay - (minutes*60))

      private$converged=(iter<niter)
            
      # Number of parameters
      private$npar = self$N + 3*self$p -1
          
      # Correction if lambdas fixed
      if("lambda" %in% fixed) private$npar = private$npar - self$p
      
      # Correction if deltas fixed
      if("delta" %in% fixed) private$npar = private$npar - self$p
      
      # Correction if taus fixed
      if("tau" %in% fixed) private$npar = private$npar - self$p

    },

    itemEMFunction = function(params, item) {

      DA = outer(rep(1,self$N),(self$theta.k - params['delta']),"*")
      x = self$x[,item]
      X = outer(x, rep(1,self$K),"*")

      m = exp(params['lambda'])
      n = exp(DA+params['tau']) + exp(-DA+params['tau'])

      if(self$tobit) {
        left  = (x <= self$eps)
        right = (x >= 1 - self$eps)
      } else {
        left  = rep(FALSE, self$N)
        right = rep(FALSE, self$N)
      }
      interior = !(left | right)

      # Log-likelihood contributions (N x K matrix)
      LL = matrix(0, self$N, self$K)

      if(any(interior)) LL[interior,] = dbeta(X[interior,,drop=FALSE], m, n[interior,,drop=FALSE], log=TRUE)
      if(any(left))     LL[left,]     = t(sapply(which(left), function(i) pbeta(self$eps, m, n[i,], log.p=TRUE)))
      if(any(right))    LL[right,]    = t(sapply(which(right), function(i) pbeta(1-self$eps, m, n[i,], lower.tail=FALSE, log.p=TRUE)))

      # Negative expected log-posterior
      nll = -sum(LL * private$posterior)
      nll = nll - dnorm(params['delta'],  self$prior.delta$mean,  self$prior.delta$sd,  log=TRUE)
      nll = nll - dnorm(params['lambda'], self$prior.lambda$mean, self$prior.lambda$sd, log=TRUE)
      nll = nll - dnorm(params['tau'],    self$prior.tau$mean,    self$prior.tau$sd,    log=TRUE)
      nll

    },

    itemEMGradient = function(params, item) {

      D = outer(rep(1,self$N),self$theta.k - params['delta'],"*")
      DA = D
      x = self$x[,item]
      X = outer(x, rep(1,self$K),"*")

      m = exp(params['lambda'])
      n = exp(DA+params['tau']) + exp(-DA+params['tau'])

      if(self$tobit) {
        left  = (x <= self$eps)
        right = (x >= 1 - self$eps)
      } else {
        left  = rep(FALSE, self$N)
        right = rep(FALSE, self$N)
      }
      interior = !(left | right)

      # --- Gradient contributions from interior observations (N x K) ---
      G.lambda = matrix(0, self$N, self$K)
      G.delta  = matrix(0, self$N, self$K)
      G.tau    = matrix(0, self$N, self$K)

      if(any(interior)) {
        ii = which(interior)
        G.lambda[ii,] = -m * (digamma(m+n[ii,,drop=FALSE])-digamma(m)+log(X[ii,,drop=FALSE]))
        G.delta[ii,]  =  2*exp(params['tau']) * sinh(DA[ii,,drop=FALSE]) * (digamma(m+n[ii,,drop=FALSE])-digamma(n[ii,,drop=FALSE])+log(1-X[ii,,drop=FALSE]))
        G.tau[ii,]    = -n[ii,,drop=FALSE] * (digamma(m+n[ii,,drop=FALSE])-digamma(n[ii,,drop=FALSE])+log(1-X[ii,,drop=FALSE]))
      }

      # --- Censored observations: finite-difference gradient of log pbeta ---
      if(self$tobit) {
        h = 1e-5
        m1 = exp(params['lambda'] + h)

        if(any(left)) {
          for(i in which(left)) {
            lp0 = pbeta(self$eps, m, n[i,], log.p=TRUE)

            lp1 = pbeta(self$eps, m1, n[i,], log.p=TRUE)
            G.lambda[i,] = -(lp1 - lp0) / h

            n1 = exp(DA[i,]+h+params['tau']) + exp(-(DA[i,]+h)+params['tau'])
            lp1 = pbeta(self$eps, m, n1, log.p=TRUE)
            G.delta[i,] = (lp1 - lp0) / h

            n1 = exp(DA[i,]+params['tau']+h) + exp(-DA[i,]+params['tau']+h)
            lp1 = pbeta(self$eps, m, n1, log.p=TRUE)
            G.tau[i,] = -(lp1 - lp0) / h
          }
        }

        if(any(right)) {
          for(i in which(right)) {
            lp0 = pbeta(1-self$eps, m, n[i,], lower.tail=FALSE, log.p=TRUE)

            lp1 = pbeta(1-self$eps, m1, n[i,], lower.tail=FALSE, log.p=TRUE)
            G.lambda[i,] = -(lp1 - lp0) / h

            n1 = exp(DA[i,]+h+params['tau']) + exp(-(DA[i,]+h)+params['tau'])
            lp1 = pbeta(1-self$eps, m, n1, lower.tail=FALSE, log.p=TRUE)
            G.delta[i,] = (lp1 - lp0) / h

            n1 = exp(DA[i,]+params['tau']+h) + exp(-DA[i,]+params['tau']+h)
            lp1 = pbeta(1-self$eps, m, n1, lower.tail=FALSE, log.p=TRUE)
            G.tau[i,] = -(lp1 - lp0) / h
          }
        }
      }

      # Weight by posterior and sum
      g.lambda = sum(G.lambda * private$posterior)
      g.delta  = sum(G.delta  * private$posterior)
      g.tau    = sum(G.tau    * private$posterior)

      # Add gradient of negative log-prior
      g.lambda = g.lambda + (params['lambda'] - self$prior.lambda$mean) / self$prior.lambda$sd^2
      g.delta  = g.delta  + (params['delta']  - self$prior.delta$mean)  / self$prior.delta$sd^2
      g.tau    = g.tau    + (params['tau']    - self$prior.tau$mean)    / self$prior.tau$sd^2

      c(lambda = g.lambda, delta = g.delta, tau = g.tau)

    },

    EM = function (lower = c(-5, -10, -5), upper = c(5, 10, 5), plot.type = "erf",
        fixed = NULL, unimodal.constraint = FALSE, save.plot = FALSE,
        save.history = TRUE, save.path = getwd(), display = TRUE, verbose = TRUE,
        mfrow = c(3, 3), start.method = "CA", prior = "gaussian", mu=0, sigma=1, theta.lim=c(-10,10), K = 81, ...) {

      # To constrain the response density, you have to fit the model uncontrained first, to detect BRD items
      if(unimodal.constraint && !private$start.config) {
        cat("Error: You have to fit the model uncontrained first, to detect BRD items.\n")
        return()
      }
      
      private$start.time = proc.time()
      private$method = "EM"
      private$logLik = NULL
      self$mfrow = mfrow
      private$fixed = fixed
      self$K = K
      private$fixed = fixed
      private$theta.lim = theta.lim

      # Clip boundary values unless tobit is enabled
      if(!self$tobit) {
        for(j in 1:self$p) {
          self$x[self$x[,j] <    private$small,j] = private$small
          self$x[self$x[,j] > (1-private$small),j] = 1 - private$small
        }
      }

      # Starting configuration
      if(!private$start.config)  {
      
        if(start.method == "CA") private$CA(verbose=verbose)
        else if(start.method == "random") private$randomParams()
      }
      
      init.delta  = self$delta
      init.lambda = self$lambda
      init.tau    = self$tau
      paramsTable = cbind(lambda=self$lambda,delta=self$delta,tau=self$tau)
      ones = rep(1,self$p)

      # Save history
      if(save.history) {
        self$history = rbind(self$history,matrix(c(1,0,self$delta,self$lambda,self$tau), 1, 3*self$p +2))
        colnames(self$history) = c("Algo","Iter",paste("delta",1:self$p,sep="."),paste("lambda",1:self$p,sep="."),paste("tau",1:self$p,sep="."))
      }

      # Set the fixed theta.k and their probability masses
      if(prior == "gaussian") private$setPrior("gaussian",mu,sigma,theta.lim,K)
      else                   private$setPrior("empirical")
      
      # Optimization options
      control = list(maxit=10)
      
      # Boundary constraints
      lb = outer(rep(1,self$p),lower,"*")
      colnames(lb) = c("lambda","delta","tau")
      ub = outer(rep(1,self$p),upper,"*")
      colnames(ub) = c("lambda","delta","tau")
      
      # Fix parameters to their current value
      if("lambda" %in% fixed) lb[,"lambda"] <- ub[,"lambda"] <- self$lambda
      if("delta" %in% fixed)  lb[,"delta"]  <- ub[,"delta"]  <- self$delta
      if("tau" %in% fixed)    lb[,"tau"]    <- ub[,"tau"] <- self$tau

      # Unimodality constraint if required
      if(unimodal.constraint) {
        is.bimodal = (self$lambda < 0) & (self$tau < -log(2))
        lambda.smallest = self$lambda >  log(2) + self$tau 
        
        select1 = is.bimodal & lambda.smallest
        select2 = is.bimodal & !lambda.smallest
        
        lb[select1,"lambda"] = 0
        lb[select2,"tau"] = -log(2)
        self$lambda[select1] = 0
        self$tau[select2] = -log(2)
      }

      # Main loop
      for(iter in 1:private$tmax) {
        
        # Plot current solution
        if(display) {

          # Displays the relationship to true values if provided
          if(!is.null(private$true.par)) self$Plot(type="simulation")

          # Otherwise displays as a function of estimates
          else {

          self$Plot(type=plot.type,mfrow=mfrow,...)

            if(save.plot && (plot.type=="3D")) {
              rgl.snapshot(paste(save.path,"/bum",formatC(t,digits=4,flag="0"),".png",sep=""))
            }
          }
        }

        # E-step
        P.ijk = private$Lijk(log=TRUE)
        private$posterior = exp(tensor(P.ijk,ones,2,1)) %*% diag(self$pi.k)
        private$logLik = c(private$logLik,sum(log(rowSums(private$posterior))))
        private$posterior = private$posterior / rowSums(private$posterior)

        # M-step
        if(prior == "empirical") self$pi.k = colMeans(private$posterior)

        ll = 0
        for(item in 1:self$p) {

          if(verbose) cat(sprintf('\rEM cycle %d/%d — item %d/%d   ', iter, private$tmax, item, self$p))
          params = paramsTable[item,]
          names(params) = c("lambda","delta","tau")

          newpar = optim(params,private$itemEMFunction,private$itemEMGradient,control=control,method="L-BFGS-B",lower=lb[item,],upper=ub[item,],item=item)
          paramsTable[item,] = newpar$par
          ll = ll + newpar$value
        }
        
        # Fix the metric (scale the theta.k) unless deltas or prior probs are fixed (not estimated)
        if(!( ("delta" %in% fixed) || (prior == "gaussian") )) {
          m = sum(self$theta.k * self$pi.k)
          self$theta.k = self$theta.k - m     
          s = sqrt(sum(self$theta.k**2 * self$pi.k))
          self$theta.k = self$theta.k / s
        }
        
        self$lambda = paramsTable[,'lambda']
        self$delta  = paramsTable[,'delta']
        self$tau    = paramsTable[,'tau']
      
        if(save.history) self$history = rbind(self$history,c(1,iter,self$delta,self$lambda,self$tau))

        # EAP estimates of theta
        self$theta = as.vector(private$posterior %*% self$theta.k)
        
        # Stopping criterion
        private$target = c(private$target,-ll)
        it = length(private$target)
        if( max(c(abs(init.delta  - self$delta),
                  abs(init.lambda - self$lambda),
                  abs(init.tau    - self$tau))) < private$tol ) break
        # if((iter>1) && (((private$target[it]-private$target[it-1])/private$target[it]) < private$reltol)) break

        init.delta  = self$delta
        init.lambda = self$lambda
        init.tau    = self$tau

      }
        
      # Store marginal loglikelihood
      P.ijk = private$Lijk(log=TRUE)
      private$logLik = c(private$logLik,sum(log(rowSums(exp(tensor(P.ijk,ones,2,1)) %*% diag(self$pi.k)))))
      
      # Computation time
      private$end.time = proc.time()
      delay = (private$end.time - private$start.time)[3]
      minutes = floor(delay/60)
      seconds = round(delay - (minutes*60))

    private$converged = (iter<private$tmax)
    
      if(verbose) {
        cat("\n")
        if(private$converged) cat("Solution converged in",iter,"iterations\n")
        else cat("Warning: Maximum number of iterations reached.\n")
        cat("Computation time:",minutes,"min.",seconds,"sec.\n\n")
      }
      
      # Theoretical number of parameters (-1 because of the location constraint on the thetas)
      private$npar = 3*self$p + (self$K-1) - 1
      
      # Correction if the pi.k are fixed and not estimated
      if(prior == "gaussian") private$npar = private$npar - (self$K-1)   # K-1 class probs fixed
      
      # Correction if lambdas fixed
      if("lambda" %in% fixed) private$npar = private$npar - self$p
      
      # Correction if deltas fixed
      if("delta" %in% fixed) private$npar = private$npar - self$p
      
      # Correction if taus fixed
      if("tau" %in% fixed) private$npar = private$npar - self$p
      
      # Correction if unimodality constraints
      if(unimodal.constraint) private$npar = private$npar - sum(is.bimodal)

    },
    Simulate = function (N = 500, p = 25, mu.theta = 0, sigma.theta = 1, delta = NULL, lambda = NULL, tau = NULL, delta.lim = c(-2,  
        2), lambda.lim = c(-1.5, 1.5), tau.lim = c(-3, 1), theta.scale = TRUE) {

      theta = sort(rnorm(N,mu.theta,sigma.theta))
      if(theta.scale) theta  = as.vector(scale(theta))
      
      if(is.null(delta))  delta  = seq(delta.lim[1],delta.lim[2],length=p)
      if(is.null(lambda)) lambda = runif(p,lambda.lim[1],lambda.lim[2])
      if(is.null(tau))    tau    = runif(p,tau.lim[1],tau.lim[2])
      
      stopifnot( (length(delta) == length(lambda)) && (length(delta) == length(lambda)) )

      p = length(delta) # Just in case specific parameter values were provided

      T = outer(rep(1,N),tau,"*")
      D = outer(theta,delta,"-")
      L = outer(rep(1,N),lambda,"*")
      
      self$x = matrix(0,N,p)
      m = exp(L)
      n = exp(D+T)+exp(-D+T)
      
      for(i in 1:N)
        for(j in 1:p)
          self$x[i,j] = rbeta(1,m[i,j],n[i,j])

      colnames(self$x) = paste("Item",1:p,sep="")
      
      private$true.par = list(theta=theta,delta=delta,lambda=lambda,tau=tau)
      self$N = N
      self$p = p
      private$logLik = NULL

    },
    jointLogLikelihood = function() {

      T = outer(rep(1,self$N),self$tau,"*")
      D = outer(self$theta,self$delta,"-")
      L = outer(rep(1,self$N),self$lambda,"*")

      m = exp(L)
      n = exp(D+T)+exp(-D+T)

      # Censored joint log-likelihood
      ll = 0
      for(j in 1:self$p) {
        cl = private$censoredLogLik(self$x[,j], m[,j], n[,j])
        ll = ll + sum(cl$ll)
      }
      ll

    },  
    GS = function (verbose = TRUE, display = FALSE) {

      # Estimate theta by correspondence analysis
      if(is.null(self$theta)) private$CA(verbose=verbose)

      delta  = seq(-4,4,len=8)
      lambda = seq(-3,3,len=8)
      tau    = seq(-3,1,len=8)
      
      gr = as.matrix(expand.grid(delta=delta,lambda=lambda,tau=tau))
      nvec = nrow(gr)
      M = matrix(1,self$p,4)
      colnames(M) = c("delta","lambda","tau")
      
      # Do a grid search over lambda and gamma2 for all items
      if(verbose) {
        cat("\n -----------\n")
        cat(" Grid search\n")
        cat(" -----------\n")
      }
      
      for(i in 1:self$p) {
        L = private$itemJointLogLikelihoodFunction(list(delta=gr[,"delta"],lambda=gr[,"lambda"],tau=gr[,"tau"],x=self$x[,i]))
        Lmax = max(L,na.rm=TRUE)
        if(is.finite(Lmax)) M[i,] = gr[which(L==Lmax,arr.ind=TRUE),]
        if(verbose) cat("   - Item",i,"done\n")
      }

      self$lambda = M[,"lambda"]
      self$delta  = M[,"delta"]
      self$tau    = M[,"tau"]

      private$JML.theta(verbose=verbose)
    
      # FIx the origin
      m = mean(self$theta)
      self$theta = (self$theta-m)
      self$delta = (self$delta-m)
    
      private$start.config = TRUE
    
      if(verbose) cat("  OK.\n")
      if(display) self$Plot()

    },

    CA = function(verbose = TRUE) {

      if(verbose) {
        cat("Initial estimates (correspondence analysis)\n")
      }
      
      library(MASS)
      options(warn=-1)
      
      # Estimate theta by correspondence analysis
      cc = corresp(self$x)

      if(is.null(self$theta)) {
        # Item locations are fixed: Set initial thetas as an weighted average
        if("delta" %in% private$fixed) self$theta = as.vector((self$x %*% self$delta) / rowSums(self$x))
        else                           self$theta = cc$rscore
      }

      # Standardize theta and rescale deltas accordingly
      if((!("delta" %in% private$fixed)) && is.null(self$delta)) {
        m = mean(self$theta)
        s = sd(self$theta)

      self$theta = (self$theta-m)/s
      self$theta[self$theta < -4.5] = -4.5
      self$theta[self$theta >  4.5] =  4.5

      self$delta = cc$cscore
      self$delta = (self$delta-m)/s
      self$delta[self$delta < -4.5] = -4.5
      self$delta[self$delta >  4.5] =  4.5
      }
      
      if(!is.null(private$true.par)) {
        # Just in case initial solution is just of the wrong sign in a simulation study
        # (loadings signs are arbitrary in CA)
        if(cor(self$theta,private$true.par$theta)<0) { 
        self$theta = -self$theta
        self$delta = -self$delta
        }
      }

      if(is.null(self$lambda)) self$lambda = rep(1,self$p)
      if(is.null(self$tau))    self$tau    = rep(0,self$p)
    
      # private$logLik = private$jointLogLikelihood()
      # private$target = private$logLik
      private$start.config = TRUE

    },

    Vcov = function () {

      D = outer(self$theta,self$delta,"-")
      DA = D
      L = outer(rep(1,self$N),self$lambda,"*")
      T = outer(rep(1,self$N),self$tau,"*")

      m = exp(L)
      n = exp(DA+T)+exp(-DA+T)

      private$dsec = array(0,dim=c(self$p,3,3),list(colnames(self$x),c("lambda","delta","tau"),c("lambda","delta","tau")))

      private$dsec[,"lambda","lambda"] =  exp(2*self$lambda)*colSums( trigamma(m+n)-trigamma(m) )
      private$dsec[,"delta","delta"]   =  4*exp(2*self$tau)*colSums( sinh(DA)**2 * (trigamma(m+n)-trigamma(n)))
      private$dsec[,"tau","tau"]       =  colSums( n**2 * (trigamma(m+n)-trigamma(n)) )

      private$dsec[,"lambda","delta"]  = -2*exp(self$tau+self$lambda)*colSums(sinh(DA)*trigamma(m+n))
      private$dsec[,"lambda","tau"]    =  exp(self$lambda)*colSums(n*trigamma(m+n))

      private$dsec[,"delta","tau"]     = -2*exp(self$tau)*colSums(n*sinh(DA)*(trigamma(m+n)-trigamma(n)))

      
      # Symmetric elements
      private$dsec[,"delta","lambda"] = private$dsec[,"lambda","delta"]
      private$dsec[,"tau","lambda"]   = private$dsec[,"lambda","tau"]
      private$dsec[,"tau","delta"]    = private$dsec[,"delta","tau"]

      # Inverting the matrices
      vcov = private$dsec
      for(i in 1:self$p) vcov[i,,] = solve(-private$dsec[i,,])
      vcov

    },

    Corr = function () {

      # Inverting the matrices
      v.cov = Vcov()
      R = v.cov
      for(i in 1:self$p) {
        s = sqrt(diag(v.cov[i,,]))
        R[i,,] = v.cov[i,,] / outer(s,s,'*')
      }
      
      R

    },

    ASE = function () {

      v.cov = private$Vcov()
      ase = matrix(0,self$p,nrow(v.cov[1,,]),dimnames=list(colnames(self$x),colnames(v.cov[1,,])))
      for(i in 1:self$p) ase[i,] = sqrt(diag(v.cov[i,,]))
      
      ase

    },  
    JMLupdateThetaDeriv = function() {

      T = outer(rep(1,self$N),self$tau,"*")
      D = outer(self$theta,self$delta,"-")
      L = outer(rep(1,self$N),self$lambda,"*")

      m = exp(L)
      n = exp(D+T)+exp(-D+T)

      if(self$tobit) {
        # Per-item loop with censoring
        h = 1e-5
        deriv1 = matrix(0, self$N, self$p)
        deriv2 = matrix(0, self$N, self$p)

        for(j in 1:self$p) {
          x = self$x[,j]
          left  = (x <= self$eps)
          right = (x >= 1 - self$eps)
          interior = !(left | right)

          if(any(interior)) {
            deriv1[interior,j] = 2*exp(self$tau[j])*sinh(D[interior,j])*(digamma(m[interior,j]+n[interior,j])-digamma(n[interior,j])+log(1-x[interior]))
            deriv2[interior,j] = 4*(exp(self$tau[j])*sinh(D[interior,j]))^2 * (trigamma(m[interior,j]+n[interior,j])-trigamma(n[interior,j]))
          }

          if(any(left)) {
            ii = which(left)
            lp0 = pbeta(self$eps, m[ii,j], n[ii,j], log.p=TRUE)
            D1 = D[ii,j] + h
            n1 = exp(D1+self$tau[j]) + exp(-D1+self$tau[j])
            lp1 = pbeta(self$eps, m[ii,j], n1, log.p=TRUE)
            deriv1[ii,j] = (lp1 - lp0) / h
            Dm = D[ii,j] - h
            nm = exp(Dm+self$tau[j]) + exp(-Dm+self$tau[j])
            lpm = pbeta(self$eps, m[ii,j], nm, log.p=TRUE)
            deriv2[ii,j] = (lp1 - 2*lp0 + lpm) / h^2
          }

          if(any(right)) {
            ii = which(right)
            lp0 = pbeta(1-self$eps, m[ii,j], n[ii,j], lower.tail=FALSE, log.p=TRUE)
            D1 = D[ii,j] + h
            n1 = exp(D1+self$tau[j]) + exp(-D1+self$tau[j])
            lp1 = pbeta(1-self$eps, m[ii,j], n1, lower.tail=FALSE, log.p=TRUE)
            deriv1[ii,j] = (lp1 - lp0) / h
            Dm = D[ii,j] - h
            nm = exp(Dm+self$tau[j]) + exp(-Dm+self$tau[j])
            lpm = pbeta(1-self$eps, m[ii,j], nm, lower.tail=FALSE, log.p=TRUE)
            deriv2[ii,j] = (lp1 - 2*lp0 + lpm) / h^2
          }
        }

        private$thetaDeriv1 = rowSums(deriv1)
        private$thetaDeriv2 = rowSums(deriv2)

      } else {
        # Standard analytical derivatives (no censoring, data already clipped)
        private$thetaDeriv1 = 2*rowSums( (exp(T)*sinh(D))*(digamma(m+n)-digamma(n)+log(1-self$x)) )
        private$thetaDeriv2 = 4*rowSums( ((exp(T)*sinh(D))**2) * (trigamma(m+n)-trigamma(n)) )
      }

    },
    JML.theta = function (c1 = 1e-04, rho = 0.5, tol = 1e-06, itmax = 30, verbose = FALSE) {

      # Initialize
      last.theta = self$theta

      for(k in 1:itmax)
      {
        # Compute first and second derivatives
      private$JMLupdateThetaDeriv()

        # Derivatives are close to zero or infinite
        ng = private$norm(private$thetaDeriv1)
        nh = private$norm(private$thetaDeriv2)
        if( !is.finite(ng) || (ng<tol) || !is.finite(nh)  || (nh<tol) ) break

        p = - private$thetaDeriv1/private$thetaDeriv2
        slope = t(-private$thetaDeriv1) %*% p

        if( (slope==0) || (!is.finite(slope)) ) break

        # Going uphill ?
        if(slope > 0)
        {
          p = private$thetaDeriv1	# Use steepest descent instead
          if(verbose) cat("Newton: Using steepest descent direction\n")
        }

        # Check decrease
        a = 1
        fval0 = - private$jointLogLikelihood()
        self$theta =  self$theta + a*p
        fval1 = - private$jointLogLikelihood()

        # Backtrack while criterion not satisfied
        while(is.na(fval1) || (fval1 > fval0 + c1*a*slope))
        {
        self$theta = self$theta - a*p           # Cancel previous move
          a = rho*a                        # Step halving
          if(a < .0625) break              # Too small a move
          self$theta =  self$theta + a*p          # Retry a smallest move
          fval1 = - private$jointLogLikelihood() # Check EM target function
        }

        if(a < .0625) break

        if(verbose) cat(k,"    Newton step: a =",a,"\n")

        if(max(abs(self$theta - last.theta))<tol) break
        last.theta = self$theta
      }

    },
    Lijk = function(log = FALSE) {

      # Array of densities
      P = array(0,c(self$N,self$p,self$K))

      for(k in 1:self$K) {

        T = outer(rep(1,self$N),self$tau,"*")
        D = outer(rep(self$theta.k[k],self$N),self$delta,"-")
        L = outer(rep(1,self$N),self$lambda,"*")

        m = exp(L)
        n = exp(D+T)+exp(-D+T)

        if(self$tobit) {
          for(j in 1:self$p) {
            x = self$x[,j]
            left  = (x <= self$eps)
            right = (x >= 1 - self$eps)
            interior = !(left | right)

            if(log) {
              if(any(interior)) P[interior,j,k] = dbeta(x[interior], m[interior,j], n[interior,j], log=TRUE)
              if(any(left))     P[left,j,k]     = pbeta(self$eps,     m[left,j],     n[left,j],     log.p=TRUE)
              if(any(right))    P[right,j,k]    = pbeta(1-self$eps,   m[right,j],    n[right,j],    lower.tail=FALSE, log.p=TRUE)
            } else {
              if(any(interior)) P[interior,j,k] = dbeta(x[interior], m[interior,j], n[interior,j])
              if(any(left))     P[left,j,k]     = pbeta(self$eps,     m[left,j],     n[left,j])
              if(any(right))    P[right,j,k]    = pbeta(1-self$eps,   m[right,j],    n[right,j],    lower.tail=FALSE)
            }
          }
        } else {
          P[,,k] = dbeta(self$x,m,n,log=log)
        }
      }

      P

    },
    Reset = function () {

      private$start.config=FALSE
      private$logLik = NULL
      self$theta = NULL
      self$delta = NULL
      self$lambda = NULL
      private$gamma = NULL
      self$tau = NULL

    },
    # Predict from a true, simulated model
    Predict2 = function (tt = NULL, conf.level = 0.9) {

      if(is.null(private$true.par)) return()
      if(is.null(tt)) tt = seq(min(private$true.par$theta),max(private$true.par$theta),length=100)
      N = length(tt)

      T = outer(rep(1,N),private$true.par$tau,"*")
      D = outer(tt,private$true.par$delta,"-")
      L = outer(rep(1,N),private$true.par$lambda,"*")

      m = exp(L)                          # Acceptation
      n = exp(D+T)+exp(-D+T)              # Refusal

      mu = m / (m+n)                      # Expectation
      colnames(mu) = colnames(self$x)

      md = (m-1) / (m+n-2)                # Mode
      colnames(md) = colnames(self$x)
      md[ (m<1) & (n>1) ] = 0
      md[ (m>1) & (n<1) ] = 1
      md[ (m<1) & (n<1) & (m<n) ] = 0
      md[ (m<1) & (n<1) & (m>n) ] = 1
      
      med = qbeta(.50,m,n)                # Median
      colnames(med) = colnames(self$x)
      q5 = qbeta((1-conf.level)/2,m,n)    # 5th percentile or lower bound of confidence interval
      colnames(q5) = colnames(self$x)
      q95 = qbeta((1+conf.level)/2,m,n)   # 95th percentile or upper bound of confidence interval
      colnames(q95) = colnames(self$x)
      sigma2 = mu * (1-mu) / (m + n + 1)  # Variance
      colnames(mu) = colnames(self$x)

      list(theta=tt,m=m,n=n,mu=mu,mode=md,median=med,q5=q5,q95=q95,sigma2=sigma2)

    },

    itemJointLogLikelihoodFunction = function(sol) {

      p = length(sol$delta)

      T = outer(rep(1,self$N),sol$tau,"*")
      D = outer(self$theta,sol$delta,"-")
      L = outer(rep(1,self$N),sol$lambda,"*")
      X = outer(sol$x,rep(1,p),"*")

      m = exp(L)
      n = exp(D+T)+exp(-D+T)

      # Censored log-likelihood per item
      result = numeric(p)
      for(j in 1:p) {
        cl = private$censoredLogLik(X[,j], m[,j], n[,j])
        result[j] = sum(cl$ll)
      }
      result

    },

    thetaPosteriorSD = function() {

      sqrt(private$posterior %*% (self$theta.k**2) - (private$posterior %*% self$theta.k)**2)

    },

    Sample = function(N = NULL, mu.theta = 0, sigma.theta = 1, theta.scale = TRUE) {

      if(is.null(true.par)) {
        cat("No parameter values provided.\n")
        return()
      }
      
      # Resample with a new sample size
      if(!is.null(N)) {
      self$N = N
      private$true.par$theta = sort(rnorm(N,mu.theta,sigma.theta))
        if(theta.scale) private$true.par$theta  = as.vector(scale(private$true.par$theta))
      }
      
      T = outer(rep(1,self$N),private$true.par$tau,"*")
      D = outer(private$true.par$theta,private$true.par$delta,"-")
      L = outer(rep(1,self$N),private$true.par$lambda,"*")
      
      m = exp(L)
      n = exp(D+T)+exp(-D+T)
      
      self$x = matrix(0,self$N,self$p) 
      for(i in 1:self$N)
        for(j in 1:self$p)
          self$x[i,j] = rbeta(1,m[i,j],n[i,j])

      names(self$x) = paste("Item",1:self$p,sep="")
      private$start.config = FALSE
      self$theta = NULL

    },

    randomParams = function(mu.theta = 0, sigma.theta = 0.5, delta.lim = c(-2, 2), lambda.lim = c(-1, 1), tau.lim = c(-3, 1), theta.scale = TRUE) {

      self$theta = rnorm(self$N,mu.theta,sigma.theta)
      if(theta.scale) self$theta  = as.vector(scale(self$theta))*sigma.theta + mu.theta
      
      self$delta  = runif(self$p,delta.lim[1],delta.lim[2])
      self$lambda = runif(self$p,lambda.lim[1],lambda.lim[2])
      self$tau    = runif(self$p,tau.lim[1],tau.lim[2])
      
      private$start.config = TRUE  

    },

    setPrior = function(type = "empirical", mu=0, sigma=4, theta.lim=c(-10,10), K=81) {
    
      if(type == "uniform") {
        self$theta.k = seq(private$theta.lim[1],private$theta.lim[2],len=self$K)
        self$pi.k = rep(1/self$K, self$K)    
      }
      else if(type == "gaussian") {

        brks = seq(theta.lim[1],theta.lim[2],len=K + 1)
        self$pi.k = diff(pnorm(brks,mu,sigma))
        self$pi.k = self$pi.k / sum(self$pi.k)
        mids = (brks[1:K]+brks[2:(K+1)])/2

        m = sum(mids * self$pi.k)
        mids = mids - m
        s = sqrt(sum(mids**2 * self$pi.k))
        mids = mids/s
        self$theta.k = mids * sigma + mu
      }
      
      else if(type == "empirical") {
        self$theta.k = seq(private$theta.lim[1],private$theta.lim[2],len=self$K)
        stopifnot(private$start.config)
        self$theta[self$theta < -49] = -49
        self$theta[self$theta >  49] =  49
        breaks = (self$theta.k[1:(self$K-1)]+self$theta.k[2:self$K])/2
        breaks = c(-50,breaks,50)
        self$pi.k = hist(self$theta,breaks=breaks,plot=F)$counts
        self$pi.k = self$pi.k / sum(self$pi.k)
        m = sum(self$theta.k * self$pi.k)
        s = sqrt(sum(self$theta.k**2 * self$pi.k) - m**2)
        self$theta.k = (self$theta.k - m)/s
      }

    },
    norm = function (v) {

      sqrt(sum(v*v))

    }, 
    # Estimation parameters
    tmax = 300, 
    npar = 0,
    method = NULL,
    tol = 0.01, 
    reltol = sqrt(.Machine$double.eps), 
    theta.lim = c(-4.5,4.5), 
    fixed = NULL,
    prior = "gaussian",
    # Model parameters
    true.par = NULL, 
    dsec = array(),
    # Estimation results
    converged = FALSE,
    start.time = 0, 
    end.time = 0, 
    target = NULL, 
    small = 1e-08, 
    posterior = NULL, 
    thetaDeriv1 = NULL, 
    thetaDeriv2 = NULL,
    start.config = FALSE, 
    unimodal.constraint = FALSE,
    fixed.slopes = TRUE,
    logLik = NULL,
    inverted = NULL
  )
)
