##################################################################
# GS_CopulaCoursePhD_R
##################################################################
# Original code by Gianfausto Salvadori (Universita del Salento)
# Edited further by Damian Smug
##################################################################

GS_CopulaCoursePhD_R <- function(choose_copula,OptMeth,TxtLbl,data_file,survival,standardise){
	# TxtLbl: text label for naming results

	# load library
	library(copula)
	library(graphics) 
	library(grDevices)
	
	# plot parameter(s) to be fixed manually
	gr.npt <- 65 # n. of contour grid points (i.e., gr.npt-1 intervals)
	gr.lin <- seq(from=0,to=1,length.out=gr.npt) # linear grid
	gr.mat <- expand.grid(gr.lin,gr.lin) # grid matrix
	Cont.lev <- seq(from=0.1,to=0.9,by=0.2) # contour-plot levels
	
	# load data
	dat <- read.csv(data_file,header=FALSE)
	dat = dat[,1:2] #take into account only two first entries of data file
	# print(dat)
	
	N <- nrow(dat)
	print(paste("Sample size=",N,sep=""))
	data.dim <- ncol(dat)
	print(paste("Number of variables=",data.dim,sep=""))

	X <- dat[,1] # tau1
	Y <- dat[,2] # tau2
	
	if (standarise){
	  X = X/mean(X)
	  Y = Y/mean(Y)
	}
	
	# data summary
	print(paste("Summary X",sep=""))
	summary(X)
	print(paste("Summary Y",sep=""))
	summary(Y)
	
	# Ties (repeated values)
	X.tie.n <- N-length(unique(X))
	print(paste("X.tie.n=",X.tie.n,sep=""))
	Y.tie.n <- N-length(unique(Y))
	print(paste("Y.tie.n=",Y.tie.n,sep=""))
	
	### Kendall Tau
	KT.test <- cor.test(x=X,y=Y,method="kendall")
	KT.est <- KT.test$estimate
	KT.pv <- KT.test$p.value
	print(paste("Kendall Tau estimate=",KT.est,sep=""))
	print(paste("Kendall Tau p-value=",KT.pv,sep=""))
	
	### Spearman Rho
	SR.test <- cor.test(x=X,y=Y,method="spearman")
	SR.est <- SR.test$estimate
	SR.pv <- SR.test$p.value
	print(paste("Spearman Rho estimate=",SR.est,sep=""))
	print(paste("Spearman Rho p-value=",SR.pv,sep=""))
	
	### pseudo-observations [result may depend on the randomization of the pobs]
	u <- pobs(cbind(X,Y),ties.method="random")
	# in the case of analysing survival copulas one need to revert the 
	if (survival){
	  u = 1-u
	}

	# plots (EPS file)
	Cop.emp <- C.n(u=as.matrix(gr.mat),X=u) # empirical copula over the grid

	setEPS() # Data
	postscript(paste("figures/",substr(TxtLbl,1,4),"_Data.eps",sep=""))
	plot(x=X,y=Y,type="p",pch=20,col="black",
		main=paste("Data",sep=""),
		xlab=expression(tau^(i)),ylab=expression(tau^(j)))
	grid(nx=NULL,ny=NULL,lwd=1)
	dev.off()

	setEPS() # Pseudo-observations
	postscript(paste("figures/",substr(TxtLbl,1,4),"_Data_PObs.eps",sep=""))
	par(pty="s")
	plot(x=u[,1],y=u[,2],type="p",pch=20,col="black",
		xlim=c(0,1),ylim=c(0,1),asp=1,
		main=paste("Pseudo-observations",sep=""),
		xlab=expression(F[X]),ylab=expression(F[Y]))
	contour(gr.lin,gr.lin,matrix(Cop.emp,nrow=gr.npt,ncol=gr.npt),levels=Cont.lev,
		col="black",lty="solid",lwd=3,drawlabels=FALSE,add=TRUE)
	grid(nx=NULL,ny=NULL,lwd=1)
	dev.off()
	
	#######################################################
	### select a copula family [uncomment the desired line]

	if (choose_copula=='amh'){
	### family: amh (Ali-Mikhail-Haq)
	Cop <- archmCopula(family="amh",dim=data.dim)
	} else if (choose_copula=='clayton'){
	# ### family: clayton
	Cop <- archmCopula(family="clayton",dim=data.dim)
	} else if (choose_copula=='frank'){
	# ### family: frank
	Cop <- archmCopula(family="frank",dim=data.dim)
	} else if (choose_copula=='gumbel'){
	# ### family: gumbel
	Cop <- archmCopula(family="gumbel",dim=data.dim)
	} else if (choose_copula=='joe'){
	# ### family: joe
	Cop <- archmCopula(family="joe",dim=data.dim)
	} else if (choose_copula=='normal'){
	# ### family: normal
	Cop <- ellipCopula(family="normal",dim=data.dim)
	} else if (choose_copula=='t'){
	# ### family: t
	Cop <- ellipCopula(family="t",dim=data.dim,df.fixed=TRUE) # df must be fixed for t copula
	} else if (choose_copula=='galambos'){
	# ### family: galambos
	Cop <- evCopula(family="galambos",dim=data.dim)
	} else if (choose_copula=='huslerReiss'){
	# ### family: huslerReiss
	Cop <- evCopula(family="huslerReiss",dim=data.dim)
	} else if (choose_copula=='tawn'){
	# ### family: tawn
	Cop <- evCopula(family="tawn",dim=data.dim)
	} else if (choose_copula=='fgm'){
	# ### family: fgm (Farlie-Gumbel-Morgenstern)
	Cop <- fgmCopula(dim=data.dim)
	} else if (choose_copula=='plackett'){
	# ### family: plackett
	Cop <- plackettCopula()
	}
	################################################
	# The Khoudraji device could also be tried: e.g.
	# Cop.k <- khoudrajiCopula(copula1=Cop,copula2=indepCopula(),shapes=c(0.5,0.5))

	#############################
	# starting par(s) estimate(s)
	s <- iTau(copula=Cop,tau=KT.est)
	
	############################
	# Fit, GoF, cross-Validation
	
	# OptMeth <- "L-BFGS-B" # choose among "Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"
	Cop.fit <- fitCopula(copula=Cop,data=u,method="mpl",start=s,optim.method=OptMeth)
	print(paste("*** Fitted copula=",sep=""))
	print(paste(attributes(Cop.fit),sep=""))
	print(paste("Par(s) Est(s)=",attributes(Cop.fit)$estimate,sep=""))
	print(paste("Par(s) Var(s)=",attributes(Cop.fit)$var.est,sep=""))

	# OptMeth <- "Nelder-Mead" # choose among "Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"
	Cop.gof <- gofCopula(copula=Cop,x=u,method="Sn",estim.method="mpl",
		start=attributes(Cop.fit)$estimate,optim.method=OptMeth)$p.value
	print(paste("GoF p-Value=",Cop.gof,sep=""))

	# # OptMeth <- "Nelder-Mead" # choose among "Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"
	# Cop.xv <- xvCopula(copula=Cop,x=u,method="mpl",
	# 	start=attributes(Cop.fit)$estimate,optim.method=OptMeth)
	# print(paste("xV=",Cop.xv,sep=""))

	### plot fitted copula (EPS file)
	Cop.fig <- pCopula(u=as.matrix(gr.mat),copula=attributes(Cop.fit)$copula) # fitted copula
	
	setEPS()
	postscript(paste("figures/",substr(TxtLbl,1,4),"_Copula_",ifelse(survival,"surv_",""),ifelse(standardise,"stndrd_",""),substr(TxtLbl,6,nchar(TxtLbl)),".eps",sep=""))
	# no "best_"
	par(pty="s")
	
	
	# plot(x=u[,1],y=u[,2],type="p",pch=20,col="black",xlim=c(0,1),ylim=c(0,1),asp=1,
	# 	main=paste("Fit: ",ifelse(survival,"survival ",""),attributes(attributes(Cop.fit)$copula)$class[1],"\np-value=",signif(Cop.gof,digits=6),"\nxV=",signif(Cop.xv,digits=6),sep=""),
	# 	xlab=expression(F[X]),ylab=expression(F[Y])) #with cross-validation
	
	plot(x=u[,1],y=u[,2],type="p",pch=20,col="black",xlim=c(0,1),ylim=c(0,1),asp=1,
	     main=paste("Fit: ",ifelse(survival,"survival ",""),attributes(attributes(Cop.fit)$copula)$class[1],"\np-value=",signif(Cop.gof,digits=6),sep=""),
	     xlab=expression(F[X]),ylab=expression(F[Y])) #without cross-validation
	
	# plot(x=u[,1],y=u[,2],type="p",pch=20,col="black",xlim=c(0,1),ylim=c(0,1),asp=1,
	#      main="",#no title at all - for publishing
	#      xlab=expression(F[X]),ylab=expression(F[Y]),
	#      cex=1,cex.lab=2,cex.axis=2) 
	
	
	contour(gr.lin,gr.lin,matrix(Cop.emp,nrow=gr.npt,ncol=gr.npt),levels=Cont.lev,
		col="blue",lty="solid",lwd=4,drawlabels=FALSE,add=TRUE) #lwd=3 or 6
	contour(gr.lin,gr.lin,matrix(Cop.fig,nrow=gr.npt,ncol=gr.npt),levels=Cont.lev,
		col="red",lty="solid",lwd=2,drawlabels=FALSE,add=TRUE) #lwd=1 or 5
	grid(nx=NULL,ny=NULL,lwd=1)
	dev.off()
		
	############
	return( list( theta=attributes(Cop.fit)$estimate, pv=Cop.gof, tau=KT.est, rho=SR.est ) )
}

