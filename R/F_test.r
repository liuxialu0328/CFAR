#' @export
F_test <- function(f,p.max=6,df_b=10,grid=1000){
	# F test: CFAR(0)vs CFAR(1), CFAR(1)vs CFAR(2), CFAR(2)vs CFAR(3) when rho is known.
	if(!is.matrix(f))f <- as.matrix(f)
	t <- nrow(f)
	n=dim(f)[2]-1
	if(is.null(p.max)){
		p.max=6
	}
	if(is.null(df_b)){
		df_b=10
	}
	if(is.null(grid)){
		grid=1000
	}
	##################################################################################
	###parameter setting
	# rho=5	### parameter for error process epsilon_t, O-U process
	# t= 1000;	### length of time
	# iter= 100;	### number of replications in the simulation
	# grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t
	# df_b=10		### number of degrees for splines
	# n=100		### number of observations for each X_t, in the paper it is 'N'
	x_grid= seq(0, 1, by= 1/grid);	### the grid points for input variables of functional time series in [0,1]
	b_grid= seq(-1, 1, by=1/grid);	### the grid points for input variables of b-spline functions in [-1,1]
	eps_grid= matrix(0, 1, grid+1);	### the grid points for input variables of error process in [0,1]
	coef_b=df_b+1;				### number of df for b-spline
	b_grid_sp= ns(b_grid, df=df_b);	###b-spline functions
	x= seq(0, 1, by=1/n);			###observations points for X_t process
	index= 1:(n+1)*(grid/n)-(grid/n)+1
	x_grid_sp= bs(x_grid, df=n, degree=1);	###basis function values if we interpolate discrete observations of x_t
	x_sp= x_grid_sp[index,];			###basis function values at each observation point
	df=n;			### number of interplolation of X_t
	coef= df+1;
	###convolutions of b-pline functions (df=df_b) for phi and basis functions (df=n) for x
	x_grid_full= cbind(rep(1,grid+1), x_grid_sp);
	b_grid_full= cbind(rep(1,2*grid+1), b_grid_sp);
	b_matrix=matrix(0,grid+1,grid+1)
	bstar_grid2= matrix(0, grid+1, coef*coef_b)
	for (j in 1:(coef_b)){
		for (i in 1:(grid+1)){
			b_matrix[i,]= b_grid_full[seq((i+grid),i, by=-1),j]/(grid+1);
		}
		bstar_grid2[, (j-1)*coef+(1:coef)]= b_matrix %*% x_grid_full
	}
	bstar2= bstar_grid2[index,]
	bstar3= matrix(0, (n+1)*coef_b, coef)
	for (i in 1:coef_b){
		bstar3[(i-1)*(n+1)+(1:(n+1)),]=bstar2[, (i-1)*coef+(1:coef)]
	}
	dsg=cbind(rep(1,n+1),x_sp);
	dsg_mat= solve(t(dsg)%*%dsg)%*%t(dsg)
	ahat= f%*%t(dsg_mat)
	indep=matrix(0,(n+1)*(t-1),coef_b)
	for (i in 2:t){	
		M_t= matrix(ahat[i-1,]%*% t(bstar3), nrow=n+1, ncol=df_b+1)
		indep[(i-2)*(n+1)+(1:(n+1)),]=M_t
	}
	######### AR(0) and AR(1)
	#########################################
	#######################
	####cat("Data set: ",q,"\n")
	### Fit cfar(1) model
	###	f= f_grid[(q-1)*tmax+(1:t)+(tmax-t), index];
	p=1
	dsg=cbind(rep(1,n+1),x_sp);
	dsg_mat= solve(t(dsg)%*%dsg)%*%t(dsg)
	ahat= f%*%t(dsg_mat)
	indep=matrix(0,(n+1)*(t-1),coef_b)
	for (i in 2:t){	
		M_t= matrix(ahat[i-1,]%*% t(bstar3), nrow=n+1, ncol=df_b+1)
		indep[(i-2)*(n+1)+(1:(n+1)),]=M_t
	}
	est=est_cfar(f,1,df_b=df_b,grid)
	phihat= est$phi_coef
	psi=exp(-est$rho/n)
	mat1= matrix(0, coef_b*p, coef_b*p);
	mat2= matrix(0, coef_b*p, 1);
	phi_matinv=matrix(0, n+1, n+1)
	phi_matinv[1,1]=1;
	phi_matinv[n+1,n+1]=1;
	phi_matinv[1,2]= -psi;
	phi_matinv[n+1,n]= -psi
	for (i in 2:n){
		phi_matinv[i,i]= 1+psi^2;
		phi_matinv[i,i+1]= -psi;
		phi_matinv[i,i-1]= -psi
	}
	predict= t(matrix(indep%*%matrix(phihat,df_b+1,1),ncol=t-1,nrow=n+1))
	ssr1= sum(diag((predict)%*%phi_matinv%*% t(predict)))
	sse1= sum(diag((f[2:t,]-predict[1:(t-1),])%*%phi_matinv%*% t(f[2:t,]-predict[1:(t-1),])))
	sse0=sum(diag(f[2:t,]%*%phi_matinv%*% t(f[2:t,])))
	statistic= (sse0-sse1)/coef_b/sse1*((n+1)*(t-2)-coef_b)
	pval=1-pf(statistic,coef_b,(t-2)*(n+1)-coef_b)
	cat("Test and p-value of Order 0 vs Order 1: ","\n")
	print(c(statistic,pval))
	p=1
	sse.pre= sum(diag((f[(p+2):t,]-predict[2:(t-p),])%*%phi_matinv%*% t(f[(p+2):t,]-predict[2:(t-p),])))
	for(p in 2:p.max){
		indep.tmp=NULL
		for(i in 1:p){
			tmp=indep[((i-1)*(n+1)+1):((n+1)*(t-p-1+i)),]
			indep.tmp=cbind(tmp,indep.tmp)
		}
		indep.p=indep.tmp
		pdf4=function(para4)
		{	
			mat1= matrix(0, coef_b*p, coef_b*p);
			mat2= matrix(0, coef_b*p, 1);
			psi= para4;
			phi_matinv=matrix(0, n+1, n+1)
			phi_matinv[1,1]=1;
			phi_matinv[n+1,n+1]=1;
			phi_matinv[1,2]= -psi;
			phi_matinv[n+1,n]= -psi
			for (i in 2:n){
				phi_matinv[i,i]= 1+psi^2;
				phi_matinv[i,i+1]= -psi;
				phi_matinv[i,i-1]= -psi
			}
			psi_matinv=phi_matinv
			for(i in (p+1):t){
				tmp.n=n+1
				tmp.index=which(!is.na(f[i,]))
				M_t=indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
				tmp_mat= t(M_t) %*% psi_matinv;
				mat1= mat1+ tmp_mat %*% M_t;
				mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
			}
			phihat= solve(mat1)%*%mat2
			ehat=0
			log.mat=0
			for (i in (p+1):t){
				tmp.n=n+1
				M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
				eps= f[i,1:tmp.n]- M_t%*% phihat
				epart= t(eps)%*%psi_matinv %*%eps
				ehat= epart+ehat
				log.mat=log.mat+ log(det(psi_matinv))
			}
			l=-((t-p)*(n+1))/2*log(ehat)+1/2*log.mat
			return(-l);
		}
		para4= exp(-5)
		result4=optim(para4,pdf4, lower=0.001, upper=0.999, method='L-BFGS-B')
		mat1= matrix(0, p*coef_b, p*coef_b);
		mat2= matrix(0, coef_b*p, 1);
		psi_invall=matrix(0,(n+1)*(t-1),(n+1))
		psi=result4$par
		phi_matinv=matrix(0, n+1, n+1)
		phi_matinv[1,1]=1;
		phi_matinv[n+1,n+1]=1;
		phi_matinv[1,2]= -psi;
		phi_matinv[n+1,n]= -psi
		psi_matinv=phi_matinv
		for (i in 2:n){
			phi_matinv[i,i]= 1+psi^2;
			phi_matinv[i,i+1]= -psi;
			phi_matinv[i,i-1]= -psi
		}
		psi_matinv=phi_matinv
		for(i in (p+1):t){
			tmp.n=n+1
			M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
			tmp_mat= t(M_t) %*% psi_matinv;
			mat1= mat1+ tmp_mat %*% M_t;
			mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
		}
		phihat= solve(mat1)%*%mat2
		predict= t(matrix(indep.p%*%phihat,ncol=t-p,nrow=n+1))
		ssr.p= sum(diag((predict)%*%phi_matinv%*% t(predict)))
		sse.p= sum(diag((f[(p+1):t,]-predict)%*%phi_matinv%*% t(f[(p+1):t,]-predict)))
		statistic= (sse.pre-sse.p)/coef_b/sse.p*((n+1)*(t-p-1)-p*coef_b)
		pval=1-pf(statistic,coef_b,(t-p)*(n+1)-p*coef_b)
		cat("Test and  p-value of Order", p-1, "vs Order",p,": ","\n")
		print(c(statistic,pval))
		sse.pre= sum(diag((f[(p+2):t,]-predict[2:(t-p),])%*%phi_matinv%*% t(f[(p+2):t,]-predict[2:(t-p),])))
	}
}
#' @export
F_test_h <- function(f,weight,p.max=3,grid=1000,df_b=5,num_obs=NULL,x_pos=NULL){
	if(!is.matrix(f))f <- as.matrix(f)
	t <- nrow(f)
	if(is.null(p.max)){
		p.max=3
	}
	if(is.null(num_obs)){
		num_obs=dim(f)[2]
		n=dim(f)[2]
	}else{
		n=max(num_obs)
	}
	if(length(num_obs)!=t){
		num_obs=rep(n,t)
	}
	if(is.null(df_b)){
		df_b=10
	}
	if(is.null(grid)){
		grid=1000
	}
	if(is.null(x_pos)){
		x_pos=matrix(rep(seq(0,1,by=1/(num_obs[1]-1)),each=t),t,num_obs[1])
	}
	# parameter setting
	# grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t\
	# coef_b=6	### k=6
	# num_obs is a t by 1 vector which records N_t for each time
	x_grid=seq(0,1,by=1/grid)
	coef_b=df_b+1
	b_grid=seq(-1,1,by=1/grid)
	b_grid_sp = ns(b_grid,df=df_b);
	b_grid_full= cbind(rep(1,2*grid+1), b_grid_sp);
	indep= matrix(0, (n+1)*(t-1),coef_b)
	bstar_grid= matrix(0, grid+1, coef_b)
	b_matrix= matrix(0, (grid+1), (grid+1))
	bstar_grid_full= matrix(0, (grid+1)*t,coef_b)
	index=matrix(0,t,n)
	for(i in 1:t){
		for(j in 1:num_obs[i]){
			index[i,j]=which(abs(x_pos[i,j]-x_grid)==min(abs(x_pos[i,j]-x_grid)))
		}
	}
	for (i in 1:(t-1)){
		rec_x= approx(x_pos[i,1:num_obs[i]],f[i,1:num_obs[i]],xout=x_grid,rule=2,method='linear')	
                         # rec_x is the interpolation of x_t
		for(j in 1:coef_b){
			for (k in 1:(grid+1)){
				b_matrix[k,]= b_grid_full[seq((k+grid),k, by=-1),j]/(grid+1);
			}
			bstar_grid[, j]= b_matrix %*% matrix(rec_x$y,ncol=1,nrow=grid+1)
		}
		bstar_grid_full[(i-1)*(grid+1)+1:(grid+1),]=bstar_grid	#convoluttion of basis spline function and x_t
		tmp=bstar_grid[index[i+1,1:num_obs[i+1]],]
		indep[(i-1)*(n+1)+1:num_obs[i+1],]=tmp
	}
	#AR(1)
	pdf4=function(para4)	# Q function in formula (12)
	{	
		mat1= matrix(0, df_b+1, df_b+1);
		mat2= matrix(0, df_b+1, 1);
		### correlation matrix of error process at observation points
		psi_invall=matrix(0,(n+1)*(t-1),(n+1))
		for(i in 2:t){
			psi= para4;
			tmp.n=num_obs[i]
			psi_mat= matrix(1, tmp.n, tmp.n)
			for(k in 1:(tmp.n)){
				for(j in 1:(tmp.n)){
					psi_mat[k,j]= psi^(abs(x_pos[i,j]-x_pos[i,k]))
				}
			}
			psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
			psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]=psi_matinv2
			M_t=indep[(i-2)*(n+1)+(1:tmp.n),]
			tmp_mat= t(M_t) %*% psi_matinv2;
			mat1= mat1+ tmp_mat %*% M_t;			# matrix under the first brackets in formula (15)
			mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];	# matrix under the second brackets in formula (15)
		}
		phihat= solve(mat1)%*%mat2
		ehat=0	# e(beta,phi) in formula (12)
		log.mat=0
		for (i in 2:t){
			tmp.n=num_obs[i]
			M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
			eps= f[i,1:tmp.n]- M_t%*% phihat
			psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
			epart= t(eps)%*%psi_matinv %*%eps
			ehat= epart+ehat
			log.mat=log.mat+ log(det(psi_matinv))
		}
		l=-sum(num_obs[2:t])/2*log(ehat)+1/2*log.mat	# the number defined in formula (14)
		return(-l);
	}
	para4= exp(-1)
	result4=optim(para4,pdf4, lower=0.001, upper=0.999, method='L-BFGS-B')
	mat1= matrix(0, df_b+1, df_b+1);
	mat2= matrix(0, df_b+1, 1);
	psi_invall=matrix(0,(n+1)*(t-1),(n+1))
	psi=result4$par
	for(i in 2:t){
		tmp.n=num_obs[i]
		psi_mat= matrix(1, tmp.n, tmp.n)
		for(k in 1:tmp.n){
			for(j in 1:tmp.n){
				psi_mat[k,j]=psi^(abs(x_pos[i,j]-x_pos[i,k]))
			}
		}
		psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))	
		psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]= psi_matinv2
		M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
		tmp_mat= t(M_t) %*% psi_matinv2;
		mat1= mat1+ tmp_mat %*% M_t;
		mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
	}
	phihat= solve(mat1)%*%mat2	
	ehat=0		# e(beta,phi) in formula (12)
	for (i in 2:t){
		tmp.n=num_obs[i]
		M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
		eps= f[i,1:tmp.n]- M_t%*% phihat
		psi_matinv=psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
		epart= t(eps)%*%psi_matinv %*%eps
		ehat= epart+ehat
	}
	rho_hat= -log(psi)
	sigma_hat=ehat/sum(num_obs[2:t])*2*rho_hat

	sigma_hat1=sigma_hat
	rho_hat1=rho_hat
	predict1= t(matrix(indep%*%phihat,ncol=t-1,nrow=n+1))
	ssr1=0
	sse1=0
	for (i in 2:t){
		tmp.n=num_obs[i]
		psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
		ssr1= ssr1+ sum((predict1[i-1,1:tmp.n] %*% psi_matinv) * (predict1[i-1,1:tmp.n]))
		sse1= sse1+ sum(((f[i,1:tmp.n]- predict1[i-1,1:tmp.n]) %*% psi_matinv) *(f[i,1:tmp.n]-predict1[i-1,1:tmp.n]))
	}
	#   No AR term
	psi=exp(-rho_hat1)
	mat1= matrix(0, df_b+1, df_b+1);
	mat2= matrix(0, df_b+1, 1);
	phihat=matrix(0,coef_b,1)
	ehat=0
	for (i in 2:t){
		tmp.n=num_obs[i]
		M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
		eps= f[i,1:tmp.n]- M_t%*% phihat
		psi_matinv=psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
		epart= t(eps)%*%psi_matinv %*%eps
		ehat= epart+ehat
	}
	rho_hat= -log(psi)
	sigma_hat=ehat/sum(num_obs[2:t])*2*rho_hat
	rho_hat0=rho_hat
	sigma_hat0=sigma_hat
	sse0=0
	test=0
	for (i in 2:t){
		tmp.n=num_obs[i]
		psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
		sse0= sse0+ sum((f[i,1:tmp.n] %*% psi_matinv) *(f[i,1:tmp.n]))
	}
	statistic= (sse0-sse1)/coef_b/sse1*(sum(num_obs[2:t])-coef_b)
	pval=1-pf(statistic,coef_b,sum(num_obs[2:t])-coef_b)
	cat("Test and  p-value of Order 0 vs Order 1: ","\n")
	print(c(statistic,pval))
	indep.p=indep
	if (p.max>1){
		for(p in 2:p.max){
			indep.pre=indep.p[(n+2):((n+1)*(t-p+1)),]
			indep.tmp=indep
			for (i in 1:(t-p)){
				for(j in 1:num_obs[i+p]){
					index[i+p,j]= which(abs(x_pos[i+p,j]-x_grid)==min(abs(x_pos[i+p,j]-x_grid)))
				}
				tmp=bstar_grid_full[(i-1)*(grid+1)+index[i+p,1:num_obs[i+p]],]
				indep.tmp[(i-1)*(n+1)+1:num_obs[i+p],]=tmp
			}
			indep.test=cbind(indep.p[(n+2):((n+1)*(t-p+1)),], indep.tmp[1:((n+1)*(t-p)),])
			indep.p=indep.test

			pdf4=function(para4)	
			{	
				mat1= matrix(0, p*(df_b+1), p*(df_b+1));
				mat2= matrix(0, p*(df_b+1), 1);
				### correlation matrix of error process at observation points
				psi_invall=matrix(0,(n+1)*(t-1),(n+1))
				for(i in (p+1):t){
					psi= para4;
					tmp.n=num_obs[i]
					psi_mat= matrix(1, tmp.n, tmp.n)
					for(k in 1:(tmp.n)){
						for(j in 1:(tmp.n)){
							psi_mat[k,j]= psi^(abs(x_pos[i,j]-x_pos[i,k]))
						}
					}
					psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
					psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]=psi_matinv2
					M_t=indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
					tmp_mat= t(M_t) %*% psi_matinv2;
					mat1= mat1+ tmp_mat %*% M_t;		
					mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
				}
				phihat= solve(mat1)%*%mat2
				ehat=0	
				log.mat=0
				for (i in (p+1):t){
					tmp.n=num_obs[i]
					M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]	
					eps= f[i,1:tmp.n]- M_t%*% phihat
					psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
					epart= t(eps)%*%psi_matinv %*%eps
					ehat= epart+ehat
					log.mat=log.mat+ log(det(psi_matinv))
				}
				l=-sum(num_obs[(p+1):t])/2*log(ehat)+1/2*log.mat	### the number defined in formula (14)
				return(-l);
			}
			para4= exp(-1)
			result4=optim(para4,pdf4, lower=0.001, upper=0.999, method='L-BFGS-B')
			mat1= matrix(0, p*(df_b+1), p*(df_b+1));
			mat2= matrix(0, p*(df_b+1), 1);
			psi_invall=matrix(0,(n+1)*(t-1),(n+1))
			psi=result4$par
			for(i in (p+1):t){
				tmp.n=num_obs[i]
				psi_mat= matrix(1, tmp.n, tmp.n)
				for(k in 1:tmp.n){
					for(j in 1:tmp.n){
						psi_mat[k,j]=psi^(abs(x_pos[i,j]-x_pos[i,k]))
					}
				}
				psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))	
				psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]= psi_matinv2
				M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
				tmp_mat= t(M_t) %*% psi_matinv2;
				mat1= mat1+ tmp_mat %*% M_t;
				mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
			}
			phihat= solve(mat1)%*%mat2	
			ehat=0		
			for (i in (p+1):t){
				tmp.n=num_obs[i]
				M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
				eps= f[i,1:tmp.n]- M_t%*% phihat
				psi_matinv=psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
				epart= t(eps)%*%psi_matinv %*%eps
				ehat= epart+ehat
			}
			rho_hat= -log(psi)
			sigma_hat=ehat/sum(num_obs[(p+1):t])*2*rho_hat

			sigma_hat1=sigma_hat
			rho_hat1=rho_hat
			predict= t(matrix(indep.p%*%phihat,ncol=t-p,nrow=n+1))
			ssr.p= 0
			sse.p=0
			for (i in (p+1):t){
				tmp.n=num_obs[i]
				psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
				ssr.p= ssr.p+ sum((predict[i-p,1:tmp.n] %*% psi_matinv) * (predict[i-p,1:tmp.n]))
				sse.p= sse.p+ sum(((f[i,1:tmp.n]- predict[i-p,1:tmp.n]) %*% psi_matinv) *(f[i,1:tmp.n]-predict[i-p,1:tmp.n]))
			}

			mat1= matrix(0, coef_b*(p-1), coef_b*(p-1));
			mat2= matrix(0, coef_b*(p-1), 1);
			psi_invall=matrix(0,(n+1)*(t-1),(n+1))
			for(i in (p+1):t){
				tmp.n=num_obs[i]
				tmp.index=which(!is.na(f[i,]))
				psi_mat= matrix(1, tmp.n, tmp.n)
				for(k in 1:tmp.n){
					for(j in 1:tmp.n){
						psi_mat[k,j]=psi^(abs(x_pos[i,j]-x_pos[i,k]))
					}
				}
				psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
				psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]= psi_matinv2
				M_t= indep.pre[(i-p-1)*(n+1)+(1:tmp.n),]
				tmp_mat= t(M_t) %*% psi_matinv2;
				mat1= mat1+ tmp_mat %*% M_t;
				mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
			}
			phihat= solve(mat1)%*%mat2
			predict= t(matrix(indep.pre%*%phihat,ncol=t-p,nrow=n+1))
			ssr.pre=0
			sse.pre=0.0000
			for (i in (p+1):t){
				tmp.n=num_obs[i]
				tmp.index=which(!is.na(f[i,]))
				psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
				ssr.pre= ssr.pre+ sum((predict[i-p,1:tmp.n] %*% psi_matinv) * (predict[i-p,1:tmp.n]))
				sse.pre= sse.pre+ sum(((f[i,1:tmp.n]- predict[i-p,1:tmp.n]) %*% psi_matinv) *(f[i,1:tmp.n]-predict[i-p,1:tmp.n]))
			}
			statistic= (sse.pre-sse.p)/coef_b/sse.p*(sum(num_obs[(p+1):t])-coef_b)
			pval=1-pf(statistic,coef_b,sum(num_obs[(p+1):t])-coef_b)
			cat("Test and  p-value of Order",p-1," vs Order", p,": ","\n")
			print(c(statistic,pval))	
		}
	}
}


