#' @export
est_cfar <- function(f,p=3,df_b=10,grid=1000){
	if(!is.matrix(f))f <- as.matrix(f)
	t <- nrow(f)
	n=dim(f)[2]-1
	if(is.null(grid)){
		grid=1000
	}
	if(is.null(df_b)){
		df_b=10
	}
	#INPUTS
	# parameter setting
	# rho=5	### parameter for error process epsilon_t, O-U process
	# t= 1000;	### length of time
	# iter= 100;	### number of replications in the simulation
	# grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t
	# df_b=10		### number of degrees for splines
	# n=100		### number of observations for each X_t, in the paper it is 'N'
	x_grid= seq(0, 1, by= 1/grid);	# the grid points for input variables of functional time series in [0,1]
	b_grid= seq(-1, 1, by=1/grid);	# the grid points for input variables of b-spline functions in [-1,1]
	eps_grid= matrix(0, 1, grid+1);	# the grid points for input variables of error process in [0,1]
	# It shows the best approximation of phi using b-spline with df=df_b
	coef_b=df_b+1;
	b_grid_sp= ns(b_grid, df=df_b);

	# Estimation
	x= seq(0, 1, by=1/n);
	index= 1:(n+1)*(grid/n)-(grid/n)+1
	x_grid_sp= bs(x_grid, df=n, degree=1);	# basis function values if we interpolate discrete observations of x_t 
	x_sp= x_grid_sp[index,];			# basis function values at each observation point 
	df=n;
	coef= df+1
	# convolutions of b-pline functions (df=df) for phi_1 and phi_2 and basis functions (df=n) for x
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
	indep=matrix(0,(n+1)*(t-1),coef_b)	# design matrix for b-spline coefficient of phi_1
	dsg=cbind(rep(1,n+1),x_sp);
	dsg_mat= solve(t(dsg)%*%dsg)%*%t(dsg)
	ahat= f%*%t(dsg_mat)	# b-spline coefficient when we interpolate x_t
	for (i in 2:t){	
		M_t= matrix(ahat[i-1,]%*% t(bstar3), nrow=n+1, ncol=df_b+1)
		indep[(i-2)*(n+1)+(1:(n+1)),]=M_t	# design matrix for b-spline coefficient of phi_1
	}
	indep.tmp=NULL
	for(i in 1:p){
		tmp=indep[((i-1)*(n+1)+1):((n+1)*(t-p-1+i)),]
		indep.tmp=cbind(tmp,indep.tmp)
	}
	indep.p=indep.tmp
	pdf4=function(para4)
	{	
		psi= para4;	# correlation between two adjacent observation point exp(-rho/n)
		psi_matinv= matrix(0, n+1, n+1)	# correlation matrix of error process at observation points
		psi_matinv[1,1]=1
		psi_matinv[n+1,n+1]=1;
		psi_matinv[1,2]=-psi;
		psi_matinv[n+1,n]=-psi;
		for (i in 2:n){
			psi_matinv[i,i]= 1+psi^2;
			psi_matinv[i,i+1]= -psi;
			psi_matinv[i,i-1]= -psi;
		}
		mat1= matrix(0, coef_b*p, coef_b*p);
		mat2= matrix(0, coef_b*p, 1)
		for(i in (p+1):t){
			tmp.n=n+1
			M_t=indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
			tmp_mat= t(M_t) %*% psi_matinv;
			mat1= mat1+ tmp_mat %*% M_t;
			mat2= mat2+ tmp_mat %*% f[i,];
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
		l=-(t-p)*n/2*log(ehat)+1/2*log.mat
		return(-l);
	}
	para4= exp(-5)
	result4=optim(para4,pdf4, lower=0.001, upper=0.9999, method='L-BFGS-B')
	psi=result4$par
	psi_matinv= matrix(0, n+1, n+1)
	psi_matinv[1,1]=1
	psi_matinv[n+1,n+1]=1;
	psi_matinv[1,2]=-psi;
	psi_matinv[n+1,n]=-psi;
	for (i in 2:n){
		psi_matinv[i,i]= 1+psi^2;
		psi_matinv[i,i+1]= -psi;
		psi_matinv[i,i-1]= -psi;
	}
	mat1= matrix(0, p*(df_b+1), p*(df_b+1));	# matrix under the first brackets in formula (15)
	mat2= matrix(0, p*(df_b+1), 1);	# matrix under the second brackets in formula (15)
	for (i in (p+1):t){	
		M_t= indep.p[(i-p-1)*(n+1)+(1:(n+1)),]
		tmp_mat= t(M_t) %*% psi_matinv;
		mat1= mat1+ tmp_mat %*% M_t;
		mat2= mat2+ tmp_mat %*% f[i,];
	}
	phihat= solve(mat1)%*%mat2
	phihat_func=matrix(0,p,(1+2*grid))
	for(i in 1:p){
		phihat_func[i,]= cbind(rep(1,2*grid+1),b_grid_sp)%*%phihat[(1+(coef_b*(i-1))):(i*coef_b)];		###b-spline coefficient of estimated phi_1
	}
	ehat=0
	for (i in (p+1):t){
		tmp.n= n+1
		M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
		eps= f[i,1:tmp.n]- M_t%*% phihat
		epart= t(eps)%*%psi_matinv %*%eps
		ehat= epart+ehat
	}
	rho_hat= -log(psi)*n
	sigma_hat=ehat/((t-p)*(n+1))*2*rho_hat/(1-psi^2)
	phi_coef=t(matrix(phihat,coef_b,p))
	est_cfar <- list(phi_coef=phi_coef,phi_func=phihat_func,rho=rho_hat,sigma2=sigma_hat)
	return(est_cfar)
}
#' @export
est_cfarh <- function(f,weight,p=2,grid=1000,df_b=5, num_obs=NULL,x_pos=NULL){
	# CFAR(2) with processes with heteroscedasticity, irregular observation locations
	# Estimation of phi(), rho and sigma, and F test 
	if(!is.matrix(f))f <- as.matrix(f)
	t <- nrow(f)
	if(is.null(num_obs)){
		num_obs=dim(f)[2]
		n=dim(f)[2]
	}else{
		n=max(num_obs)
	}
	if(length(num_obs)!=t){
		num_obs=rep(n,t)
	}
	if(is.null(x_pos)){
		x_pos=matrix(rep(seq(0,1,by=1/(num_obs[1]-1)),each=t),t,num_obs[1])
	}
	if(is.null(df_b)){
		df_b=10
	}
	if(is.null(grid)){
		grid=1000
	}
	# parameter setting
	# rho=1	### parameter for error process epsilon_t, O-U process
	# tmax= 1001;	### maximum length of time to generated data
	# t=1000	### length of time
	# iter= 100;	### number of replications in the simulation
	# grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t\
	# coef_b=6	### k=6
	# min_obs=40	### the minimal number of observations at each time
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
                          
		for(j in 1:coef_b){
			for (k in 1:(grid+1)){
				b_matrix[k,]= b_grid_full[seq((k+grid),k, by=-1),j]/(grid+1);
			}
			bstar_grid[, j]= b_matrix %*% matrix(rec_x$y,ncol=1,nrow=grid+1)
		}
		bstar_grid_full[(i-1)*(grid+1)+1:(grid+1),]=bstar_grid	
		tmp=bstar_grid[index[i+1,1:num_obs[i+1]],]
		indep[(i-1)*(n+1)+1:num_obs[i+1],]=tmp
	}
	if(p==1){
		# AR(1)
		pdf4=function(para4)	
		{	
			mat1= matrix(0, df_b+1, df_b+1);
			mat2= matrix(0, df_b+1, 1);
			# correlation matrix of error process at observation points
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
				mat1= mat1+ tmp_mat %*% M_t;		
				mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
			}
			phihat= solve(mat1)%*%mat2
			ehat=0	
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
			l=-sum(num_obs[2:t])/2*log(ehat)+1/2*log.mat	
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
		phihat_func=matrix(0,p,(1+2*grid))
		for(i in 1:p){
			phihat_func[i,]= cbind(rep(1,2*grid+1),b_grid_sp)%*%phihat[(1+(coef_b*(i-1))):(i*coef_b)];		
		}
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

	}
	indep.p=indep	
	if(p>1){
		for(q in 2:p){
			indep.tmp=indep
			for (i in 1:(t-q)){
				for(j in 1:num_obs[i+q]){
					index[i+q,j]= which(abs(x_pos[i+q,j]-x_grid)==min(abs(x_pos[i+q,j]-x_grid)))
				}
				tmp=bstar_grid_full[(i-1)*(grid+1)+index[i+q,1:num_obs[i+q]],]
				indep.tmp[(i-1)*(n+1)+1:num_obs[i+q],]=tmp
			}
			indep.test=cbind(indep.p[(n+2):((n+1)*(t-q+1)),], indep.tmp[1:((n+1)*(t-q)),])
			indep.p=indep.test

		}
	
		pdf4=function(para4)	
		{	
			mat1= matrix(0, p*(df_b+1), p*(df_b+1));
			mat2= matrix(0, p*(df_b+1), 1);
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
			l=-sum(num_obs[(p+1):t])/2*log(ehat)+1/2*log.mat	
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
		phihat_func=matrix(0,p,(1+2*grid))
		for(i in 1:p){
			phihat_func[i,]= cbind(rep(1,2*grid+1),b_grid_sp)%*%phihat[(1+(coef_b*(i-1))):(i*coef_b)];		
		}	
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
	}
	est_cfarh <- list(phi_coef=t(matrix(phihat,df_b+1,p)),phi_func=phihat_func, rho=rho_hat, sigma2=sigma_hat)
	return(est_cfarh)
}



