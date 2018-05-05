#' Prediction
#' @export
p_cfar<-function(model, f, m=3){
	# p_cfar: this function gives us the prediction of functional time series. There are 3 input variables:
	# 1. model is the estimated model obtained from est_cfar function
	# 2. f is the matrix which contains the data
	# 3. m is the forecasting horizon
	# parameter setting
	# t= 100;length of time
	# grid=1000; the number of grid points used to construct the functional time series X_t and epsilon_t
	# p is the cfar order
	# n=100, number of observations for each X_t, in the paper it is 'N'
	if(!is.matrix(f))f = as.matrix(f)
	# Estimation
	n=dim(f)[2]-1
	t=dim(f)[1]
	p=dim(model$phi_coef)[1]
	grid=(dim(model$phi_func)[2]-1)/2
	pred=matrix(0,t+m,grid+1)
	x_grid=seq(0,1,by=1/grid)
	x=seq(0,1,by=1/n)
	for(i in (t-p):t){
		pred[i,]= approx(x,f[i,],xout=x_grid,rule=2,method='linear')$y
	}
	phihat_func=model$phi_func
	for (j in 1:m){
		for (i in 1:(grid+1)){
			pred[j+t,i]= sum(phihat_func[1:p,seq((i+grid),i, by=-1)]*pred[(j+t-1):(j+t-p),])/grid
		}
	}
	pred_cfar=pred[(t+1):(t+m),]
	return(pred_cfar)
}
#' partial prediction
#' Assume that we have t curves and want to predit the curve at time t+1, but we know the first n observations in the curve, to predict the n+1 observation.
#' @export
p_cfar_part<-function(model, f, new.obs){
	# p-cfar_part: this function gives us the prediction of x_t(s), s>s_0, give x_1,\ldots, x_{t-1}, and x(s), s <s_0. There are three input variables.
	# 1. model is the estimated model obtained from est_cfar function
	# 2. f is the matrix which contains the data
	# 3. new.obs is x_t(s), s>s_0.
	t=dim(f)[1]
	p=dim(model$phi_coef)[1]
	f_grid=matrix(0,t,grid+1)
	grid=(dim(model$phi_func)[2]-1)/2
	x_grid=seq(0,1,by=1/grid)
	n=dim(f)[2]-1
	index=seq(1,grid+1,by=grid/n)
	x=seq(0,1,by=1/n)
	for(i in (t-p):t){
		f_grid[i,]=approx(x,f[i,],xout=x_grid,rule=2,method='linear')$y
	}
	pred=matrix(0,grid+1,1)
	phihat_func=model$phi_func
	for (i in 1:(grid+1)){
		pred[i]= sum(phihat_func[1:p,seq((i+grid),i, by=-1)]*f_grid[t:(t-p+1),])/grid
	}
	pred_new=matrix(0,n+1,0)
	pred_new[1]=pred[1]
	rho=model$rho
	pred_new[2:(n+1)]=pred[index[1:n]]+(new.obs[1:n]-pred[index[1:n]])*exp(-rho/n)
	return(pred_new)
}

