#' CFAR program
#' based on the original codes of Xialu Liu of SDSU.
#' @export
g_cfar1 <- function(tmax=1001,rho=5,phi_func=NULL,grid=1000,sigma=1){
	#Simulate CFAR(1) processes
	#parameter setting
	#rho parameter for error process epsilon_t, O-U process
	#sigma is the standard deviation of OU-process
	#tmax length of time
	#iter number of replications in the simulation, iter=1
	#grid the number of grid points used to construct the functional time series X_t and epsilon_t
        # phi_func:   User supplied convolution function phi(.). The following default is used if not specified.
	if(is.null(phi_func)){
		phi_func= function(x){
			return(dnorm(x,mean=0,sd=0.1))
            }
    	}
	if(is.null(grid)){
		grid=1000
	}
	if(is.null(sigma)){
		sigma=1
	}
	# OUTPUTS
	# A matrix f_grid is generated, with (tmax) rows and (grid+1) columns. 
	x_grid<- seq(0, 1, by= 1/grid); # the grid points for input variables of functional time series in [0,1]
	b_grid<- seq(-1, 1, by=1/grid);	# the grid points for input variables of b-spline functions in [-1,1]
	eps_grid <- matrix(0, 1, grid+1);	# the grid points for input variables of error process in [0,1]
	fi <- exp(-rho/grid);		# the correlation of two adjacent points for an O-U process
	phi_b <- phi_func(b_grid)		# phi(s) when s in [-1,1]
	#the convolutional function values phi(s-u), for different s in [0,1]
	phi_matrix <- matrix(0, (grid+1), (grid+1));
	for (i in 1:(grid+1)){
		phi_matrix[i,]= phi_b[seq((i+grid),i,by=-1)]/(grid+1)
	}
	# Generate data
	# functional time series f_grid is X_t in the paper
	# A matrix f_grid is generated, with (tmax*iter) rows and (grid+1) columns. 
	# It contains iter CFAR(1) processes. The i-th process is in rows (i-1)*tmax+1:tmax.
	f_grid=matrix(0,tmax,grid+1)
	i=1
	eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt((1-fi^2))*rnorm(n,0,1))  ###error process
	f_grid[1,]= eps_grid;	
	for (i in 2:tmax){
		eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
		f_grid[i,]= phi_matrix%*%f_grid[i-1,]+ eps_grid;
	}
	g_cfar1 <- list(cfar1=f_grid*sigma,epsilon=eps_grid*sigma)
	return(g_cfar1)
}
#' @export
g_cfar2 <- function(tmax=1001,rho=5,phi_func1=NULL, phi_func2=NULL,grid=1000,sigma=1){
	#Simulate CFAR(2) processes
	# OUTPUTS
	# A matrix f_grid is generated, with (tmax) rows and (grid+1) columns.
	#parameter setting
	# rho=5 parameter for error process epsilon_t, O-U process
	#tmax= 1001; length of time
	#iter= 100; number of replications in the simulation
	#grid=1000; the number of grid points used to construct the functional time series X_t and epsilon_t
        # phi_func1 and phi_func2: User specified convolution functions. The following defaults are used if not given
	if(is.null(phi_func1)){
		phi_func1= function(x){
			return(0.5*x^2+0.5*x+0.13)
    		}
	}
	if(is.null(phi_func2)){
	    phi_func2= function(x){
		  return(0.7*x^4-0.1*x^3-0.15*x)
	    }
	}
	if(is.null(grid)){
		grid=1000
	}
	if(is.null(sigma)){
		sigma=1
	}
	x_grid= seq(0, 1, by= 1/grid);	# the grid points for input variables of functional time series in [0,1]
	b_grid= seq(-1, 1, by=1/grid);	# the grid points for input variables of b-spline functions in [-1,1]
	eps_grid= matrix(0, 1, grid+1);	# the grid points for input variables of error process in [0,1]
	fi = exp(-rho/grid);		# the correlation of two adjacent points for an O-U process
	phi_b1= phi_func1(b_grid)	# phi_1(s) when s in [-1,1]
	phi_b2= phi_func2(b_grid)	# phi_2(s) when s in [-1,1]

	# the convolutional function values phi_1(s-u) and phi_2(s-u), for different s in [0,1]
	phi_matrix1= matrix(0, (grid+1), (grid+1));
	for (i in 1:(grid+1)){
		phi_matrix1[i,]= phi_b1[seq((i+grid),i,by=-1)]/(grid+1)
	}
	phi_matrix2= matrix(0, (grid+1), (grid+1));
	for (i in 1:(grid+1)){
		phi_matrix2[i,]= phi_b2[seq((i+grid),i,by=-1)]/(grid+1)
	}
	# Generate data
	# functional time series f_grid is X_t in the paper
	# A matrix f_grid is generated, with (tmax*iter) rows and (grid+1) columns. 
	#  It contains iter CFAR(2) processes. The i-th process is in rows (i-1)*tmax+1:tmax.
	f_grid=matrix(0,tmax,grid+1)
	i=1
	eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
	f_grid[1,]= eps_grid;
	i=2
	eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
	f_grid[i,]= phi_matrix1%*%as.matrix(f_grid[i-1,])+ eps_grid;
	for (i in 3:tmax){
		eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
		f_grid[i,]= phi_matrix1%*%as.matrix(f_grid[i-1,])+ phi_matrix2%*%as.matrix(f_grid[i-2,])+ eps_grid;
	}
	# The last iteration of epsilon_t is returned.
	g_cfar2 <- list(cfar2=f_grid*sigma,epsilon=eps_grid*sigma)
}
#' @export
g_cfar <- function(tmax=1001,rho=5,phi_list=NULL, grid=1000,sigma=1){
	#Simulate CFAR(p) processes
	# OUTPUTS
	# A matrix f_grid is generated, with (tmax) rows and (grid+1) columns.
	# parameter setting
	# rho=5 parameter for error process epsilon_t, O-U process
	# tmax= 1001;	# length of time
	# iter= 100;	# number of replications in the simulation
	#grid=1000;	# the number of grid points used to construct the functional time series X_t and epsilon_t
       # phi_func User specified convolution functions. The following defaults are used if not given
	if(is.null(phi_list)){
		phi_func1= function(x){
			return(dnorm(x,mean=0,sd=0.1))
    		}
		phi_list=phi_func1
	}
	if(is.null(grid)){
		grid=1000
	}
	if(is.null(sigma)){
		sigma=1
	}
	x_grid= seq(0, 1, by= 1/grid);	
	b_grid= seq(-1, 1, by=1/grid);	
	eps_grid= matrix(0, 1, grid+1);	
	fi = exp(-rho/grid);		
	p=length(phi_list)
	if(is.list(phi_list)){
		phi_b=sapply(phi_list,mapply,b_grid)	
	}else{
		phi_b=matrix(phi_list(b_grid),2*grid+1,1)
	}
	# the convolutional function values phi_1(s-u) and phi_2(s-u), for different s in [0,1]
	phi_matrix= matrix(0, (grid+1), p*(grid+1));
	for(j in 1:p){
		for (i in 1:(grid+1)){
			phi_matrix[i,((j-1)*(grid+1)+1):(j*(grid+1))]= phi_b[seq((i+grid),i,by=-1),j]/(grid+1)
		}
	}
	# Generate data
	# functional time series f_grid is X_t in the paper
	# A matrix f_grid is generated, with (tmax*iter) rows and (grid+1) columns. 
	#  It contains iter CFAR(2) processes. The i-th process is in rows (i-1)*tmax+1:tmax.
	f_grid=matrix(0,tmax,grid+1)
	i=1
	eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
	f_grid[1,]= eps_grid;
	if(p>1){
		for(i in 2:p){
			eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
			f_grid[i,]= phi_matrix[,1:((i-1)*(grid+1))]%*%matrix(t(f_grid[(i-1):1,]),(i-1)*(grid+1),1)+ eps_grid;
		}
	}
	for (i in (p+1):tmax){
		eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
		f_grid[i,]= phi_matrix%*%matrix(t(f_grid[(i-p):(i-1),]),p*(grid+1),1)+ eps_grid;
	}
	# The last iteration of epsilon_t is returned.
	g_cfar <- list(cfar=f_grid*sigma,epsilon=eps_grid*sigma)
}
#' @export
g_cfar2h <- function(tmax=1001,grid=1000,rho=1,min_obs=40, pois=5,phi_func1=NULL, phi_func2=NULL, weight=NULL){
	#Simulate CFAR(2) processes with heteroscedasticity, irregular observation locations
	# We need to have f_grid to record the functional time series
	# We need to have x_pos_full to record the observation positions
	# parameter setting
	# rho=1	### parameter for error process epsilon_t, O-U process
	# tmax= 1001;	### length of time
	# iter= 100;	### number of replications in the simulation
	# grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t\
	# min_obs=40	### the minimal number of observations at each time
      # phi_func1, phi_func2, weight: User specified convolution functions and weight function. The following defaults 
      #  are used if not given
	if(is.null(phi_func1)){
		phi_func1= function(x){
			return(0.5*x^2+0.5*x+0.13)
    		}
	}
	if(is.null(phi_func2)){
		phi_func2= function(x){
			return(0.7*x^4-0.1*x^3-0.15*x)
        	}
	}
	if(is.null(weight)){
		###heteroscedasticity weight function
		x_grid <- seq(0,1,by=1/grid)
		weight0= function(w){
			return((w<=0.6)*exp(-10*w)+(w>0.6)*(exp(-6)+0.2*(w-0.6))+0.1);
		}
		const=sum(weight0(x_grid)/(grid+1))
		weight= function(w){
			return(((w<=0.6)*exp(-10*w)+(w>0.6)*(exp(-6)+0.2*(w-0.6))+0.1)/const)
 		}
	}
	num_full=rpois(tmax,pois)+min_obs	# number of observations at time t follows a Poisson distribution 
	x_grid= seq(0, 1, by= 1/grid);	# the grid points for input variables of functional time series in [0,1]
	b_grid= seq(-1, 1, by=1/grid);	# the grid points for input variables of b-spline functions in [-1,1]
	eps_grid= matrix(0, 1, grid+1);	# the grid points for input variables of error process in [0,1]
	fi = exp(-rho/grid);		# the correlation of two adjacent points for an O-U process
	phi_b1= phi_func1(b_grid)	# phi_1(s) when s in [-1,1]
	phi_b2= phi_func2(b_grid)	# phi_2(s) when s in [-1,1]
	# the convolutional function values phi_1(s-u) and phi_2(s-u), for different s in [0,1]
	phi_matrix1= matrix(0, (grid+1), (grid+1));
	for (i in 1:(grid+1)){
		phi_matrix1[i,]= phi_b1[seq((i+grid),i,by=-1)]/(grid+1)
	}
	phi_matrix2= matrix(0, (grid+1), (grid+1));
	for (i in 1:(grid+1)){
		phi_matrix2[i,]= phi_b2[seq((i+grid),i,by=-1)]/(grid+1)	
	}
	n=max(num_full)
	x_pos_full=matrix(0,tmax,n)	# observation positions
	for(k in 1:(tmax)){
		x_pos_full[k,1:num_full[k]]=sort(runif(num_full[k]))
	}
	f_grid=matrix(0,tmax,grid+1)
	i=1
	eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
	f_grid[1,]= eps_grid*weight(x_grid);
	i=2
	f_grid[i,]= phi_matrix1 %*% as.matrix(f_grid[i-1,])+ eps_grid*weight(x_grid);	
	for (i in 3:tmax){
		eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
		f_grid[i,]= phi_matrix1%*%as.matrix(f_grid[i-1,]) + phi_matrix2%*%as.matrix(f_grid[i-2,])+ eps_grid*weight(x_grid);
	}

	# Functional time series, number of observations at different time of periods, 
	# observation points at different time of periods, the last iteration of epsilon_t is returned.
	g_cfar2h <- list(cfar2h=f_grid, num_obs=num_full, x_pos=x_pos_full,epsilon=eps_grid)
}


