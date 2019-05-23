//
#include <cmath>
#include <array>
#include <tuple>
#include "vector_maths.h"
#include "constants.h"
using namespace std;

array<array<long double, 3>, N*2> evolve(const int N, array<long double,N> M, array<array<long double, 3>, N*2> X, long double T, int n=1, long double eta = 0.02)
{
	//Current time
	long double t = 0.0;
	long double dt_max, dt;
	tuple<array<array<long double, 3>, N>, long double> result;
	while (t<T){
		dt_max = T - t;
		result = singleTimestep(N, X, M, n, eta, dt_max=dt_max);
		X = get<0>(result);
		dt = get<1>(result);
		t += dt;
	}
	return X;
}

//Evolution of binary through integration
tuple<array<>, long double> singleTimestep(const int N, array<array<long double, 3>, N*2> X, array<long double, N> M, int n, long double eta, long double dt_max=-1.0)
{
    //dt is step size
    //X = [x1, x2, x3,..., v1, v2, v3,...]
    //N is number of objects
    //M = [m1, m2, m3,...]
    //n is number of iterations of estimate and correct
	array<array<long double,3>, N> X_0, V_0, A_0, J_0, X_1, V_1, A_1, J_1;
	long double dt;
	for (int i=0; i<N; ++i){
		for (int j=0; j<3; ++j){
			X_0[i][j] = X[i][j];
			V_0[i][j] = X[i+N][j];
		}
	} 
    //Hermite scheme (Dehnen and Read 2011) with Aarseth criterion timesteps  
    //Find initial acceleration and jerk and timestep
    ajt = acc_jerk_and_timestep(N, X_0, V_0, M, eta); 
    A_0 = get<0>(ajt);
    J_0 = get<1>(ajt);
    dt = get<2>(ajt);
    if ((dt_max>0.0) && (dt_max<dt)){
    	dt = dt_max;
    }
    //Predict positions and velocities
    for (int i=0; i<N; ++i){
    	for (int j=0; j<3; ++j){
    		X_1[i][j] = X_0[i][j] + V_0[i][j]*dt + A_0[i][j]*dt*dt/2.0 + J_0[i][j]*pow(dt,3.0)/6.0;
    		V_1[i][j] = V_0[i][j] + A_0[i][j]*dt + J_0[i][j]*dt*dt/2.0;
    	}
    }
    //Iterate n times
    for (int k=0; k<n; ++k){   
        //Estimate acceleration and jerk
        aj = acc_and_jerk(N, X_1, V_1, M);
        A_1 = get<0>(aj);
        J_1 = get<1>(aj);     
        //Obtain corrected positions and velocities
        for (int i=0; i<N; ++i){
        	for (int j=0; j<3; ++j){
        		X_1[i][j] = X_0[i][j] + (V_1[i][j]+V_0[i][j])*dt/2.0 + (A_0[i][j]-A_1[i][j])*dt*dt/12.0;
	        	V_1[i][j] = V_0[i][j] + (A_1[i][j]+A_0[i][j])*dt/2.0 + (J_0[i][j]-J_1[i][j])*dt*dt/12.0;
        	}
        }
    }
    for (int i=0; i<N; ++i){
		for (int j=0; j<3; ++j){
			X[i][j] = X_1[i][j];
			X[i+N][j] = V_1[i][j];
		}
	} 
    return make_tuple(X, dt);
}

//See Nitadori and Makino 2008 for higher derivatives of acceleration
array<array<long double,N>, N> calc_alpha(const int N, array<array<array<long double,3>, N>, N> x, array<array<array<long double,3>, N>, N> v){
	array<array<long double,N>, N> alpha;
	for (int i; i<N; ++i){
		for (int j=0; j<i; ++j){
			alpha[i][j] = dot(x[i][j], v[i][j])/dot(x[i][j], x[i][j]);
		}
		alpha[i][i] = 0.0;
		for (int j=i+1; j<N; ++j){
			alpha[i][j] = dot(x[i][j], v[i][j])/dot(x[i][j], x[i][j]);
		}
	}
	return alpha;
}

array<array<long double,N>, N> calc_beta(const int N, array<array<array<long double,3>, N>, N> x, array<array<array<long double,3>, N>, N> v, array<array<array<long double,3>, N>, N> a, array<array<long double,N>, N> alpha){
	array<array<long double,N>, N> alpha;
	for (int i; i<N; ++i){
		for (int j=0; j<i; ++j){
			beta[i][j] = (dot(v[i][j], v[i][j]) + dot(x[i][j], a[i][j]))/dot(x[i][j], x[i][j]) + alpha[i][j]*alpha[i][j];
		}
		beta[i][i] = 0.0;
		for (int j=i+1; j<N; ++j){
			beta[i][j] = (dot(v[i][j], v[i][j]) + dot(x[i][j], a[i][j]))/dot(x[i][j], x[i][j]) + alpha[i][j]*alpha[i][j];
		}
	}
    return beta;
}

array<array<long double,N>, N> calc_gamma(const int N, array<array<array<long double,3>, N>, N> x, array<array<array<long double,3>, N>, N> v, array<array<array<long double,3>, N>, N> a, array<array<array<long double,3>, N>, N> j, array<array<long double,N>, N> alpha, array<array<long double,N>, N> beta){
	array<array<long double,N>, N> gamma;
	for (int i; i<N; ++i){
		for (int k=0; k<i; ++k){
			gamma[i][k] = (3.0*dot(v[i][k], a[i][k]) + dot(x[i][k], j[i][k]))/dot(x[i][k], x[i][k]) + alpha[i][k]*(3.0*beta[i][k] - 4.0*alpha[i][k]*alpha[i][k]);
		}
		gamma[i][i] = 0.0;
		for (int k=i+1; k<N; ++k){
			gamma[i][k] = (3.0*dot(v[i][k], a[i][k]) + dot(x[i][k], j[i][k]))/dot(x[i][k], x[i][k]) + alpha[i][k]*(3.0*beta[i][k] - 4.0*alpha[i][k]*alpha[i][k]);
		}
	}
	return gamma;
}
        
tuple<array<array<long double,3>, N>, array<array<array<long double,3>, N>, N>, array<array<array<long double,3>, N>, N>> calc_acc(const int N, array<array<long double,3>, N> X, array<long double, N> M){
	array<array<long double,3>, N> acc;
	array<array<array<long double,3>, N>, N> x, A;
	//Calculate x
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
					x[i][j][k] = X[j][k]- X[i][k];
			}
			x[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
					x[i][j][k] = X[j][k]- X[i][k];
			}
		}
	}
	//Calculate A
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
					A[i][j][k] = M[j]*x[i][j][k]/(pow(norm(x[i,j]), 3.0));
			}
			A[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
					A[i][j][k] = M[j]*x[i][j][k]/(pow(norm(x[i,j]), 3.0));
			}
		}
	}
	//Calculate acc
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			acc[i][k] = 0.0;
			for (int j=0; j<N; ++j){	
					acc[i][k] += G*A[i][j][k]
			}
		}
	}		
	return make_tuple(acc, x, A);
}

tuple<array<array<long double,3>, N>, array<array<array<long double,3>, N>, N>, array<array<array<long double,3>, N>, N>, array<array<long double,N>, N>> calc_jerk(const int N, array<array<array<long double,3>, N>, N> x, array<array<long double,3>, N> V, array<long double, N> M, array<array<array<long double,3>, N>, N> A){
	array<array<long double,3>, N> jerk;
	array<array<array<long double,3>, N>, N> v, J;
	//Calculate v
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
					v[i][j][k] = V[j][k] - V[i][k];
			}
			v[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
					v[i][j][k] = V[j][k] - V[i][k];
			}
		}
	}
	//Calculate alpha
	array<array<long double,N>, N> alpha = calc_alpha(N, x, v);
	//Calculate J
	for (int k; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
				J[i][j][k] = M[j]*v[i][j][k]/pow(norm(x[i,j]), 3.0) - 3.0*alpha[i][j][k]*A[i][j][k]; 
			}
			J[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
				J[i][j][k] = M[j]*v[i][j][k]/pow(norm(x[i,j]), 3.0) - 3.0*alpha[i][j][k]*A[i][j][k];
			}
		}
	}
	//Calculate jerk
	for (int k=0; k<3 ++k){
		for (int i; i<N; ++i){
			jerk[i][k] = 0.0;
			for (int j=0; j<N; ++j){	
					jerk[i][k] += G*J[i][j][k];
			}
		}
	}		
	return make_tuple(jerk, v, J, alpha);
}
        
tuple<array<array<long double,3>, N>, array<array<array<long double,3>, N>, N>, array<array<array<long double,3>, N>, N>, array<array<long double,N>, N>> calc_snap(const int N, array<array<array<long double,3>, N>, N> x, array<array<array<long double,3>, N>, N> v, array<long double, N> M, array<array<array<long double,3>, N>, N> A, array<array<array<long double,3>, N>, N> J, array<array<long double,N>, N> alpha, array<array<long double,3>, N> acc){
	array<array<long double,3>, N> snap;
	array<array<array<long double,3>, N>, N> a, S;
	//Calculate a
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
					a[i][j][k] = acc[j][k] - acc[i][k];
			}
			a[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
					a[i][j][k] = acc[j][k] - acc[i][k];
			}
		}
	}
	//Calculate beta
	array<array<long double,N>, N> beta = calc_beta(N, x, v, a, alpha);
	//Calculate S
	for(int k; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
				S[i][j][k] = M[j]*a[i][j][k]/pow(norm(x[i,j]), 3.0) - 6.0*alpha[i][j]*J[i][j][k] - 3.0*beta[i,j]*A[i][j][k];
			}
			S[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
				S[i][j][k] = M[j]*a[i][j][k]/pow(norm(x[i,j]), 3.0) - 6.0*alpha[i][j]*J[i][j][k] - 3.0*beta[i,j]*A[i][j][k];
			}
		}
	}
	//Calculate snap
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			snap[i][k] = 0.0;
			for (int j=0; j<N; ++j){	
					snap[i][k] += G*S[i][j][k];
			}
		}		
	}
	return make_tuple(snap, a, S, beta);
}

array<array<long double,3>, N> calc_crackle(const int N, array<array<array<long double,3>, N>, N> x, array<array<array<long double,3>, N>, N> v, array<array<array<long double,3>, N>, N> a, array<long double, N> M, array<array<array<long double,3>, N>, N> A, array<array<array<long double,3>, N>, N> J, array<array<array<long double,3>, N>, N> S, array<array<long double,N>, N> alpha, array<array<long double,N>, N> beta, array<array<long double,3>, N> jerk){
	array<array<long double,3>, N> crackle;
	array<array<array<long double,3>, N>, N> j, C;
	//Calculate j
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int l=0; l<i; ++l){
					j[i][l][k] = jerk[l][k] - jerk[i][k];
			}
			j[i][i][k] = 0.0;
			for (int l=i+1; l<N; ++l){
					j[i][l][k] = jerk[l][k] - jerk[i][k];
			}
		}
	}
	//Calculate gamma
	array<array<long double,N>, N> gamma = calc_gamma(N, x, v, a, j, alpha, beta);
	//Calculate C
	for(int k; k<3; ++k){
		for (int i; i<N; ++i){
			for (int l=0; l<i; ++l){
				C[i][l][k] = M[l]*j[i][l][k]/pow(norm(x[i,l]), 3.0) - 9.0*alpha[i][l]*S[i][l][k] - 9.0*beta[i,l]*J[i][l][k] - 3.0*gamma[i,l]*A[i][l][k];
			}
			C[i][i][k] = 0.0;
			for (int l=i+1; l<N; ++j){
				C[i][l][k] = M[l]*j[i][l][k]/pow(norm(x[i,l]), 3.0) - 9.0*alpha[i][l]*S[i][l][k] - 9.0*beta[i,l]*J[i][l][k] - 3.0*gamma[i,l]*A[i][l][k];
			}
		}
	}
	//Calculate crackle
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			crackle[i][k] = 0.0;
			for (int j=0; j<N; ++j){	
					crackle[i][k] += G*C[i][j][k];
			}
		}		
	}
	return crackle;
}
  
//Aarseth criterion
long double timestep(const int N, array<array<long double,3>, N> acc, array<array<long double,3>, N> jerk, array<array<long double,3>, N> snap, array<array<long double,3>, N> crackle, long double eta){
	array<long double, N> timesteps;
	for (int i=0; i<N; ++i){
		timesteps[i] = sqrt(eta*(norm(acc[i])*norm(snap[i]) + pow(norm(jerk[i]), 2.0))/(norm(jerk[i])*norm(crackle[i]) + pow(norm(snap[i]), 2.0)));
	}
	return array_min(timesteps);
}

tuple<array<array<long double,3>, N>, array<array<long double,3>, N>, long double> acc_jerk_and_timestep(const int N, array<array<long double,3>, N> X, array<array<long double,3>, N> V, array<long double, N> M, long double eta){
	array<array<long double,3>, N> acc, jerk, snap, crackle;
	array<array<array<long double,3>, N>, N> x, A, v, J, a, S;
	array<array<long double,N>, N> alpha, beta;
	//Find acceleration
	tuple<array<array<long double,3>, N>, array<array<array<long double,3>, N>, N>, array<array<array<long double,3>, N>, N>> accxA = calc_acc(N, X, M);
	acc = get<0>(accxA);
	x = get<1>(accxA);
	A = get<2>(accxA);
	//Find jerk
	tuple<array<array<long double,3>, N>, array<array<array<long double,3>, N>, N>, array<array<array<long double,3>, N>, N>, array<array<long double,N>, N>> jerkvJalpha = calc_jerk(N, x, V, M, A);
	jerk = get<0>(jerkvJalpha);
	v = get<1>(jerkvJalpha);
	J = get<2>(jerkvJalpha);
	alpha = get<3>(jerkvJalpha);
	//Find snap
	tuple<array<array<long double,3>, N>, array<array<array<long double,3>, N>, N>, array<array<array<long double,3>, N>, N>, array<array<long double,N>, N>> snapaSbeta = calc_snap(N, x, v, M, A, J, alpha, acc);
	snap = get<0>(snapaSbeta);
	a = get<1>(snapaSbeta);
	S = get<2>(snapaSbeta);
	beta = get<3>(snapaSbeta);
	//Find crackle
	crackle = calc_crackle(N, x, v, a, M, A, J, S, alpha, beta, jerk);
	//Find timestep
	long double dt = timestep(N, acc, jerk, snap, crackle, eta);
	return make_tuple(acc,jerk,dt);
}

tuple<array<array<long double,3>, N>, array<array<long double,3>, N>> acc_and_jerk(const int N, array<array<long double,3>, N> X, array<array<long double,3>, N> V, array<long double, N> M)
{
	array<array<long double,3>, N> acc, jerk;
	array<array<array<long double,3>, N>, N> x, A;
	//Find acceleration
	tuple<array<array<long double,3>, N>, array<array<array<long double,3>, N>, N>, array<array<array<long double,3>, N>, N>> accxA = calc_acc(N, X, M);
	acc = get<0>(accxA);
	x = get<1>(accxA);
	A = get<2>(accxA);
	//Find jerk
	tuple<array<array<long double,3>, N>, array<array<array<long double,3>, N>, N>, array<array<array<long double,3>, N>, N>, array<array<long double,N>, N>> jerkvJalpha = calc_jerk(N, x, V, M, A);
	jerk = get<0>(jerkvJalpha);
	return make_tuple(acc, jerk);
}          