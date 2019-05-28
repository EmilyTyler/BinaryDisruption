//
#include <cmath>
#include <array>
#include <tuple>
#include <vector>
#include <iostream>
#include "vector_maths.h"
#include "constants.h"
using namespace std;



//See Nitadori and Makino 2008 for higher derivatives of acceleration
vector<vector<long double>> calc_alpha(int N, vector<vector<array<long double,3>>> x, vector<vector<array<long double,3>>> v){
	vector<vector<long double>> alpha;
	alpha.resize(N);
	for (int i=0; i<N; ++i){
		alpha[i].resize(N);
	}
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

vector<vector<long double>> calc_beta(int N, vector<vector<array<long double,3>>> x, vector<vector<array<long double,3>>> v, vector<vector<array<long double,3>>> a, vector<vector<long double>> alpha){
	vector<vector<long double>> beta;
	beta.resize(N);
	for (int i=0; i<N; ++i){
		beta[i].resize(N);
	}
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

vector<vector<long double>> calc_gamma(int N, vector<vector<array<long double,3>>> x, vector<vector<array<long double,3>>> v, vector<vector<array<long double,3>>> a, vector<vector<array<long double,3>>> j, vector<vector<long double>> alpha, vector<vector<long double>> beta){
	vector<vector<long double>> gamma;
	gamma.resize(N);
	for (int i=0; i<N; ++i){
		gamma[i].resize(N);
	}
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
        
tuple<vector<array<long double,3>>, vector<vector<array<long double,3>>>, vector<vector<array<long double,3>>>> calc_acc(int N, vector<array<long double,3>> X, vector<long double> M){
	vector<array<long double,3>> acc;
	acc.resize(N);
	vector<vector<array<long double,3>>> x, A;
	x.resize(N);
	for (int i=0; i<N; ++i){
		x[i].resize(N);
	}
	A.resize(N);
	for (int i=0; i<N; ++i){
		A[i].resize(N);
	}
	//Calculate x
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<N; ++j){
					x[i][j][k] = X[j][k]- X[i][k];
			}
		}
	}
	//Calculate A
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
					A[i][j][k] = M[j]*x[i][j][k]/(pow(norm(x[i][j]), 3.0));
			}
			A[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
					A[i][j][k] = M[j]*x[i][j][k]/(pow(norm(x[i][j]), 3.0));
			}
		}
	}
	//Calculate acc
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			acc[i][k] = 0.0;
			for (int j=0; j<N; ++j){	
					acc[i][k] += G*A[i][j][k];
			}
		}
	}		
	return make_tuple(acc, x, A);
}

tuple<vector<array<long double,3>>, vector<vector<array<long double,3>>>, vector<vector<array<long double,3>>>, vector<vector<long double>>> calc_jerk(int N, vector<vector<array<long double,3>>> x, vector<array<long double,3>> V, vector<long double> M, vector<vector<array<long double,3>>> A){
	vector<array<long double,3>> jerk;
	jerk.resize(N);
	vector<vector<array<long double,3>>> v, J;
	v.resize(N);
	for (int i=0; i<N; ++i){
		v[i].resize(N);
	}
	J.resize(N);
	for (int i=0; i<N; ++i){
		J[i].resize(N);
	}
	//Calculate v
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<N; ++j){
					v[i][j][k] = V[j][k] - V[i][k];
			}
		}
	}
	//Calculate alpha
	vector<vector<long double>> alpha = calc_alpha(N, x, v);
	//Calculate J
	for (int k; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
				J[i][j][k] = M[j]*v[i][j][k]/pow(norm(x[i][j]), 3.0) - 3.0*alpha[i][j]*A[i][j][k]; 
			}
			J[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
				J[i][j][k] = M[j]*v[i][j][k]/pow(norm(x[i][j]), 3.0) - 3.0*alpha[i][j]*A[i][j][k];
			}
		}
	}
	//Calculate jerk
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			jerk[i][k] = 0.0;
			for (int j=0; j<N; ++j){	
					jerk[i][k] += G*J[i][j][k];
			}
		}
	}		
	return make_tuple(jerk, v, J, alpha);
}
        
tuple<vector<array<long double,3>>, vector<vector<array<long double,3>>>, vector<vector<array<long double,3>>>, vector<vector<long double>>> calc_snap(int N, vector<vector<array<long double,3>>> x, vector<vector<array<long double,3>>> v, vector<long double> M, vector<vector<array<long double,3>>> A, vector<vector<array<long double,3>>> J, vector<vector<long double>> alpha, vector<array<long double,3>> acc){
	vector<array<long double,3>> snap;
	snap.resize(N);
	vector<vector<array<long double,3>>> a, S;
	a.resize(N);
	for (int i=0; i<N; ++i){
		a[i].resize(N);
	}
	S.resize(N);
	for (int i=0; i<N; ++i){
		S[i].resize(N);
	}
	cout << "Here11" << endl;
	//Calculate a
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<N; ++j){
					a[i][j][k] = acc[j][k] - acc[i][k];
			}
		}
	}
	cout << "Here12" << endl;
	//Calculate beta
	vector<vector<long double>> beta = calc_beta(N, x, v, a, alpha);
	//Calculate S
	for(int k; k<3; ++k){
		for (int i; i<N; ++i){
			for (int j=0; j<i; ++j){
				S[i][j][k] = M[j]*a[i][j][k]/pow(norm(x[i][j]), 3.0) - 6.0*alpha[i][j]*J[i][j][k] - 3.0*beta[i][j]*A[i][j][k];
			}
			S[i][i][k] = 0.0;
			for (int j=i+1; j<N; ++j){
				S[i][j][k] = M[j]*a[i][j][k]/pow(norm(x[i][j]), 3.0) - 6.0*alpha[i][j]*J[i][j][k] - 3.0*beta[i][j]*A[i][j][k];
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
	cout << "Here22" << endl;
	return make_tuple(snap, a, S, beta);
}

vector<array<long double,3>> calc_crackle(int N, vector<vector<array<long double,3>>> x, vector<vector<array<long double,3>>> v, vector<vector<array<long double,3>>> a, vector<long double> M, vector<vector<array<long double,3>>> A, vector<vector<array<long double,3>>> J, vector<vector<array<long double,3>>> S, vector<vector<long double>> alpha, vector<vector<long double>> beta, vector<array<long double,3>> jerk){
	vector<array<long double,3>> crackle;
	crackle.resize(N);
	vector<vector<array<long double,3>>> j, C;
	j.resize(N);
	for (int i=0; i<N; ++i){
		j[i].resize(N);
	}
	C.resize(N);
	for (int i=0; i<N; ++i){
		C[i].resize(N);
	}
	//Calculate j
	for (int k=0; k<3; ++k){
		for (int i; i<N; ++i){
			for (int l=0; l<N; ++l){
					j[i][l][k] = jerk[l][k] - jerk[i][k];
			}
		}
	}
	//Calculate gamma
	vector<vector<long double>> gamma = calc_gamma(N, x, v, a, j, alpha, beta);
	//Calculate C
	for(int k; k<3; ++k){
		for (int i; i<N; ++i){
			for (int l=0; l<i; ++l){
				C[i][l][k] = M[l]*j[i][l][k]/pow(norm(x[i][l]), 3.0) - 9.0*alpha[i][l]*S[i][l][k] - 9.0*beta[i][l]*J[i][l][k] - 3.0*gamma[i][l]*A[i][l][k];
			}
			C[i][i][k] = 0.0;
			for (int l=i+1; l<N; ++l){
				C[i][l][k] = M[l]*j[i][l][k]/pow(norm(x[i][l]), 3.0) - 9.0*alpha[i][l]*S[i][l][k] - 9.0*beta[i][l]*J[i][l][k] - 3.0*gamma[i][l]*A[i][l][k];
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
long double timestep(int N, vector<array<long double,3>> acc, vector<array<long double,3>> jerk, vector<array<long double,3>> snap, vector<array<long double,3>> crackle, long double eta){
	vector<long double> timesteps;
	timesteps.resize(N);
	for (int i=0; i<N; ++i){
		timesteps[i] = sqrt(eta*(norm(acc[i])*norm(snap[i]) + pow(norm(jerk[i]), 2.0))/(norm(jerk[i])*norm(crackle[i]) + pow(norm(snap[i]), 2.0)));
	}
	long double min_dt = timesteps[0];
	for (int i=1; i<N; ++i){
		min_dt = min(timesteps[i], timesteps[i-1]);
	}
	return min_dt;
}

tuple<vector<array<long double,3>>, vector<array<long double,3>>, long double> acc_jerk_and_timestep(int N, vector<array<long double,3>> X, vector<array<long double,3>> V, vector<long double> M, long double eta){
	vector<array<long double,3>> acc, jerk, snap, crackle;
	vector<vector<array<long double,3>>> x, A, v, J, a, S;
	vector<vector<long double>> alpha, beta;
	//Find acceleration
	tuple<vector<array<long double,3>>, vector<vector<array<long double,3>>>, vector<vector<array<long double,3>>>> accxA = calc_acc(N, X, M);
	acc = get<0>(accxA);
	x = get<1>(accxA);
	A = get<2>(accxA);
	//Find jerk
	tuple<vector<array<long double,3>>, vector<vector<array<long double,3>>>, vector<vector<array<long double,3>>>, vector<vector<long double>>> jerkvJalpha = calc_jerk(N, x, V, M, A);
	jerk = get<0>(jerkvJalpha);
	v = get<1>(jerkvJalpha);
	J = get<2>(jerkvJalpha);
	alpha = get<3>(jerkvJalpha);
	//Find snap
	cout << "Here1" << endl;
	tuple<vector<array<long double,3>>, vector<vector<array<long double,3>>>, vector<vector<array<long double,3>>>, vector<vector<long double>>> snapaSbeta = calc_snap(N, x, v, M, A, J, alpha, acc);
	cout << "Here2" << endl;
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

tuple<vector<array<long double,3>>, vector<array<long double,3>>> acc_and_jerk(int N, vector<array<long double,3>> X, vector<array<long double,3>> V, vector<long double> M)
{
	vector<array<long double,3>> acc, jerk;
	vector<vector<array<long double,3>>> x, A;
	//Find acceleration
	tuple<vector<array<long double,3>>, vector<vector<array<long double,3>>>, vector<vector<array<long double,3>>>> accxA = calc_acc(N, X, M);
	acc = get<0>(accxA);
	x = get<1>(accxA);
	A = get<2>(accxA);
	//Find jerk
	tuple<vector<array<long double,3>>, vector<vector<array<long double,3>>>, vector<vector<array<long double,3>>>, vector<vector<long double>>> jerkvJalpha = calc_jerk(N, x, V, M, A);
	jerk = get<0>(jerkvJalpha);
	return make_tuple(acc, jerk);
}

//Evolution of binary through integration
long double singleTimestep(int N, vector<array<long double, 3>> &X, vector<long double> M, int n, long double eta, long double dt_max=-1.0)
{
    //dt is step size
    //X = [x1, x2, x3,..., v1, v2, v3,...]
    //N is number of objects
    //M = [m1, m2, m3,...]
    //n is number of iterations of estimate and correct
	vector<array<long double,3>> X_0, V_0, A_0, J_0, X_1, V_1, A_1, J_1;
	X_0.resize(N);
	V_0.resize(N);
	A_0.resize(N);
	J_0.resize(N);
	X_1.resize(N);
	V_1.resize(N);
	A_1.resize(N);
	J_1.resize(N);
	X_0.shrink_to_fit();
	V_0.shrink_to_fit();
	A_0.shrink_to_fit();
	J_0.shrink_to_fit();
	X_1.shrink_to_fit();
	V_1.shrink_to_fit();
	A_1.shrink_to_fit();
	J_1.shrink_to_fit();

	long double dt;
	for (int i=0; i<N; ++i){
		for (int j=0; j<3; ++j){
			X_0[i][j] = X[i][j];
			V_0[i][j] = X[i+N][j];
		}
	} 
    //Hermite scheme (Dehnen and Read 2011) with Aarseth criterion timesteps  
    //Find initial acceleration and jerk and timestep
    tuple<vector<array<long double,3>>, vector<array<long double,3>>, long double> ajt = acc_jerk_and_timestep(N, X_0, V_0, M, eta); 
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
        tuple<vector<array<long double,3>>, vector<array<long double,3>>> aj = acc_and_jerk(N, X_1, V_1, M);
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
    return dt;
}

vector<array<long double, 3>> evolve(int N, vector<long double> M, vector<array<long double, 3>> X, long double T, int n=1, long double eta = 0.02)
{
	//Current time
	long double t = 0.0;
	long double dt_max, dt;
	while (t<T){
		dt_max = T - t;
		dt = singleTimestep(N, X, M, n, eta, dt_max=dt_max);
		t += dt;
	}
	return X;
}          