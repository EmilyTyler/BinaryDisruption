//
#include <cmath>
#include <array>
#include <tuple>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vector_maths.h"
#include "constants.h"
using namespace std;

vector<vector<long double>> alpha;
vector<vector<long double>> beta;
vector<vector<long double>> gama;

//x[i][j][k] = X[j][k] - X[i][k]
vector<vector<array<long double,3>>> x;
vector<vector<array<long double,3>>> v;
vector<vector<array<long double,3>>> a;
vector<vector<array<long double,3>>> j;

//Acceleration etc
vector<array<long double,3>> acc;
vector<array<long double,3>> jerk;
vector<array<long double,3>> snap;
vector<array<long double,3>> crackle;

//
vector<vector<array<long double,3>>> A;
vector<vector<array<long double,3>>> J;
vector<vector<array<long double,3>>> S;
vector<vector<array<long double,3>>> C;


//See Nitadori and Makino 2008 for higher derivatives of acceleration
//Depends on x and v
void calc_alpha(int N){
	for (int i=0; i<N; ++i){
		for (int l=0; l<i; ++l){
			alpha[i][l] = dot(x[i][l], v[i][l])/dot(x[i][l], x[i][l]);
		}
		alpha[i][i] = 0.0;
		for (int l=i+1; l<N; ++l){
			alpha[i][l] = dot(x[i][l], v[i][l])/dot(x[i][l], x[i][l]);
		}
	}
}

//Depends on x, v, a and alpha
void calc_beta(int N){
	for (int i=0; i<N; ++i){
		for (int l=0; l<i; ++l){
			beta[i][l] = (dot(v[i][l], v[i][l]) + dot(x[i][l], a[i][l]))/dot(x[i][l], x[i][l]) + alpha[i][l]*alpha[i][l];
		}
		beta[i][i] = 0.0;
		for (int l=i+1; l<N; ++l){
			beta[i][l] = (dot(v[i][l], v[i][l]) + dot(x[i][l], a[i][l]))/dot(x[i][l], x[i][l]) + alpha[i][l]*alpha[i][l];
		}
	}
}

//Depends on x, v, a, j, alpha, beta
void calc_gama(int N){
	for (int i=0; i<N; ++i){
		for (int k=0; k<i; ++k){
			gama[i][k] = (3.0*dot(v[i][k], a[i][k]) + dot(x[i][k], j[i][k]))/dot(x[i][k], x[i][k]) + alpha[i][k]*(3.0*beta[i][k] - 4.0*alpha[i][k]*alpha[i][k]);
		}
		gama[i][i] = 0.0;
		for (int k=i+1; k<N; ++k){
			gama[i][k] = (3.0*dot(v[i][k], a[i][k]) + dot(x[i][k], j[i][k]))/dot(x[i][k], x[i][k]) + alpha[i][k]*(3.0*beta[i][k] - 4.0*alpha[i][k]*alpha[i][k]);
		}
	}
}
        
void calc_acc(int N, vector<array<long double,3>> X, vector<long double> M){
	//Calculate x
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			for (int l=0; l<N; ++l){
					x[i][l][k] = X[l][k]- X[i][k];
			}
		}
	}	
	//Calculate A
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			for (int l=0; l<i; ++l){
					A[i][l][k] = M[l]*x[i][l][k]/(pow(norm(x[i][l]), 3.0));
			}
			A[i][i][k] = 0.0;
			for (int l=i+1; l<N; ++l){
					A[i][l][k] = M[l]*x[i][l][k]/(pow(norm(x[i][l]), 3.0));
			}
		}
	}
	//Calculate acc
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			acc[i][k] = 0.0;
			for (int l=0; l<N; ++l){	
					acc[i][k] += G*A[i][l][k];
			}
		}
	}	
	return;
}

//Depends on x, A
void calc_jerk(int N, vector<array<long double,3>> V, vector<long double> M){
	//Calculate v
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			for (int l=0; l<N; ++l){
					v[i][l][k] = V[l][k] - V[i][k];
			}
		}
	}
	//Calculate alpha
	calc_alpha(N);
	//Calculate J
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			for (int l=0; l<i; ++l){
				J[i][l][k] = M[l]*v[i][l][k]/pow(norm(x[i][l]), 3.0) - 3.0*alpha[i][l]*A[i][l][k]; 
			}
			J[i][i][k] = 0.0;
			for (int l=i+1; l<N; ++l){
				J[i][l][k] = M[l]*v[i][l][k]/pow(norm(x[i][l]), 3.0) - 3.0*alpha[i][l]*A[i][l][k];
			}
		}
	}
	//Calculate jerk
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			jerk[i][k] = 0.0;
			for (int l=0; l<N; ++l){	
					jerk[i][k] += G*J[i][l][k];
			}
		}
	}		
}

//Depends on x, v, A, J, acc, alpha        
void calc_snap(int N, vector<long double> M){
	//Calculate a
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			for (int l=0; l<N; ++l){
					a[i][l][k] = acc[l][k] - acc[i][k];
			}
		}
	}
	//Calculate beta
	calc_beta(N);
	//Calculate S
	for(int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			for (int l=0; l<i; ++l){
				S[i][l][k] = M[l]*a[i][l][k]/pow(norm(x[i][l]), 3.0) - 6.0*alpha[i][l]*J[i][l][k] - 3.0*beta[i][l]*A[i][l][k];
			}
			S[i][i][k] = 0.0;
			for (int l=i+1; l<N; ++l){
				S[i][l][k] = M[l]*a[i][l][k]/pow(norm(x[i][l]), 3.0) - 6.0*alpha[i][l]*J[i][l][k] - 3.0*beta[i][l]*A[i][l][k];
			}
		}
	}
	//Calculate snap
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			snap[i][k] = 0.0;
			for (int l=0; l<N; ++l){	
					snap[i][k] += G*S[i][l][k];
			}
		}		
	}
}

//Depends on x, v, a, A, J, S, alpha, beta, jerk
void calc_crackle(int N, vector<long double> M){
	//Calculate j
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			for (int l=0; l<N; ++l){
					j[i][l][k] = jerk[l][k] - jerk[i][k];
			}
		}
	}
	//Calculate gama
	calc_gama(N);
	//Calculate C
	for(int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			for (int l=0; l<i; ++l){
				C[i][l][k] = M[l]*j[i][l][k]/pow(norm(x[i][l]), 3.0) - 9.0*alpha[i][l]*S[i][l][k] - 9.0*beta[i][l]*J[i][l][k] - 3.0*gama[i][l]*A[i][l][k];
			}
			C[i][i][k] = 0.0;
			for (int l=i+1; l<N; ++l){
				C[i][l][k] = M[l]*j[i][l][k]/pow(norm(x[i][l]), 3.0) - 9.0*alpha[i][l]*S[i][l][k] - 9.0*beta[i][l]*J[i][l][k] - 3.0*gama[i][l]*A[i][l][k];
			}
		}
	}
	//Calculate crackle
	for (int k=0; k<3; ++k){
		for (int i=0; i<N; ++i){
			crackle[i][k] = 0.0;
			for (int l=0; l<N; ++l){	
					crackle[i][k] += G*C[i][l][k];
			}
		}		
	}
}

//Aarseth criterion
//Depends on acc, jerk, snap, crackle
long double timestep(int N, long double eta){
	vector<long double> timesteps;
	timesteps.resize(N);
	for (int i=0; i<N; ++i){
		//cout << sqrt(eta*(norm(acc[i])*norm(snap[i]) + pow(norm(jerk[i]), 2.0))/(norm(jerk[i])*norm(crackle[i]) + pow(norm(snap[i]), 2.0))) << endl;
		timesteps[i] = sqrt(eta*(norm(acc[i])*norm(snap[i]) + pow(norm(jerk[i]), 2.0))/(norm(jerk[i])*norm(crackle[i]) + pow(norm(snap[i]), 2.0)));
	}
	long double min_dt = timesteps[0];
	for (int i=1; i<N; ++i){
		min_dt = min(timesteps[i], timesteps[i-1]);
	}
	return min_dt;
}

long double acc_jerk_and_timestep(int N, vector<array<long double,3>> X, vector<array<long double,3>> V, vector<long double> M, long double eta){
	//Calculate acc, x, A
	calc_acc(N, X, M);
	//Calculate jerk, v, J, alpha
	calc_jerk(N, V, M);
	//Calculate snap, a, S, beta
	calc_snap(N, M);
	//Calculate crackle
	calc_crackle(N, M);
	//Find timestep
	long double dt = timestep(N, eta);
	//long double dt=0.0000002*giga*year/time_scale;
	return dt;
}

void acc_and_jerk(int N, vector<array<long double,3>> X, vector<array<long double,3>> V, vector<long double> M)
{
	//Calculate acc, x, A
	calc_acc(N, X, M);
	//Calculate jerk, v, J, alpha
	calc_jerk(N, V, M);
	return;
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
		for (int l=0; l<3; ++l){
			X_0[i][l] = X[i][l];
			V_0[i][l] = X[i+N][l];
		}
	} 
    //Hermite scheme (Dehnen and Read 2011) with Aarseth criterion timesteps  
    //Find initial acceleration and jerk and timestep
    //cout << "Total energy1, internal units = " << M[0]*M[1]*(dot(v[1][0], v[1][0])/(2*(M[0]+M[1])) - G/(norm(x[1][0]))) << endl;
    dt = acc_jerk_and_timestep(N, X_0, V_0, M, eta); 
    //acc_and_jerk(N, X_0, V_0, M);
    //dt = 7*pow(10.0, -6.0);
    //cout << "Total energy2, internal units = " << M[0]*M[1]*(dot(v[1][0], v[1][0])/(2*(M[0]+M[1])) - G/(norm(x[1][0]))) << endl;
    A_0 = acc;
    J_0 = jerk;
    if ((dt_max>0.0) && (dt_max<dt)){
    	dt = dt_max;
    }
    //Predict positions and velocities
    for (int i=0; i<N; ++i){
    	for (int l=0; l<3; ++l){
    		X_1[i][l] = X_0[i][l] + V_0[i][l]*dt + A_0[i][l]*dt*dt/2.0 + J_0[i][l]*pow(dt,3.0)/6.0;
    		V_1[i][l] = V_0[i][l] + A_0[i][l]*dt + J_0[i][l]*dt*dt/2.0;
    	}
    }
    //Iterate n times
    for (int k=0; k<n; ++k){   
        //Estimate acceleration and jerk
        acc_and_jerk(N, X_1, V_1, M);
        //cout << "Total energy3, internal units = " << M[0]*M[1]*(dot(v[1][0], v[1][0])/(2*(M[0]+M[1])) - G/(norm(x[1][0]))) << endl;
        A_1 = acc;
        J_1 = jerk;     
        //Obtain corrected positions and velocities
        for (int i=0; i<N; ++i){
        	for (int l=0; l<3; ++l){
        		X_1[i][l] = X_0[i][l] + (V_1[i][l]+V_0[i][l])*dt/2.0 + (A_0[i][l]-A_1[i][l])*dt*dt/12.0;
	        	V_1[i][l] = V_0[i][l] + (A_1[i][l]+A_0[i][l])*dt/2.0 + (J_0[i][l]-J_1[i][l])*dt*dt/12.0;
        	}
        }
    }
    for (int i=0; i<N; ++i){
		for (int l=0; l<3; ++l){
			X[i][l] = X_1[i][l];
			X[i+N][l] = V_1[i][l];
		}
	} 
	//cout << "Total energy4, internal units = " << M[0]*M[1]*(dot(v[1][0], v[1][0])/(2*(M[0]+M[1])) - G/(norm(x[1][0]))) << endl;
	//cin.ignore();
    return dt;
}

void initialise_arrays(int N)
{
	alpha.resize(N);
	alpha.shrink_to_fit();
	for (int i=0; i<N; ++i){
		alpha[i].resize(N);
		alpha[i].shrink_to_fit();
	}
	beta.resize(N);
	beta.shrink_to_fit();
	for (int i=0; i<N; ++i){
		beta[i].resize(N);
		beta[i].shrink_to_fit();
	}
	gama.resize(N);
	gama.shrink_to_fit();
	for (int i=0; i<N; ++i){
		gama[i].resize(N);
		gama[i].shrink_to_fit();
	}
	acc.resize(N);
	acc.shrink_to_fit();
	x.resize(N);
	x.shrink_to_fit();
	for (int i=0; i<N; ++i){
		x[i].resize(N);
		x[i].shrink_to_fit();
	}
	A.resize(N);
	A.shrink_to_fit();
	for (int i=0; i<N; ++i){
		A[i].resize(N);
		A[i].shrink_to_fit();
	}
	jerk.resize(N);
	jerk.shrink_to_fit();
	v.resize(N);
	v.shrink_to_fit();
	for (int i=0; i<N; ++i){
		v[i].resize(N);
		v[i].shrink_to_fit();
	}
	J.resize(N);
	J.shrink_to_fit();
	for (int i=0; i<N; ++i){
		J[i].resize(N);
		J[i].shrink_to_fit();
	}
	snap.resize(N);
	snap.shrink_to_fit();
	a.resize(N);
	a.shrink_to_fit();
	for (int i=0; i<N; ++i){
		a[i].resize(N);
		a[i].shrink_to_fit();
	}
	S.resize(N);
	S.shrink_to_fit();
	for (int i=0; i<N; ++i){
		S[i].resize(N);
		S[i].shrink_to_fit();
	}
	crackle.resize(N);
	crackle.shrink_to_fit();
	j.resize(N);
	j.shrink_to_fit();
	for (int i=0; i<N; ++i){
		j[i].resize(N);
		j[i].shrink_to_fit();
	}
	C.resize(N);
	C.shrink_to_fit();
	for (int i=0; i<N; ++i){
		C[i].resize(N);
		C[i].shrink_to_fit();
	}
	return;
}

vector<array<long double, 3>> evolve(int N, vector<long double> M, vector<array<long double, 3>> X, long double T, int n=1, long double eta = 0.02, bool ini_arrays = true)
{
	//Current time
	long double t = 0.0;
	long double dt_max, dt;
	if (ini_arrays){
		//cout << "Initialising arrays" << endl;
		initialise_arrays(N);
	}
	n=10;
	eta = 0.00002;
	//ofstream myfile;
	//myfile.open("test_nbody.csv");
	while (t<T){
		//cout << setprecision(16) << "Current time, Gyr = " << t*time_scale/(giga*year) << endl;
		dt_max = T - t;
		dt = singleTimestep(N, X, M, n, eta, dt_max=dt_max);
		//cout << "Timestep, Gyr = " << dt*time_scale/(giga*year) << endl;
		t += dt;
		//myfile << setprecision(16) << X[0][0] << ", " << X[0][1] << ", " << X[0][2] << ", " << X[1][0] << ", " << X[1][1] << ", " << X[1][2] << endl;
	}
	//myfile.close();
	return X;
}          