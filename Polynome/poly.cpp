#include <iostream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

double factorial(int n){
    double fac = 1.;
    for(int i = 1; i <= n; ++i){
        fac *= i;
    }
    return fac;
}

void Hermite(int n,double h[]){
    //double h[n+1] = {};
    h[0] = 1.;
    for(int i = 1; i <= n; ++i){
        for(int j = i%2; j <= i; j+=2){
            if(j == 0){
                h[j] = (j+1.)*h[j+1];
            }
            if(j == i){
                h[j] = -2.*h[j-1];
            }
            if(j != 0 && j != i){
                h[j] = (j+1.)*h[j+1] - 2.*h[j-1];
            }
        }
    }

}

double Hermite_poly(double x, double w, double m, int n, double h[]){
    double x_prime = sqrt(m*w)*x;
    double phi = h[n] * x_prime;

    for(int i = n-2; i > 0; i -= 2){
        phi += h[i];
        phi *= x_prime*x_prime;
    }
    if((n+1)%2) phi += h[0];

    phi *= pow(m*w/m_Pi,0.25)*exp(-x_prime*x_prime/2.)/(pow(2.,n/2.)*sqrt(fac(n)));
    return phi;

}
/*
Andere Funktion so schnell, Zeiterspanis irrelavant
void Hermite(unsigned int n, unsigned int k, double h[]){ //h muss bereits richtige Länge haben, vor Funktionsaufruf klären
	if(n == 0) h[0] = 1.;

	if(n < k){
		for(int i = n+1; i <= k; ++i){
        	for(int j = i%2; j <= i; j+=2){
 	           if(j == 0){
    	            h[j] = (j+1.)*h[j+1];
	            }
	            if(j == i){
    	            h[j] = -2.*h[j-1];
        	    }
            	if(j != 0 && j != i){
   	             h[j] = (j+1.)*h[j+1] - 2.*h[j-1];
    	        }
    	    }
   		}
	}
	if(n > k){

	}
}
*/
void Laguerre(int n, int k){
    double l[n][2] = {};
    l[0][0] = 1.;
    for(int i = 1; i <= n; ++i){
        for(int j = 0; j <= i+1; ++j){
            if(j == 0){
                l[j][i%2] = -l[j][(i+1)%2];
            }
            if(j == i){
                l[j][i%2] = (n+k-(j-1))*l[j-1][(i+1)%2];
            }
            if(j != 0 && j != i){
                l[j][i%2] = (n+k-(j-1))*l[j-1][(i+1)%2] - l[j][(i+1)%2];
            }
        }
    }

}

void Legendre(int l, int m){
    double c [l+m+1][l+1] = {};
    c[0][0]=1.;
    for(int n = 1; n <= l+m; ++n){
        for(int k = 0; k <= min(n,l); ++k){ //es reicht vermutlich bis l zu gehen, da dann die Einträge außen 0 werden
            if(k <= n && n%2 == k%2){
                if(k == 0){
                   c[n][k] = c[n-1][k+1]; //*(0+1)
                }
                if(k == n){
                    c[n][k] = pow(2.,n)*factorial(l)/factorial(l-n);
                }
                if(k != 0 && k!= n){
                    c[n][k] = 2.*(l-(k+n)/2.+1)*c[n-1][k-1] + (k+1)*c[n-1][k+1];
                }
            }
        }
    }

    for(int i = 0; i <= l+m; ++i){
        for(int j = 0; j <= min(i, l); ++j ){
            cout<<c[i][j]<<" ";
        }
        cout<<endl;
    }

    //c[l+m+1] = *koef;
}



int main(){
    int l = 250, n = 100;
    
    double h[(l/50+1)*50];

	Hermite(3,h);
	cout<<h[3];
/*
    for( int i = 0; i<=l/2; ++i){
    	if(l%2) cout<<h[l-2*i]<<" ";
    	else cout<<h[l-2*i]<<" ";	
    }
    //save vector with coefficients as private variables of ui and change them, when the state is changed!
*/
}
