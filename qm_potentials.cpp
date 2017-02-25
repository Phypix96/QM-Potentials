#include "qm_pot_gui.h"

double m_Pi = 3.141592653;

/*
QVector<double> bisection(double k){

    int l = sqrt(2.*k)/m_Pi+1;
    QVector<double> modes(2*l);
    double a1, a2, b1, b2, x;

    for(int i = 0; i < l; ++i ){
        a1 = m_Pi*i;
        b1 = m_Pi*(i+1);
        a2 = a1;
        b2 = b1;

        for(int j=0; j<40; j++){
            x = (a1+b1)/2.;
            if(tan(x)-sqrt(2.*k/(x*x)-1.) < 0) a1=x;
            else b1=x;
        }
        for(int j=0; j<40; j++){
            x = (a2+b2)/2.;
            if(-1./tan(x)-sqrt(2.*k/(x*x)-1.) < 0) a2=x;
            else b2=x;
        }
        modes[2*i] = (a1+b1)/2.;
        modes[2*i+1] = (a2+b2)/2.;
    }

    return modes;
}
*/

