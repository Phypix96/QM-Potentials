#include "qm_pot_gui.h"
#include "ui_qm_pot_gui.h"

double const m_Pi = 3.141596253;

double fac(double l){
    if(l >= 0){
        double fac_n = 1.;
        for(int i = 1; i <= l; ++i){
            fac_n *= i;
        }
        return fac_n;
    }
    //error handling missing!!!!!
    else{
        return 0.;
    }
}

QVector<double> bisection(double a, double k0){

    int l = sqrt(k0/2.)*a/m_Pi;
    QVector<double> modes(l*2+2);
    double a1, a2, b1, b2, x;

    for(int i = 0; i <= l; ++i ){
        a1 = m_Pi*(i-1./2.);
        b1 = m_Pi*(i+1./2.);
        a2 = m_Pi*i;
        b2 = m_Pi*(i+1.);
        for(int j=0; j<40; j++){
            x = (a1+b1)/2.;
            if(tan(x)-sqrt(k0*a*a/(2.*x*x)-1.) < 0) a1=x;
            else b1=x;
        }
        for(int j=0; j<40; j++){
            x = (a2+b2)/2.;
            if(1./tan(x)+sqrt(k0*a*a/(2.*x*x)-1.) > 0) a2=x;
            else b2=x;
        }
        modes[2*i] = (a1+b1)/a;
        modes[2*i+1] = (a2+b2)/a;
    }
    return modes;
}



void Hermite(int n,double h[]){
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
    double x_prime_square = m*w*x*x;
    double phi = h[n] * x_prime_square;

    for(int i = n-2; i > 0; i -= 2){
        phi += h[i];
        phi *= x_prime_square;
    }
    if((n+1)%2) phi += h[0];
    else phi /= (sqrt(m*w)*x);

    phi *= pow(m*w/m_Pi,0.25)*exp(-x_prime_square/2.)/(pow(2.,n/2.)*sqrt(fac(n)));
    return phi;

}

//////////////////////////////////////
/// \brief QM_Pot_GUI::QM_Pot_GUI
/// \param parent
///


QM_Pot_GUI::QM_Pot_GUI(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::QM_Pot_GUI)
{
    ui->setupUi(this);

    connect(ui->actionQuit, SIGNAL(triggered(bool)), this, SLOT(close()));
    connect(ui->actionPDF, SIGNAL(triggered(bool)), this, SLOT(print_PDF()));
    connect(ui->actionPNG, SIGNAL(triggered(bool)), this, SLOT(print_PNG()));

    infinite_init();
    finite_init();
    harmonic_init();

}


QM_Pot_GUI::~QM_Pot_GUI()
{
    delete ui;
}

void QM_Pot_GUI::infinite_init()
{
    connect(ui->infinite_mode_slide, SIGNAL(valueChanged(int)), ui->infinite_mode_spin, SLOT(setValue(int)));
    connect(ui->infinite_mode_spin, SIGNAL(valueChanged(int)), ui->infinite_mode_slide, SLOT(setValue(int)));
    connect(ui->infinite_width_slide, SIGNAL(valueChanged(int)), ui->infinite_width_spin, SLOT(setValue(int)));
    connect(ui->infinite_width_spin, SIGNAL(valueChanged(int)), ui->infinite_width_slide, SLOT(setValue(int)));

    QPen pen;
    pen.setColor(QColor(0, 20, 200, 100));
    pen.setWidthF(0);
    ui->infinite_graph->yAxis->setRange(-1.1,1.1);
    ui->infinite_graph->xAxis->setRange(-25,25);
    ui->infinite_graph->yAxis2->setRange(0,1);
    ui->infinite_graph->addGraph();
    ui->infinite_graph->addGraph(ui->infinite_graph->xAxis, ui->infinite_graph->yAxis2);
    ui->infinite_graph->graph(1)->setPen(pen);
    ui->infinite_graph->graph(1)->setBrush(QBrush(QColor(0,20,200,100)));

    connect(ui->infinite_width_spin, SIGNAL(valueChanged(int)), this, SLOT(infinite_graph_update()));
    connect(ui->infinite_mode_spin, SIGNAL(valueChanged(int)), this, SLOT(infinite_graph_update()));
    connect(ui->actionProbability, SIGNAL(changed()), this, SLOT(infinite_graph_update()));
    connect(ui->actionMultiple_Modes, SIGNAL(changed()), this, SLOT(infinite_multiple_modes()));

    infinite_graph_update();

}

void QM_Pot_GUI::finite_init()
{
    /*
    connect(ui->finite_width_slide, &QSlider::valueChanged, ui->finite_width_spin, &QDoubleSpinBox::setValue);
    connect(ui->finite_depth_slide, &QSlider::valueChanged, ui->finite_depth_spin, &QDoubleSpinBox::setValue);
    connect(ui->finite_mass_slide, &QSlider::valueChanged, ui->finite_mass_spin, &QDoubleSpinBox::setValue);

    connect(ui->finite_width_spin, static_cast<void (QSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), ui->finite_width_slide, &QSlider::setValue);
    connect(ui->finite_depth_spin, static_cast<void (QSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), ui->finite_depth_slide, &QSlider::setValue);
    connect(ui->finite_mass_spin, static_cast<void (QSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), ui->finite_mass_slide, &QSlider::setValue);
    */
    connect(ui->finite_width_slide, SIGNAL(valueChanged(int)), ui->finite_width_spin, SLOT(setValue(int)));
    connect(ui->finite_width_spin, SIGNAL(valueChanged(int)), ui->finite_width_slide, SLOT(setValue(int)));
    connect(ui->finite_depth_slide, SIGNAL(valueChanged(int)), ui->finite_depth_spin, SLOT(setValue(int)));
    connect(ui->finite_depth_spin, SIGNAL(valueChanged(int)), ui->finite_depth_slide, SLOT(setValue(int)));
    connect(ui->finite_mass_slide, SIGNAL(valueChanged(int)), ui->finite_mass_spin, SLOT(setValue(int)));
    connect(ui->finite_mass_spin, SIGNAL(valueChanged(int)), ui->finite_mass_slide, SLOT(setValue(int)));

    QVector<double> x(101), y(101); // initialize with entries 0..100
    for (int i=0; i<101; ++i)
    {
      x[i] = i/25.0 - 2; // x goes from -1 to 1
      y[i] = pow(sin(x[i]),2); // let's plot a quadratic function
    }
    ui->finite_graph->addGraph();
    ui->finite_graph->graph(0)->addData(x,y);
    ui->finite_graph->xAxis->setRange(-40,40);
    ui->finite_graph->yAxis->setRange(-1,1);

    connect(ui->finite_width_spin, SIGNAL(valueChanged(QString)), this, SLOT(finite_graph_update()));
    connect(ui->finite_depth_spin, SIGNAL(valueChanged(QString)), this, SLOT(finite_graph_update()));
    connect(ui->finite_mass_spin, SIGNAL(valueChanged(QString)), this, SLOT(finite_graph_update()));

    finite_graph_update();

}

void QM_Pot_GUI::harmonic_init(){
    QVector<double> x(1001), y(1001); // initialize with entries 0..100
    double h[250];
    Hermite(0,h);
    for (int i=0; i<1001; ++i)
    {
      x[i] = double(i)/100.0 - 5;
      y[i] = Hermite_poly(x[i], 1., 2., 0, h);
    }
    ui->harmonic_graph->addGraph();
    ui->harmonic_graph->graph(0)->addData(x,y);
    ui->harmonic_graph->xAxis->setRange(-10,10);
    ui->harmonic_graph->yAxis->setRange(-1,1);
}

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
///public slots


void QM_Pot_GUI::infinite_graph_update()
{
    double width = double(ui->infinite_width_spin->value());
    int mode = ui->infinite_mode_spin->value();

    double number_points = 50*width;

    QVector<double> x(number_points+1), y(number_points+1, sqrt(2/width)), inf_x(8), inf_y(8);
    inf_x << -25. << -width/2.-1e-6 << -width/2. << 0.-1e-5 << 0. << width/2. << width/2.+1e-6 << 25.;
    inf_y << 1.1 << 1.1 << 0.5 << 0.5 << 0.5 << 0.5 << 1.1 << 1.1;

    if(ui->actionProbability->isChecked()){
        for(int i = 0; i <= number_points; ++i){
            x[i] = (double(i)/number_points-0.5)*width;
            if(mode%2){
                y[i] = 2/width*pow(cos(m_Pi*mode/width*x[i]), 2);
            }
            else y[i] = 2/width*pow(sin(m_Pi*mode/width*x[i]), 2);
        }
    }
    else{
        for(int i = 0; i <= number_points; ++i){
            x[i] = (double(i)/number_points-0.5)*width;
            if(mode%2){
                y[i] *= cos(m_Pi*mode/width*x[i]);
            }
            else y[i] *= sin(m_Pi*mode/width*x[i]);
        }
    }

    ui->infinite_graph->graph(0)->setData(x,y);

    ui->infinite_graph->graph(1)->setData(inf_x, inf_y);
    ui->infinite_graph->replot();
}

void QM_Pot_GUI::infinite_multiple_modes()
{
    if(ui->actionMultiple_Modes->isChecked()){
        ui->infinite_mode->setEnabled(0);
        ui->infinite_mode_list->setEnabled(1);
    }
    if(!ui->actionMultiple_Modes->isChecked()){
        ui->infinite_mode->setEnabled(1);
        ui->infinite_mode_list->setEnabled(0);
    }
}




void QM_Pot_GUI::finite_graph_update()
{
    double width = double(ui->finite_width_spin->value());
    double depth = double(ui->finite_depth_spin->value())/10.;
    double mass = double(ui->finite_mass_spin->value())/1000.;


    QVector<double> mode = bisection(width, mass*depth);

    //for(int i = 0; i < mode.size(); ++i){
    unsigned int i = 1;
    double k = mode[i];
        double kappa = sqrt(2*mass*depth-k*k);

        finite_graph_plot(k, kappa, i, width );
    //}


}

void QM_Pot_GUI::finite_graph_plot(double k, double kappa, unsigned int mode, double width ){
    double position, number_points = 1.5*ui->finite_graph->width();
    QVector<double> x(number_points+1), y(number_points+1);

    if(mode%2){
        double A = 1./sqrt( pow(sin(k*width/2.),2)/kappa + 1./2.*( width - 1/k*sin(k*width) ) );
        double B = A*sin(k*width/2.)*exp(kappa*width/2.);
        double AA = A*A, BB = B*B;

        if(ui->actionProbability->isChecked()){
            for(int i = 0; i <= number_points; ++i){
                position = -75. + 150./double(number_points)*double(i);
                x[i] = position;
                if(position < -width/2.){
                    y[i] = BB*exp(2.*kappa*position);
                }
                if(position <= width/2. && position >= -width/2.){
                    y[i] = AA*pow(sin(k*position),2);
                }
                if(position > width/2.){
                    y[i] = BB*exp(-2.*kappa*position);
                }
            }
        }

        else{
            for(int i = 0; i <= number_points; ++i){
                position = -75. + 150./double(number_points)*double(i);
                x[i] = position;
                if(position < -width/2.){
                    y[i] = -B*exp(kappa*position);
                }
                if(position <= width/2. && position >= -width/2.){
                    y[i] = A*sin(k*position);
                }
                if(position > width/2.){
                    y[i] = B*exp(-kappa*position);
                }
            }
        }

    }

    else{
        double A = 1./sqrt( pow(cos(k*width/2.),2)/kappa + 1./2.*( width + 1/k*sin(k*width) ) );
        double B = A*cos(k*width/2.)*exp(kappa*width/2.);
        double AA = A*A, BB = B*B;
        if(ui->actionProbability->isChecked()){
            for(int i = 0; i <= number_points; ++i){
                position = -75. + 150./double(number_points)*double(i);
                x[i] = position;
                if(position < -width/2.){
                    y[i] = BB*exp(2.*kappa*position);
                }
                if(position <= width/2. && position >= -width/2.){
                    y[i] = AA*pow(cos(k*position),2);
                }
                if(position > width/2.){
                    y[i] = BB*exp(-2.*kappa*position);
                }
            }
        }


        else{
            for(int i = 0; i <= number_points; ++i){
                position = -75. + 150./double(number_points)*double(i);
                x[i] = position;
                if(position < -width/2.){
                    y[i] = B*exp(kappa*position);
                }
                if(position <= width/2. && position >= -width/2.){
                    y[i] = A*cos(k*position);
                }
                if(position > width/2.){
                    y[i] = B*exp(-kappa*position);
                }
            }
        }

    }


    ui->finite_graph->graph(0)->setData(x, y);
    ui->finite_graph->replot();
}

void QM_Pot_GUI::print_PDF()
{

}


void QM_Pot_GUI::print_PNG()
{

}

