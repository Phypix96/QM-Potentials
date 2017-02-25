#ifndef QM_POT_GUI_H
#define QM_POT_GUI_H

#include <QMainWindow>
#include <vector>

//#include <QtDataVisualization>

//using namespace QtDataVisualization;

namespace Ui {

QVector<double> bisection(double);

class QM_Pot_GUI;
}

class QM_Pot_GUI : public QMainWindow
{
    Q_OBJECT

public:
    explicit QM_Pot_GUI(QWidget *parent = 0);
    ~QM_Pot_GUI();

public slots:
    void infinite_graph_update();
    void infinite_multiple_modes();
    void finite_graph_update();
    void finite_graph_plot(double, double, unsigned int, double);
    void print_PDF();
    void print_PNG();

private:
    Ui::QM_Pot_GUI *ui;
    void infinite_init();
    void finite_init();
    void harmonic_init();
};

#endif // QM_POT_GUI_H
