#include "qm_pot_gui.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QM_Pot_GUI w;
    w.show();

    return a.exec();
}
