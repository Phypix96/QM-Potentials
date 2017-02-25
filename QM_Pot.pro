#-------------------------------------------------
#
# Project created by QtCreator 2016-10-30T08:25:28
#
#-------------------------------------------------

QT       += core gui


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport


TARGET = QM_Pot
TEMPLATE = app


SOURCES += main.cpp\
        qm_pot_gui.cpp \
    qcustomplot.cpp \
    qm_potentials.cpp

HEADERS  += qm_pot_gui.h \
    qcustomplot.h

FORMS    += qm_pot_gui.ui
