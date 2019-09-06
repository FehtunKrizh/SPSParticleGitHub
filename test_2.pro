TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

HEADERS += \
    MathVector3d.h \
    Particle.h \
    Box.h \
    Matrix.h \
    Array3DBox.h \
    Press.h
CONFIG += -O3
QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS+= -fopenmp


LIBS+= -fopenmp -lgomp
