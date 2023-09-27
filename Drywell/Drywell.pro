TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

DEFINES += _arma

SOURCES += \
    ../../Utilities/Distribution.cpp \
    ../../Utilities/Matrix.cpp \
    ../../Utilities/Matrix_arma.cpp \
    ../../Utilities/Matrix_arma_sp.cpp \
    ../../Utilities/QuickSort.cpp \
    ../../Utilities/Utilities.cpp \
    ../../Utilities/Vector.cpp \
    ../../Utilities/Vector_arma.cpp \
    cell.cpp \
    grid.cpp \
    main.cpp

HEADERS += \
    ../../Utilities/BTC.h \
    ../../Utilities/BTC.hpp \
    ../../Utilities/BTCSet.h \
    ../../Utilities/BTCSet.hpp \
    ../../Utilities/Distribution.h \
    ../../Utilities/Matrix.h \
    ../../Utilities/Matrix_arma.h \
    ../../Utilities/Matrix_arma_sp.h \
    ../../Utilities/QuickSort.h \
    ../../Utilities/Utilities.h \
    ../../Utilities/Vector.h \
    ../../Utilities/Vector_arma.h \
    cell.h \
    grid.h


INCLUDEPATH += ../../Utilities/

win32 {
    LAPACK_INCLUDE = $$PWD/include
    #64 bits build
    contains(QMAKE_TARGET.arch, x86_64) {
        #debug
        CONFIG(debug, debug|release) {
            LAPACK_LIB_DIR = $$PWD/libs/lapack-blas_lib_win64/debug
            LIBS +=  -L$${LAPACK_LIB_DIR} -llapack_win64_MTd \
                    -lblas_win64_MTd
        }
        #release
        CONFIG(release, debug|release) {
            LAPACK_LIB_DIR = $$PWD/libs/lapack-blas_lib_win64/release
            LIBS +=  -L$${LAPACK_LIB_DIR} -llapack_win64_MT \
                    -lblas_win64_MT
        }
    }

    INCLUDEPATH += $${LAPACK_INCLUDE}

    DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS

}

linux {
    #sudo apt-get install libblas-dev liblapack-dev
     DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS
     LIBS += -larmadillo -llapack -lblas -lgsl
}

CONFIG(debug, debug|release) {
    message(Building in debug mode)
    !macx: QMAKE_CXXFLAGS *= "-Xpreprocessor -fopenmp"
    !macx: QMAKE_LFLAGS +=  -fopenmp
    !macx: LIBS += -lgomp -lpthread
    LIBS += -lpthread
    DEFINES += NO_OPENMP DEBUG

} else {
    message(Building in release mode)
    !macx: QMAKE_CXXFLAGS *= "-Xpreprocessor -fopenmp"
    !macx: QMAKE_LFLAGS +=  -fopenmp
    # QMAKE_CFLAGS+=-pg
    # QMAKE_CXXFLAGS+=-pg
    # QMAKE_LFLAGS+=-pg
    # macx: DEFINES += NO_OPENMP
    ! macx: LIBS += -lgomp -lpthread
    macx: LIBS += -lpthread
}
