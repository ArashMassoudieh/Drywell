TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

DEFINES += _arma
DEFINES += use_VTK
CONFIG += use_VTK
#VTKBUILDPATH = /home/arash/Projects/VTK-9.2.0/
VTKBUILDPATH = /home/arash/Projects/VTK/VTK-build
VTKHEADERPATH = /home/arash/Projects/VTK
VTK_V = -9.0

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


use_VTK {
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkChartsCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonColor$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonComputationalGeometry$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonDataModel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonExecutionModel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonMath$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonMisc$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonSystem$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonTransforms$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkDICOMParser$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkexpat$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersAMR$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersExtraction$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersFlowPaths$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeneral$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeneric$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeometry$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersHybrid$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersHyperTree$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersImaging$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersModeling$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersParallel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersParallelImaging$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersPoints$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersProgrammable$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSelection$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSMP$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSources$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersStatistics$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersTexture$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersTopology$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersVerdict$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkfreetype$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkGeovisCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkgl2ps$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkglew$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkhdf5$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkhdf5_hl$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingColor$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingFourier$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingGeneral$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingHybrid$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingMath$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingMorphological$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingSources$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingStatistics$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingStencil$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInfovisCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInfovisLayout$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionImage$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionStyle$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionWidgets$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOAMR$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOEnSight$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOExodus$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOGeometry$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOImage$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOImport$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOInfovis$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOLegacy$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOLSDyna$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOMINC$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOMovie$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIONetCDF$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOParallel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOParallelXML$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOPLY$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOSQL$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOTecplotTable$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOVideo$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOXML$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOXMLParser$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkjpeg$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkjsoncpp$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtklibharu$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtklibxml2$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtklz4$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkmetaio$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkParallelCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkpng$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingAnnotation$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingContext2D$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingFreeType$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingGL2PSOpenGL2$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingImage$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingLabel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingLOD$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingOpenGL2$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingVolume$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingVolumeOpenGL2$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtksqlite$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtksys$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtktiff$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkverdict$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsContext2D$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsInfovis$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkzlib$$VTK_V
    LIBS += -L"/usr/local/lib/ -lsuperlu.so"

    #VTK Include files
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/Core
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Color
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/Color
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/DataModel
    INCLUDEPATH +=$${VTKBUILDPATH}/Utilities/KWIML
    INCLUDEPATH +=$${VTKHEADERPATH}/Utilities/KWIML
    INCLUDEPATH +=$${VTKHEADERPATH}/Rendering/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Rendering/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Core
    INCLUDEPATH +=$${VTKHEADERPATH}/Charts/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Charts/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Filters/General
    INCLUDEPATH +=$${VTKBUILDPATH}/Rendering/Context2D
    INCLUDEPATH +=$${VTKHEADERPATH}/Rendering/Context2D
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/DataModel
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Math
    INCLUDEPATH +=$${VTKHEADERPATH}/Views/Context2D
    INCLUDEPATH +=$${VTKBUILDPATH}/Views/Context2D
    INCLUDEPATH +=$${VTKBUILDPATH}/Views/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Interaction/Widgets
    INCLUDEPATH +=$${VTKHEADERPATH}/Views/Core
    INCLUDEPATH +=$${VTKHEADERPATH}/Interaction/Style
    INCLUDEPATH +=$${VTKBUILDPATH}/Interaction/Style
    INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Modeling
    INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Modeling
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/ExecutionModel
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/ExecutionModel
    INCLUDEPATH +=$${VTKHEADERPATH}/Interaction/Widgets/
    INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Core/
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Misc/
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/Misc
    INCLUDEPATH +=$${VTKHEADERPATH}/IO/XML/
    INCLUDEPATH +=$${VTKBUILDPATH}/IO/XML
    INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Sources
    INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Sources
    INCLUDEPATH +=$${VTKHEADERPATH}/Filters/General
    INCLUDEPATH +=$${VTKHEADERPATH}/IO/Image
    INCLUDEPATH +=$${VTKBUILDPATH}/IO/Image
    INCLUDEPATH +=$${VTKHEADERPATH}/Imaging/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Imaging/Core

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
