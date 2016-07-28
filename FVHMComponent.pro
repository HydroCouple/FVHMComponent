#-------------------------------------------------
#
# Project created by QtCreator 2016-06-06T10:49:44
#
#-------------------------------------------------

TEMPLATE = lib
TARGET = FVHMComponent
QT       -= gui

DEFINES += FVHMCOMPONENT_LIBRARY

INCLUDEPATH += .\
               ./include \
               ./../HydroCouple/include \
               ./../HydroCoupleSDK/include

PRECOMPILED_HEADER = ./include/stdafx.h

HEADERS += ./include/stdafx.h \
           ./include/fvhmcomponent.h \
           ./include/fvhmcomponent_global.h \
           ./include/fvhmcomponentinfo.h \
           ./include/core/consts.h \
           ./include/core/datetime.h \
           ./include/core/enums.h \
           ./include/core/error.h \
           ./include/core/exfil.h \
           ./include/core/findroot.h \
           ./include/core/funcs.h \
           ./include/core/globals.h \
           ./include/core/hash.h \
           ./include/core/headers.h \
           ./include/core/infil.h \
           ./include/core/keywords.h \
           ./include/core/lid.h \
           ./include/core/macros.h \
           ./include/core/mathexpr.h \
           ./include/core/mempool.h \
           ./include/core/objects.h \
           ./include/core/odesolve.h \
           ./include/core/swmm5.h \
           ./include/core/text.h \
           ./include/core/dataexchangecache.h \
           ./include/swmmobjectitems.h \
           ./include/swmmtimeseriesexchangeitems.h \
           ./include/fvhmobjectitems.h \
           ./include/fvhmtimeseriesexchangeitems.h

SOURCES += ./src/stdafx.cpp \
           ./src/fvhmcomponent.cpp \
           ./src/fvhmcomponentinfo.cpp \
           ./src/core/climate.c \
           ./src/core/controls.c \
           ./src/core/culvert.c \
           ./src/core/datetime.c \
           ./src/core/dwflow.c \
           ./src/core/dynwave.c \
           ./src/core/error.c \
           ./src/core/exfil.c \
          ./src/core/findroot.c \
          ./src/core/flowrout.c \
          ./src/core/forcmain.c \
          ./src/core/gage.c \
          ./src/core/gwater.c \
          ./src/core/hash.c \
          ./src/core/hotstart.c \
          ./src/core/iface.c \
          ./src/core/infil.c \
          ./src/core/inflow.c \
          ./src/core/input.c \
          ./src/core/inputrpt.c \
          ./src/core/keywords.c \
          ./src/core/kinwave.c \
          ./src/core/landuse.c \
          ./src/core/lid.c \
          ./src/core/lidproc.c \
          ./src/core/link.c \
          ./src/core/massbal.c \
          ./src/core/mathexpr.c \
          ./src/core/mempool.c \
          ./src/core/node.c \
          ./src/core/odesolve.c \
          ./src/core/output.c \
          ./src/core/project.c \
          ./src/core/qualrout.c \
          ./src/core/rain.c \
          ./src/core/rdii.c \
          ./src/core/report.c \
          ./src/core/roadway.c \
          ./src/core/routing.c \
          ./src/core/runoff.c \
          ./src/core/shape.c \
          ./src/core/snow.c \
          ./src/core/stats.c \
          ./src/core/statsrpt.c \
          ./src/core/subcatch.c \
          ./src/core/surfqual.c \
          ./src/core/swmm5.c \
          ./src/core/table.c \
          ./src/core/toposort.c \
          ./src/core/transect.c \
          ./src/core/treatmnt.c \
          ./src/core/xsect.c \
          ./src/core/dataexchangecache.cpp \
          ./src/fvhmtimeseriesexchangeitems.cpp



macx{

    INCLUDEPATH += /usr/local/include \
                   /usr/local/include/libiomp

#    QMAKE_CC = clang-omp
#    QMAKE_CXX = clang-omp++
#    QMAKE_LINK = $$QMAKE_CXX

#    QMAKE_CFLAGS = -fopenmp
#    QMAKE_LFLAGS = -fopenmp
#    QMAKE_CXXFLAGS = -fopenmp
#    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
#    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS

    LIBS += -L/usr/local/lib/ -liomp5
 }

CONFIG(debug, debug|release) {

   DESTDIR = ./build/debug
   OBJECTS_DIR = $$DESTDIR/.obj
   MOC_DIR = $$DESTDIR/.moc
   RCC_DIR = $$DESTDIR/.qrc
   UI_DIR = $$DESTDIR/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK.1.0.0
   }
}

CONFIG(release, debug|release) {

    DESTDIR = lib
    RELEASE_EXTRAS = ./build/release
    OBJECTS_DIR = $$RELEASE_EXTRAS/.obj
    MOC_DIR = $$RELEASE_EXTRAS/.moc
    RCC_DIR = $$RELEASE_EXTRAS/.qrc
    UI_DIR = $$RELEASE_EXTRAS/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/lib -lHydroCoupleSDK.1.0.0
   }
}
