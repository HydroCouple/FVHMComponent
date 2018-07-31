#Author Caleb Amoa Buahin
#Email caleb.buahin@gmail.com
#Date 2015-2018
#License GNU Lesser General Public License (see <http: //www.gnu.org/licenses/> for details).


TARGET = FVHMComponent
QT       -= gui

DEFINES += FVHMCOMPONENT_LIBRARY
DEFINES += UTAH_CHPC
DEFINES += USE_OPENMP
DEFINES += USE_MPI
DEFINES += USE_HYPRE_OPENMP
#DEFINES += USE_SUITESPARSE

CONFIG += c++11
CONFIG += debug_and_release

contains(DEFINES,FVHMCOMPONENT_LIBRARY){
  TEMPLATE = lib
  message("Compiling as library")
} else {
  TEMPLATE = app
  CONFIG-=app_bundle
  message("Compiling as application")
}


INCLUDEPATH += .\
               ./include \
               ./../HydroCouple/include \
               ./../HydroCoupleSDK/include

PRECOMPILED_HEADER = ./include/stdafx.h

HEADERS += ./include/stdafx.h \
           ./include/fvhmcomponent.h \
           ./include/fvhmcomponent_global.h \
           ./include/fvhmcomponentinfo.h \
           ./include/controlvolume.h \
           ./include/sparsemat.h \
           ./include/inletoutletflowbc.h \
           ./include/wsebc.h \
           ./include/sourcesinkbc.h \
           ./include/edgebc.h \
           ./include/qrsolve.h \
           ./include/macros.h \
           ./include/outletwseslope.h \
           ./include/precipbc.h \
           ./include/boundarycondition.h \
           ./include/inflowinput.h \
           ./include/cvwseoutput.h \
           ./include/cvdepthoutput.h \
           ./include/initialwsebc.h \
           ./include/criticaldepthoutflowbc.h \
           ./include/edgefluxesio.h


SOURCES += ./src/stdafx.cpp \
           ./src/fvhmcomponent.cpp \
           ./src/fvhmcomponentinfo.cpp \
           ./src/main.cpp \
           ./src/controlvolume.cpp \
           ./src/sparsemat.cpp \
           ./src/inletoutletflowbc.cpp \
           ./src/wsebc.cpp \
           ./src/fvhmcomponentcompute.cpp \
           ./src/fvhmcomponentcompgeom.cpp \
           ./src/edgebc.cpp \
           ./src/fvhmcomponentio.cpp \
           ./src/qrsolve.cpp \
           ./src/outletwseslope.cpp \
           ./src/precipbc.cpp \
           ./src/inflowinput.cpp \
           ./src/cvwseoutput.cpp \
           ./src/cvdepthoutput.cpp \
           ./src/initialwsebc.cpp \
           ./src/criticaldepthoutflowbc.cpp \
           ./src/edgefluxesio.cpp


macx{

    INCLUDEPATH += /usr/local/include \
                   ./../HYPRE/hypre-2.11.2/include \
                   /usr/local/include/libiomp

    LIBS += -L/usr/local/lib -lnetcdf-cxx4 \
            -L./../HYPRE/hypre-2.11.2/lib -lHYPRE

    contains(DEFINES,USE_OPENMP){

        QMAKE_CC = /usr/local/opt/llvm/bin/clang
        QMAKE_CXX = /usr/local/opt/llvm/bin/clang++
        QMAKE_LINK = /usr/local/opt/llvm/bin/clang++

        QMAKE_CFLAGS+= -fopenmp
        QMAKE_LFLAGS+= -fopenmp
        QMAKE_CXXFLAGS+= -fopenmp

        INCLUDEPATH += /usr/local/opt/llvm/lib/clang/5.0.0/include
        LIBS += -L /usr/local/opt/llvm/lib -lomp


      message("OpenMP enabled")
     } else {
      message("OpenMP disabled")
     }


    contains(DEFINES,USE_HYPRE_OPENMP){

        QMAKE_CC = /usr/local/opt/llvm/bin/clang
        QMAKE_CXX = /usr/local/opt/llvm/bin/clang++
        QMAKE_LINK = /usr/local/opt/llvm/bin/clang++

        QMAKE_CFLAGS+= -fopenmp
        QMAKE_LFLAGS+= -fopenmp
        QMAKE_CXXFLAGS+= -fopenmp

        INCLUDEPATH += /usr/local/opt/llvm/lib/clang/5.0.0/include
        LIBS += -L /usr/local/opt/llvm/lib -lomp


      message("OpenMP enabled")
     } else {
      message("OpenMP disabled")
     }

    contains(DEFINES,USE_MPI){

        QMAKE_CC = /usr/local/bin/mpicc
        QMAKE_CXX = /usr/local/bin/mpicxx
        QMAKE_LINK = /usr/local/bin/mpicxx

        QMAKE_CFLAGS += $$system(/usr/local/bin/mpicc --showme:compile)
        QMAKE_CXXFLAGS += $$system(/usr/local/bin/mpic++ --showme:compile)
        QMAKE_LFLAGS += $$system(/usr/local/bin/mpic++ --showme:link)

        LIBS += -L/usr/local/lib/ -lmpi

      message("MPI enabled")
     } else {
      message("MPI disabled")
     }
}

linux{

INCLUDEPATH += /usr/include \
               ../gdal/include

LIBS += -L/usr/lib/ogdi -lgdal \
        -L../gdal/lib -lgdal

    contains(DEFINES,UTAH_CHPC){

         INCLUDEPATH += /uufs/chpc.utah.edu/sys/installdir/hdf5/1.8.17-c7/include \
                        /uufs/chpc.utah.edu/sys/installdir/netcdf-c/4.4.1/include \
                        ../netcdf-cxx4-4.3.0/installdir/include \
                        ./../HYPRE/hypre-2.11.2/include

         LIBS += -L/uufs/chpc.utah.edu/sys/installdir/hdf5/1.8.17-c7/lib -lhdf5 \
                 -L/uufs/chpc.utah.edu/sys/installdir/netcdf-cxx/4.3.0-c7/lib -lnetcdf_c++4 \
                 -L./../HYPRE/hypre-2.11.2/lib -lHYPRE

         message("Compiling on CHPC")
     }


    contains(DEFINES,USE_MPI){

        QMAKE_CC = mpicc
        QMAKE_CXX = mpic++
        QMAKE_LINK = mpic++

        QMAKE_CFLAGS += $$system(/usr/local/bin/mpicc --showme:compile)
        QMAKE_CXXFLAGS += $$system(/usr/local/bin/mpic++ --showme:compile)
        QMAKE_LFLAGS += $$system(/usr/local/bin/mpic++ --showme:link)

        LIBS += -L/usr/local/lib/ -lmpi

      message("MPI enabled")
     } else {
      message("MPI disabled")
     }

    contains(DEFINES,USE_OPENMP){

    QMAKE_CFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
    QMAKE_CXXFLAGS += -fopenmp

    LIBS += -L/usr/lib/x86_64-linux-gnu -lgomp

      message("OpenMP enabled")
     } else {
      message("OpenMP disabled")
     }
}

CONFIG(debug, debug|release) {

   DESTDIR = ./build/debug
   OBJECTS_DIR = $$DESTDIR/.obj
   MOC_DIR = $$DESTDIR/.moc
   RCC_DIR = $$DESTDIR/.qrc
   UI_DIR = $$DESTDIR/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK.1.0.0
    QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK* ./build/debug/";
     }

   linux{
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK.so.1.0.0
    QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK* ./build/debug/";
     }

   win32{
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK1
    QMAKE_POST_LINK += "copy ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK* ./build/debug/";
     }

   message("Debug mode...")
}

CONFIG(release, debug|release) {

    RELEASE_EXTRAS = ./build/release
    OBJECTS_DIR = $$RELEASE_EXTRAS/.obj
    MOC_DIR = $$RELEASE_EXTRAS/.moc
    RCC_DIR = $$RELEASE_EXTRAS/.qrc
    UI_DIR = $$RELEASE_EXTRAS/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/lib/macx -lHydroCoupleSDK.1.0.0
     }

   linux{
    LIBS += -L./../HydroCoupleSDK/lib/linux -lHydroCoupleSDK.1.0.0
     }

   win32{
    LIBS += -L./../HydroCoupleSDK/lib/win32 -lHydroCoupleSDK1
     }


     contains(DEFINES,FVHMCOMPONENT_LIBRARY){
         #MacOS
         macx{
             DESTDIR = lib/macx
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/macx/*HydroCoupleSDK* ./lib/macx";
          }

         #Linux
         linux{
             DESTDIR = lib/linux
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/linux/*HydroCoupleSDK* ./lib/linux";
          }

         #Windows
         win32{
             DESTDIR = lib/win32
             QMAKE_POST_LINK += "xcopy ./../HydroCoupleSDK/lib/win32/*HydroCoupleSDK* ./lib/win32";
          }
     } else {
         #MacOS
         macx{
             DESTDIR = bin/macx
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/macx/*HydroCoupleSDK* ./bin/macx";
          }

         #Linux
         linux{
             DESTDIR = bin/linux
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/linux/*HydroCoupleSDK* ./bin/linux";
          }

         #Windows
         win32{
             DESTDIR = bin/win32
             QMAKE_POST_LINK += "xcopy ./../HydroCoupleSDK/lib/win32/*HydroCoupleSDK* ./bin/win32";
          }
     }

   message("Release mode...")
}


contains(DEFINES,USE_SUITESPARSE){

INCLUDEPATH += ../SuiteSparse/include

LIBS += -L../SuiteSparse/lib -lspqr \
        -L../SuiteSparse/lib -lsuitesparseconfig \
        -L../SuiteSparse/lib -lcholmod \
        -L../SuiteSparse/lib -lamd \
        -L../SuiteSparse/lib -lcolamd
}
