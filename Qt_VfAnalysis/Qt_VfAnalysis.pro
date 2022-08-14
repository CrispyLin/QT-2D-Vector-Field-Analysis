QT       += core gui openglwidgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    Analysis/CalSepAttpPts.cpp \
    Analysis/DirGraphOps.cpp \
    Analysis/FixedPtExtraction.cpp \
    Analysis/MCG&ECG/BuildECG.cpp \
    Analysis/MCG&ECG/BuildMCG.cpp \
    Analysis/MorseDecomp.cpp \
    Analysis/PeriodicOrbitExtraction.cpp \
    BuildGeom/BuildCorners.cpp \
    BuildGeom/BuildEdges.cpp \
    BuildGeom/Geometry.cpp \
    Analysis/CalSepAttpPts.cpp \
    ConleyIndex.cpp \
    Contour_2D.cpp \
    EdgeSamplePts.cpp \
    FileLoader/PlyLoader.cpp \
    FileLoader/vfply_io.cpp \
    Numerical.cpp \
    Others/TraceBall.cpp \
    Others/common_routines.cpp \
    RegionTauMap.cpp \
    StreamlineCalculate/EvenStreamlines.cpp \
    StreamlineCalculate/Streamline_cal.cpp \
    main.cpp \
    MainWindow.cpp \
    OpenGLWindow.cpp \
    utility_functions.cpp

HEADERS += \
    Analysis/MorseDecomp.h \
    Analysis/PeriodicOrbitExtraction.h \
    BuildGeom/Geometry.h \
    BuildGeom/Point.h \
    ConleyIndex.h \
    Edit/RotateReflectField.h \
    ExternalDependencies/icMatrix.h \
    ExternalDependencies/icVector.h \
    ExternalDependencies/nr.h \
    ExternalDependencies/nrtypes_nr.h \
    ExternalDependencies/nrutil_nr.h \
    FileLoader/PlyLoader.h \
    FileLoader/vfply_io.h \
    GL_LIB/glut.h \
    Graph2D.h \
    Numerical.h \
    Others/TraceBall.h \
    Others/common_routines.h \
    RegionTauMap.h \
    Resource.h \
    RotateReflectField.h \
    StreamlineCalculate/EvenStreamlines.h \
    VField.h \
    MainWindow.h \
    OpenGLWindow.h \
    Predefined.h \
    GL_LIB/glew.h

FORMS += \
    MainWindow.ui

LIBS += -lglu32 -lopengl32


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    Datasets/bunny1.ply \
    Datasets/simple_fld2.ply \
    Datasets/sphere1.ply \
    Datasets/sphere2.ply \
    Datasets/torus1.ply \
    Datasets/torus2.ply

