# QT-2D-Vector-Field-Analysis
The original code is retrived from [`2D Vector Field Analysis`](http://www2.cs.uh.edu/~chengu/Publications/MorseDecomp/vfAnalysis2D.html). It contains a series of great papers about 2D vector field analysis using Morse Decompsition. However, the code they provided has to run in VS2019 and was written many years ago. 
My goal is to migrate the code from VS2019 to QT in order to fully understand the Morse Decompsition and 2D vector field analysis. Also, I want a better UI. This work should be beneficial to my own visualization research.

[`QT Creator`](https://www.qt.io/product/development-tools): is a mandatory IDE. It can be downloaded from QT company's office website. It's open source and free.
Note, I am currently using Qt 6.3.1. If you do not see the right version, you may need to configure my code to work for the lastest version.

#How to run (migrating is not finished)
1. Open the .pro file using QT creator (if the version is different, you need to configure the code to your setting)
2. Either click run or debug.


#Issues
1. Cannot use glut.h.
2. Since glew.h has to be defined before gl.h, everytime you build (newly) the project, you need to modify the ui_MainWindow.h as follow:
    Move 
    '#include <VectorFieldWindow.h>' 
    to the top of the include section (anywhere before the QOpenGLWidget).
