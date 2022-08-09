/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created by: Qt User Interface Compiler version 6.3.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <OpenGLWindow.h>
#include <QtCore/QVariant>
#include <QtOpenGLWidgets/QOpenGLWidget>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_21;
    OpenGLWindow *openGLWindow;
    QVBoxLayout *verticalLayout_16;
    QGroupBox *groupBox_12;
    QHBoxLayout *horizontalLayout_22;
    QOpenGLWidget *MCG_Window;
    QGroupBox *groupBox_13;
    QHBoxLayout *horizontalLayout_23;
    QOpenGLWidget *ECG_Window;
    QWidget *layoutWidget1;
    QHBoxLayout *horizontalLayout_24;
    QVBoxLayout *verticalLayout_15;
    QGroupBox *groupBox_10;
    QVBoxLayout *verticalLayout_9;
    QHBoxLayout *horizontalLayout_18;
    QLabel *label_17;
    QLineEdit *lineEdit_15;
    QPushButton *pushButton_6;
    QHBoxLayout *horizontalLayout_19;
    QLabel *label_18;
    QLineEdit *lineEdit_16;
    QPushButton *pushButton_7;
    QPushButton *pushButton_8;
    QGroupBox *groupBox_9;
    QVBoxLayout *verticalLayout_8;
    QCheckBox *checkBox_22;
    QHBoxLayout *horizontalLayout_14;
    QLabel *label_14;
    QLineEdit *lineEdit_12;
    QHBoxLayout *horizontalLayout_17;
    QHBoxLayout *horizontalLayout_15;
    QLabel *label_15;
    QLineEdit *lineEdit_13;
    QHBoxLayout *horizontalLayout_16;
    QLabel *label_16;
    QLineEdit *lineEdit_14;
    QCheckBox *checkBox_23;
    QPushButton *pushButton_4;
    QPushButton *pushButton_5;
    QHBoxLayout *horizontalLayout_20;
    QGroupBox *groupBox_7;
    QVBoxLayout *verticalLayout_5;
    QHBoxLayout *horizontalLayout_10;
    QLabel *label_10;
    QLineEdit *lineEdit_8;
    QHBoxLayout *horizontalLayout_11;
    QLabel *label_11;
    QLineEdit *lineEdit_9;
    QHBoxLayout *horizontalLayout_12;
    QLabel *label_12;
    QLineEdit *lineEdit_10;
    QHBoxLayout *horizontalLayout_13;
    QLabel *label_13;
    QLineEdit *lineEdit_11;
    QCheckBox *checkBox_18;
    QVBoxLayout *verticalLayout_14;
    QGroupBox *groupBox_8;
    QVBoxLayout *verticalLayout_12;
    QVBoxLayout *verticalLayout_7;
    QCheckBox *checkBox_19;
    QCheckBox *checkBox_20;
    QCheckBox *checkBox_21;
    QPushButton *pushButton_3;
    QGroupBox *groupBox_2;
    QVBoxLayout *verticalLayout_11;
    QHBoxLayout *horizontalLayout_6;
    QCheckBox *checkBox_11;
    QCheckBox *checkBox_13;
    QCheckBox *checkBox_14;
    QHBoxLayout *horizontalLayout_7;
    QCheckBox *checkBox_12;
    QCheckBox *checkBox_15;
    QCheckBox *checkBox_16;
    QCheckBox *checkBox_17;
    QGroupBox *groupBox_6;
    QVBoxLayout *verticalLayout_4;
    QHBoxLayout *horizontalLayout_9;
    QLabel *label_6;
    QLineEdit *lineEdit_6;
    QLabel *label_7;
    QHBoxLayout *horizontalLayout_8;
    QLabel *label_8;
    QLineEdit *lineEdit_7;
    QLabel *label_9;
    QGroupBox *groupBox_11;
    QVBoxLayout *verticalLayout_6;
    QTextEdit *debug_window;
    QVBoxLayout *verticalLayout_13;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout_3;
    QCheckBox *checkBox;
    QCheckBox *checkBox_2;
    QCheckBox *checkBox_3;
    QCheckBox *checkBox_4;
    QCheckBox *checkBox_6;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_3;
    QLineEdit *lineEdit_3;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_4;
    QLineEdit *lineEdit_4;
    QCheckBox *checkBox_7;
    QCheckBox *checkBox_8;
    QCheckBox *checkBox_9;
    QVBoxLayout *verticalLayout_10;
    QGroupBox *groupBox_3;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_5;
    QPushButton *pushButton;
    QGroupBox *groupBox_5;
    QVBoxLayout *verticalLayout;
    QRadioButton *radioButton;
    QRadioButton *radioButton_2;
    QRadioButton *radioButton_3;
    QGroupBox *Refine_Morse_Set_Layout;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QLineEdit *lineEdit;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QLineEdit *lineEdit_2;
    QPushButton *pushButton_2;
    QButtonGroup *IntergratorsBG;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(2283, 882);
        MainWindow->setMinimumSize(QSize(0, 0));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        centralwidget->setMinimumSize(QSize(0, 0));
        layoutWidget = new QWidget(centralwidget);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(20, 30, 1339, 813));
        horizontalLayout_21 = new QHBoxLayout(layoutWidget);
        horizontalLayout_21->setObjectName(QString::fromUtf8("horizontalLayout_21"));
        horizontalLayout_21->setContentsMargins(0, 0, 0, 0);
        openGLWindow = new OpenGLWindow(layoutWidget);
        openGLWindow->setObjectName(QString::fromUtf8("openGLWindow"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(80);
        sizePolicy.setVerticalStretch(80);
        sizePolicy.setHeightForWidth(openGLWindow->sizePolicy().hasHeightForWidth());
        openGLWindow->setSizePolicy(sizePolicy);
        openGLWindow->setMinimumSize(QSize(800, 800));
        openGLWindow->setMaximumSize(QSize(800, 800));

        horizontalLayout_21->addWidget(openGLWindow);

        verticalLayout_16 = new QVBoxLayout();
        verticalLayout_16->setObjectName(QString::fromUtf8("verticalLayout_16"));
        groupBox_12 = new QGroupBox(layoutWidget);
        groupBox_12->setObjectName(QString::fromUtf8("groupBox_12"));
        QSizePolicy sizePolicy1(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(groupBox_12->sizePolicy().hasHeightForWidth());
        groupBox_12->setSizePolicy(sizePolicy1);
        horizontalLayout_22 = new QHBoxLayout(groupBox_12);
        horizontalLayout_22->setObjectName(QString::fromUtf8("horizontalLayout_22"));
        MCG_Window = new QOpenGLWidget(groupBox_12);
        MCG_Window->setObjectName(QString::fromUtf8("MCG_Window"));
        QSizePolicy sizePolicy2(QSizePolicy::Maximum, QSizePolicy::Maximum);
        sizePolicy2.setHorizontalStretch(50);
        sizePolicy2.setVerticalStretch(50);
        sizePolicy2.setHeightForWidth(MCG_Window->sizePolicy().hasHeightForWidth());
        MCG_Window->setSizePolicy(sizePolicy2);
        MCG_Window->setMinimumSize(QSize(500, 350));
        MCG_Window->setMaximumSize(QSize(500, 350));

        horizontalLayout_22->addWidget(MCG_Window);


        verticalLayout_16->addWidget(groupBox_12);

        groupBox_13 = new QGroupBox(layoutWidget);
        groupBox_13->setObjectName(QString::fromUtf8("groupBox_13"));
        sizePolicy1.setHeightForWidth(groupBox_13->sizePolicy().hasHeightForWidth());
        groupBox_13->setSizePolicy(sizePolicy1);
        horizontalLayout_23 = new QHBoxLayout(groupBox_13);
        horizontalLayout_23->setObjectName(QString::fromUtf8("horizontalLayout_23"));
        ECG_Window = new QOpenGLWidget(groupBox_13);
        ECG_Window->setObjectName(QString::fromUtf8("ECG_Window"));
        QSizePolicy sizePolicy3(QSizePolicy::Maximum, QSizePolicy::Maximum);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(ECG_Window->sizePolicy().hasHeightForWidth());
        ECG_Window->setSizePolicy(sizePolicy3);
        ECG_Window->setMinimumSize(QSize(500, 350));
        ECG_Window->setMaximumSize(QSize(500, 350));

        horizontalLayout_23->addWidget(ECG_Window);


        verticalLayout_16->addWidget(groupBox_13);


        horizontalLayout_21->addLayout(verticalLayout_16);

        layoutWidget1 = new QWidget(centralwidget);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(1350, 30, 911, 811));
        horizontalLayout_24 = new QHBoxLayout(layoutWidget1);
        horizontalLayout_24->setObjectName(QString::fromUtf8("horizontalLayout_24"));
        horizontalLayout_24->setContentsMargins(0, 0, 0, 0);
        verticalLayout_15 = new QVBoxLayout();
        verticalLayout_15->setObjectName(QString::fromUtf8("verticalLayout_15"));
        groupBox_10 = new QGroupBox(layoutWidget1);
        groupBox_10->setObjectName(QString::fromUtf8("groupBox_10"));
        QSizePolicy sizePolicy4(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(groupBox_10->sizePolicy().hasHeightForWidth());
        groupBox_10->setSizePolicy(sizePolicy4);
        verticalLayout_9 = new QVBoxLayout(groupBox_10);
        verticalLayout_9->setObjectName(QString::fromUtf8("verticalLayout_9"));
        horizontalLayout_18 = new QHBoxLayout();
        horizontalLayout_18->setObjectName(QString::fromUtf8("horizontalLayout_18"));
        horizontalLayout_18->setSizeConstraint(QLayout::SetMinimumSize);
        label_17 = new QLabel(groupBox_10);
        label_17->setObjectName(QString::fromUtf8("label_17"));
        sizePolicy4.setHeightForWidth(label_17->sizePolicy().hasHeightForWidth());
        label_17->setSizePolicy(sizePolicy4);

        horizontalLayout_18->addWidget(label_17);

        lineEdit_15 = new QLineEdit(groupBox_10);
        lineEdit_15->setObjectName(QString::fromUtf8("lineEdit_15"));
        QSizePolicy sizePolicy5(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(lineEdit_15->sizePolicy().hasHeightForWidth());
        lineEdit_15->setSizePolicy(sizePolicy5);

        horizontalLayout_18->addWidget(lineEdit_15);

        pushButton_6 = new QPushButton(groupBox_10);
        pushButton_6->setObjectName(QString::fromUtf8("pushButton_6"));
        sizePolicy5.setHeightForWidth(pushButton_6->sizePolicy().hasHeightForWidth());
        pushButton_6->setSizePolicy(sizePolicy5);

        horizontalLayout_18->addWidget(pushButton_6);


        verticalLayout_9->addLayout(horizontalLayout_18);

        horizontalLayout_19 = new QHBoxLayout();
        horizontalLayout_19->setObjectName(QString::fromUtf8("horizontalLayout_19"));
        label_18 = new QLabel(groupBox_10);
        label_18->setObjectName(QString::fromUtf8("label_18"));
        sizePolicy4.setHeightForWidth(label_18->sizePolicy().hasHeightForWidth());
        label_18->setSizePolicy(sizePolicy4);

        horizontalLayout_19->addWidget(label_18);

        lineEdit_16 = new QLineEdit(groupBox_10);
        lineEdit_16->setObjectName(QString::fromUtf8("lineEdit_16"));
        sizePolicy5.setHeightForWidth(lineEdit_16->sizePolicy().hasHeightForWidth());
        lineEdit_16->setSizePolicy(sizePolicy5);

        horizontalLayout_19->addWidget(lineEdit_16);

        pushButton_7 = new QPushButton(groupBox_10);
        pushButton_7->setObjectName(QString::fromUtf8("pushButton_7"));
        sizePolicy5.setHeightForWidth(pushButton_7->sizePolicy().hasHeightForWidth());
        pushButton_7->setSizePolicy(sizePolicy5);

        horizontalLayout_19->addWidget(pushButton_7);


        verticalLayout_9->addLayout(horizontalLayout_19);

        pushButton_8 = new QPushButton(groupBox_10);
        pushButton_8->setObjectName(QString::fromUtf8("pushButton_8"));
        sizePolicy5.setHeightForWidth(pushButton_8->sizePolicy().hasHeightForWidth());
        pushButton_8->setSizePolicy(sizePolicy5);
        pushButton_8->setMaximumSize(QSize(1000, 16777215));

        verticalLayout_9->addWidget(pushButton_8);


        verticalLayout_15->addWidget(groupBox_10);

        groupBox_9 = new QGroupBox(layoutWidget1);
        groupBox_9->setObjectName(QString::fromUtf8("groupBox_9"));
        sizePolicy4.setHeightForWidth(groupBox_9->sizePolicy().hasHeightForWidth());
        groupBox_9->setSizePolicy(sizePolicy4);
        verticalLayout_8 = new QVBoxLayout(groupBox_9);
        verticalLayout_8->setObjectName(QString::fromUtf8("verticalLayout_8"));
        verticalLayout_8->setSizeConstraint(QLayout::SetMinimumSize);
        checkBox_22 = new QCheckBox(groupBox_9);
        checkBox_22->setObjectName(QString::fromUtf8("checkBox_22"));

        verticalLayout_8->addWidget(checkBox_22);

        horizontalLayout_14 = new QHBoxLayout();
        horizontalLayout_14->setObjectName(QString::fromUtf8("horizontalLayout_14"));
        label_14 = new QLabel(groupBox_9);
        label_14->setObjectName(QString::fromUtf8("label_14"));

        horizontalLayout_14->addWidget(label_14);

        lineEdit_12 = new QLineEdit(groupBox_9);
        lineEdit_12->setObjectName(QString::fromUtf8("lineEdit_12"));
        QSizePolicy sizePolicy6(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(lineEdit_12->sizePolicy().hasHeightForWidth());
        lineEdit_12->setSizePolicy(sizePolicy6);

        horizontalLayout_14->addWidget(lineEdit_12);


        verticalLayout_8->addLayout(horizontalLayout_14);

        horizontalLayout_17 = new QHBoxLayout();
        horizontalLayout_17->setObjectName(QString::fromUtf8("horizontalLayout_17"));
        horizontalLayout_15 = new QHBoxLayout();
        horizontalLayout_15->setObjectName(QString::fromUtf8("horizontalLayout_15"));
        label_15 = new QLabel(groupBox_9);
        label_15->setObjectName(QString::fromUtf8("label_15"));

        horizontalLayout_15->addWidget(label_15);

        lineEdit_13 = new QLineEdit(groupBox_9);
        lineEdit_13->setObjectName(QString::fromUtf8("lineEdit_13"));
        sizePolicy6.setHeightForWidth(lineEdit_13->sizePolicy().hasHeightForWidth());
        lineEdit_13->setSizePolicy(sizePolicy6);

        horizontalLayout_15->addWidget(lineEdit_13);


        horizontalLayout_17->addLayout(horizontalLayout_15);

        horizontalLayout_16 = new QHBoxLayout();
        horizontalLayout_16->setObjectName(QString::fromUtf8("horizontalLayout_16"));
        label_16 = new QLabel(groupBox_9);
        label_16->setObjectName(QString::fromUtf8("label_16"));

        horizontalLayout_16->addWidget(label_16);

        lineEdit_14 = new QLineEdit(groupBox_9);
        lineEdit_14->setObjectName(QString::fromUtf8("lineEdit_14"));
        sizePolicy6.setHeightForWidth(lineEdit_14->sizePolicy().hasHeightForWidth());
        lineEdit_14->setSizePolicy(sizePolicy6);

        horizontalLayout_16->addWidget(lineEdit_14);


        horizontalLayout_17->addLayout(horizontalLayout_16);


        verticalLayout_8->addLayout(horizontalLayout_17);

        checkBox_23 = new QCheckBox(groupBox_9);
        checkBox_23->setObjectName(QString::fromUtf8("checkBox_23"));

        verticalLayout_8->addWidget(checkBox_23);

        pushButton_4 = new QPushButton(groupBox_9);
        pushButton_4->setObjectName(QString::fromUtf8("pushButton_4"));
        QSizePolicy sizePolicy7(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy7.setHorizontalStretch(0);
        sizePolicy7.setVerticalStretch(0);
        sizePolicy7.setHeightForWidth(pushButton_4->sizePolicy().hasHeightForWidth());
        pushButton_4->setSizePolicy(sizePolicy7);

        verticalLayout_8->addWidget(pushButton_4);

        pushButton_5 = new QPushButton(groupBox_9);
        pushButton_5->setObjectName(QString::fromUtf8("pushButton_5"));

        verticalLayout_8->addWidget(pushButton_5);


        verticalLayout_15->addWidget(groupBox_9);

        horizontalLayout_20 = new QHBoxLayout();
        horizontalLayout_20->setObjectName(QString::fromUtf8("horizontalLayout_20"));
        groupBox_7 = new QGroupBox(layoutWidget1);
        groupBox_7->setObjectName(QString::fromUtf8("groupBox_7"));
        sizePolicy1.setHeightForWidth(groupBox_7->sizePolicy().hasHeightForWidth());
        groupBox_7->setSizePolicy(sizePolicy1);
        verticalLayout_5 = new QVBoxLayout(groupBox_7);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
        label_10 = new QLabel(groupBox_7);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        horizontalLayout_10->addWidget(label_10);

        lineEdit_8 = new QLineEdit(groupBox_7);
        lineEdit_8->setObjectName(QString::fromUtf8("lineEdit_8"));
        sizePolicy5.setHeightForWidth(lineEdit_8->sizePolicy().hasHeightForWidth());
        lineEdit_8->setSizePolicy(sizePolicy5);

        horizontalLayout_10->addWidget(lineEdit_8);


        verticalLayout_5->addLayout(horizontalLayout_10);

        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setObjectName(QString::fromUtf8("horizontalLayout_11"));
        label_11 = new QLabel(groupBox_7);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        horizontalLayout_11->addWidget(label_11);

        lineEdit_9 = new QLineEdit(groupBox_7);
        lineEdit_9->setObjectName(QString::fromUtf8("lineEdit_9"));
        sizePolicy5.setHeightForWidth(lineEdit_9->sizePolicy().hasHeightForWidth());
        lineEdit_9->setSizePolicy(sizePolicy5);

        horizontalLayout_11->addWidget(lineEdit_9);


        verticalLayout_5->addLayout(horizontalLayout_11);

        horizontalLayout_12 = new QHBoxLayout();
        horizontalLayout_12->setObjectName(QString::fromUtf8("horizontalLayout_12"));
        label_12 = new QLabel(groupBox_7);
        label_12->setObjectName(QString::fromUtf8("label_12"));

        horizontalLayout_12->addWidget(label_12);

        lineEdit_10 = new QLineEdit(groupBox_7);
        lineEdit_10->setObjectName(QString::fromUtf8("lineEdit_10"));
        sizePolicy5.setHeightForWidth(lineEdit_10->sizePolicy().hasHeightForWidth());
        lineEdit_10->setSizePolicy(sizePolicy5);

        horizontalLayout_12->addWidget(lineEdit_10);


        verticalLayout_5->addLayout(horizontalLayout_12);

        horizontalLayout_13 = new QHBoxLayout();
        horizontalLayout_13->setObjectName(QString::fromUtf8("horizontalLayout_13"));
        label_13 = new QLabel(groupBox_7);
        label_13->setObjectName(QString::fromUtf8("label_13"));

        horizontalLayout_13->addWidget(label_13);

        lineEdit_11 = new QLineEdit(groupBox_7);
        lineEdit_11->setObjectName(QString::fromUtf8("lineEdit_11"));
        sizePolicy6.setHeightForWidth(lineEdit_11->sizePolicy().hasHeightForWidth());
        lineEdit_11->setSizePolicy(sizePolicy6);

        horizontalLayout_13->addWidget(lineEdit_11);


        verticalLayout_5->addLayout(horizontalLayout_13);

        checkBox_18 = new QCheckBox(groupBox_7);
        checkBox_18->setObjectName(QString::fromUtf8("checkBox_18"));

        verticalLayout_5->addWidget(checkBox_18);


        horizontalLayout_20->addWidget(groupBox_7);


        verticalLayout_15->addLayout(horizontalLayout_20);


        horizontalLayout_24->addLayout(verticalLayout_15);

        verticalLayout_14 = new QVBoxLayout();
        verticalLayout_14->setObjectName(QString::fromUtf8("verticalLayout_14"));
        groupBox_8 = new QGroupBox(layoutWidget1);
        groupBox_8->setObjectName(QString::fromUtf8("groupBox_8"));
        QSizePolicy sizePolicy8(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        sizePolicy8.setHorizontalStretch(0);
        sizePolicy8.setVerticalStretch(0);
        sizePolicy8.setHeightForWidth(groupBox_8->sizePolicy().hasHeightForWidth());
        groupBox_8->setSizePolicy(sizePolicy8);
        groupBox_8->setMaximumSize(QSize(400, 16777215));
        verticalLayout_12 = new QVBoxLayout(groupBox_8);
        verticalLayout_12->setObjectName(QString::fromUtf8("verticalLayout_12"));
        verticalLayout_7 = new QVBoxLayout();
        verticalLayout_7->setObjectName(QString::fromUtf8("verticalLayout_7"));
        checkBox_19 = new QCheckBox(groupBox_8);
        checkBox_19->setObjectName(QString::fromUtf8("checkBox_19"));

        verticalLayout_7->addWidget(checkBox_19);

        checkBox_20 = new QCheckBox(groupBox_8);
        checkBox_20->setObjectName(QString::fromUtf8("checkBox_20"));

        verticalLayout_7->addWidget(checkBox_20);

        checkBox_21 = new QCheckBox(groupBox_8);
        checkBox_21->setObjectName(QString::fromUtf8("checkBox_21"));

        verticalLayout_7->addWidget(checkBox_21);

        pushButton_3 = new QPushButton(groupBox_8);
        pushButton_3->setObjectName(QString::fromUtf8("pushButton_3"));

        verticalLayout_7->addWidget(pushButton_3);


        verticalLayout_12->addLayout(verticalLayout_7);

        groupBox_2 = new QGroupBox(groupBox_8);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        QSizePolicy sizePolicy9(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy9.setHorizontalStretch(0);
        sizePolicy9.setVerticalStretch(0);
        sizePolicy9.setHeightForWidth(groupBox_2->sizePolicy().hasHeightForWidth());
        groupBox_2->setSizePolicy(sizePolicy9);
        verticalLayout_11 = new QVBoxLayout(groupBox_2);
        verticalLayout_11->setObjectName(QString::fromUtf8("verticalLayout_11"));
        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        checkBox_11 = new QCheckBox(groupBox_2);
        checkBox_11->setObjectName(QString::fromUtf8("checkBox_11"));
        sizePolicy9.setHeightForWidth(checkBox_11->sizePolicy().hasHeightForWidth());
        checkBox_11->setSizePolicy(sizePolicy9);

        horizontalLayout_6->addWidget(checkBox_11);

        checkBox_13 = new QCheckBox(groupBox_2);
        checkBox_13->setObjectName(QString::fromUtf8("checkBox_13"));
        sizePolicy9.setHeightForWidth(checkBox_13->sizePolicy().hasHeightForWidth());
        checkBox_13->setSizePolicy(sizePolicy9);

        horizontalLayout_6->addWidget(checkBox_13);

        checkBox_14 = new QCheckBox(groupBox_2);
        checkBox_14->setObjectName(QString::fromUtf8("checkBox_14"));
        sizePolicy9.setHeightForWidth(checkBox_14->sizePolicy().hasHeightForWidth());
        checkBox_14->setSizePolicy(sizePolicy9);

        horizontalLayout_6->addWidget(checkBox_14);


        verticalLayout_11->addLayout(horizontalLayout_6);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        checkBox_12 = new QCheckBox(groupBox_2);
        checkBox_12->setObjectName(QString::fromUtf8("checkBox_12"));
        sizePolicy9.setHeightForWidth(checkBox_12->sizePolicy().hasHeightForWidth());
        checkBox_12->setSizePolicy(sizePolicy9);

        horizontalLayout_7->addWidget(checkBox_12);

        checkBox_15 = new QCheckBox(groupBox_2);
        checkBox_15->setObjectName(QString::fromUtf8("checkBox_15"));
        sizePolicy9.setHeightForWidth(checkBox_15->sizePolicy().hasHeightForWidth());
        checkBox_15->setSizePolicy(sizePolicy9);

        horizontalLayout_7->addWidget(checkBox_15);

        checkBox_16 = new QCheckBox(groupBox_2);
        checkBox_16->setObjectName(QString::fromUtf8("checkBox_16"));
        sizePolicy9.setHeightForWidth(checkBox_16->sizePolicy().hasHeightForWidth());
        checkBox_16->setSizePolicy(sizePolicy9);

        horizontalLayout_7->addWidget(checkBox_16);


        verticalLayout_11->addLayout(horizontalLayout_7);

        checkBox_17 = new QCheckBox(groupBox_2);
        checkBox_17->setObjectName(QString::fromUtf8("checkBox_17"));
        sizePolicy9.setHeightForWidth(checkBox_17->sizePolicy().hasHeightForWidth());
        checkBox_17->setSizePolicy(sizePolicy9);

        verticalLayout_11->addWidget(checkBox_17);

        groupBox_6 = new QGroupBox(groupBox_2);
        groupBox_6->setObjectName(QString::fromUtf8("groupBox_6"));
        QSizePolicy sizePolicy10(QSizePolicy::Preferred, QSizePolicy::Expanding);
        sizePolicy10.setHorizontalStretch(0);
        sizePolicy10.setVerticalStretch(0);
        sizePolicy10.setHeightForWidth(groupBox_6->sizePolicy().hasHeightForWidth());
        groupBox_6->setSizePolicy(sizePolicy10);
        verticalLayout_4 = new QVBoxLayout(groupBox_6);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
        label_6 = new QLabel(groupBox_6);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        QSizePolicy sizePolicy11(QSizePolicy::MinimumExpanding, QSizePolicy::Preferred);
        sizePolicy11.setHorizontalStretch(0);
        sizePolicy11.setVerticalStretch(0);
        sizePolicy11.setHeightForWidth(label_6->sizePolicy().hasHeightForWidth());
        label_6->setSizePolicy(sizePolicy11);
        label_6->setAlignment(Qt::AlignCenter);

        horizontalLayout_9->addWidget(label_6);

        lineEdit_6 = new QLineEdit(groupBox_6);
        lineEdit_6->setObjectName(QString::fromUtf8("lineEdit_6"));
        QSizePolicy sizePolicy12(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
        sizePolicy12.setHorizontalStretch(0);
        sizePolicy12.setVerticalStretch(0);
        sizePolicy12.setHeightForWidth(lineEdit_6->sizePolicy().hasHeightForWidth());
        lineEdit_6->setSizePolicy(sizePolicy12);

        horizontalLayout_9->addWidget(lineEdit_6);

        label_7 = new QLabel(groupBox_6);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        sizePolicy9.setHeightForWidth(label_7->sizePolicy().hasHeightForWidth());
        label_7->setSizePolicy(sizePolicy9);
        label_7->setAlignment(Qt::AlignCenter);

        horizontalLayout_9->addWidget(label_7);


        verticalLayout_4->addLayout(horizontalLayout_9);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        label_8 = new QLabel(groupBox_6);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        sizePolicy11.setHeightForWidth(label_8->sizePolicy().hasHeightForWidth());
        label_8->setSizePolicy(sizePolicy11);
        label_8->setAlignment(Qt::AlignCenter);

        horizontalLayout_8->addWidget(label_8);

        lineEdit_7 = new QLineEdit(groupBox_6);
        lineEdit_7->setObjectName(QString::fromUtf8("lineEdit_7"));
        sizePolicy12.setHeightForWidth(lineEdit_7->sizePolicy().hasHeightForWidth());
        lineEdit_7->setSizePolicy(sizePolicy12);

        horizontalLayout_8->addWidget(lineEdit_7);

        label_9 = new QLabel(groupBox_6);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        sizePolicy11.setHeightForWidth(label_9->sizePolicy().hasHeightForWidth());
        label_9->setSizePolicy(sizePolicy11);
        label_9->setAlignment(Qt::AlignCenter);

        horizontalLayout_8->addWidget(label_9);


        verticalLayout_4->addLayout(horizontalLayout_8);


        verticalLayout_11->addWidget(groupBox_6);


        verticalLayout_12->addWidget(groupBox_2);


        verticalLayout_14->addWidget(groupBox_8);

        groupBox_11 = new QGroupBox(layoutWidget1);
        groupBox_11->setObjectName(QString::fromUtf8("groupBox_11"));
        sizePolicy8.setHeightForWidth(groupBox_11->sizePolicy().hasHeightForWidth());
        groupBox_11->setSizePolicy(sizePolicy8);
        groupBox_11->setMaximumSize(QSize(400, 400));
        verticalLayout_6 = new QVBoxLayout(groupBox_11);
        verticalLayout_6->setObjectName(QString::fromUtf8("verticalLayout_6"));
        debug_window = new QTextEdit(groupBox_11);
        debug_window->setObjectName(QString::fromUtf8("debug_window"));
        debug_window->setEnabled(false);

        verticalLayout_6->addWidget(debug_window);


        verticalLayout_14->addWidget(groupBox_11);


        horizontalLayout_24->addLayout(verticalLayout_14);

        verticalLayout_13 = new QVBoxLayout();
        verticalLayout_13->setObjectName(QString::fromUtf8("verticalLayout_13"));
        groupBox = new QGroupBox(layoutWidget1);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        sizePolicy1.setHeightForWidth(groupBox->sizePolicy().hasHeightForWidth());
        groupBox->setSizePolicy(sizePolicy1);
        verticalLayout_3 = new QVBoxLayout(groupBox);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        checkBox = new QCheckBox(groupBox);
        checkBox->setObjectName(QString::fromUtf8("checkBox"));

        verticalLayout_3->addWidget(checkBox);

        checkBox_2 = new QCheckBox(groupBox);
        checkBox_2->setObjectName(QString::fromUtf8("checkBox_2"));

        verticalLayout_3->addWidget(checkBox_2);

        checkBox_3 = new QCheckBox(groupBox);
        checkBox_3->setObjectName(QString::fromUtf8("checkBox_3"));

        verticalLayout_3->addWidget(checkBox_3);

        checkBox_4 = new QCheckBox(groupBox);
        checkBox_4->setObjectName(QString::fromUtf8("checkBox_4"));

        verticalLayout_3->addWidget(checkBox_4);

        checkBox_6 = new QCheckBox(groupBox);
        checkBox_6->setObjectName(QString::fromUtf8("checkBox_6"));

        verticalLayout_3->addWidget(checkBox_6);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        horizontalLayout_3->addWidget(label_3);

        lineEdit_3 = new QLineEdit(groupBox);
        lineEdit_3->setObjectName(QString::fromUtf8("lineEdit_3"));
        sizePolicy6.setHeightForWidth(lineEdit_3->sizePolicy().hasHeightForWidth());
        lineEdit_3->setSizePolicy(sizePolicy6);

        horizontalLayout_3->addWidget(lineEdit_3);


        verticalLayout_3->addLayout(horizontalLayout_3);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        horizontalLayout_4->addWidget(label_4);

        lineEdit_4 = new QLineEdit(groupBox);
        lineEdit_4->setObjectName(QString::fromUtf8("lineEdit_4"));
        sizePolicy6.setHeightForWidth(lineEdit_4->sizePolicy().hasHeightForWidth());
        lineEdit_4->setSizePolicy(sizePolicy6);

        horizontalLayout_4->addWidget(lineEdit_4);


        verticalLayout_3->addLayout(horizontalLayout_4);

        checkBox_7 = new QCheckBox(groupBox);
        checkBox_7->setObjectName(QString::fromUtf8("checkBox_7"));
        checkBox_7->setChecked(true);

        verticalLayout_3->addWidget(checkBox_7);

        checkBox_8 = new QCheckBox(groupBox);
        checkBox_8->setObjectName(QString::fromUtf8("checkBox_8"));

        verticalLayout_3->addWidget(checkBox_8);

        checkBox_9 = new QCheckBox(groupBox);
        checkBox_9->setObjectName(QString::fromUtf8("checkBox_9"));

        verticalLayout_3->addWidget(checkBox_9);


        verticalLayout_13->addWidget(groupBox);

        verticalLayout_10 = new QVBoxLayout();
        verticalLayout_10->setObjectName(QString::fromUtf8("verticalLayout_10"));
        groupBox_3 = new QGroupBox(layoutWidget1);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        sizePolicy1.setHeightForWidth(groupBox_3->sizePolicy().hasHeightForWidth());
        groupBox_3->setSizePolicy(sizePolicy1);
        groupBox_3->setMaximumSize(QSize(10001, 100));
        horizontalLayout_5 = new QHBoxLayout(groupBox_3);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        label_5 = new QLabel(groupBox_3);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        QSizePolicy sizePolicy13(QSizePolicy::Preferred, QSizePolicy::Minimum);
        sizePolicy13.setHorizontalStretch(0);
        sizePolicy13.setVerticalStretch(0);
        sizePolicy13.setHeightForWidth(label_5->sizePolicy().hasHeightForWidth());
        label_5->setSizePolicy(sizePolicy13);

        horizontalLayout_5->addWidget(label_5);

        pushButton = new QPushButton(groupBox_3);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        sizePolicy1.setHeightForWidth(pushButton->sizePolicy().hasHeightForWidth());
        pushButton->setSizePolicy(sizePolicy1);

        horizontalLayout_5->addWidget(pushButton);


        verticalLayout_10->addWidget(groupBox_3);

        groupBox_5 = new QGroupBox(layoutWidget1);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        sizePolicy1.setHeightForWidth(groupBox_5->sizePolicy().hasHeightForWidth());
        groupBox_5->setSizePolicy(sizePolicy1);
        verticalLayout = new QVBoxLayout(groupBox_5);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        radioButton = new QRadioButton(groupBox_5);
        IntergratorsBG = new QButtonGroup(MainWindow);
        IntergratorsBG->setObjectName(QString::fromUtf8("IntergratorsBG"));
        IntergratorsBG->addButton(radioButton);
        radioButton->setObjectName(QString::fromUtf8("radioButton"));
        radioButton->setChecked(true);

        verticalLayout->addWidget(radioButton);

        radioButton_2 = new QRadioButton(groupBox_5);
        IntergratorsBG->addButton(radioButton_2);
        radioButton_2->setObjectName(QString::fromUtf8("radioButton_2"));

        verticalLayout->addWidget(radioButton_2);

        radioButton_3 = new QRadioButton(groupBox_5);
        IntergratorsBG->addButton(radioButton_3);
        radioButton_3->setObjectName(QString::fromUtf8("radioButton_3"));

        verticalLayout->addWidget(radioButton_3);


        verticalLayout_10->addWidget(groupBox_5);

        Refine_Morse_Set_Layout = new QGroupBox(layoutWidget1);
        Refine_Morse_Set_Layout->setObjectName(QString::fromUtf8("Refine_Morse_Set_Layout"));
        sizePolicy1.setHeightForWidth(Refine_Morse_Set_Layout->sizePolicy().hasHeightForWidth());
        Refine_Morse_Set_Layout->setSizePolicy(sizePolicy1);
        verticalLayout_2 = new QVBoxLayout(Refine_Morse_Set_Layout);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetMinimumSize);
        label = new QLabel(Refine_Morse_Set_Layout);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);

        lineEdit = new QLineEdit(Refine_Morse_Set_Layout);
        lineEdit->setObjectName(QString::fromUtf8("lineEdit"));
        sizePolicy6.setHeightForWidth(lineEdit->sizePolicy().hasHeightForWidth());
        lineEdit->setSizePolicy(sizePolicy6);

        horizontalLayout->addWidget(lineEdit);


        verticalLayout_2->addLayout(horizontalLayout);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        label_2 = new QLabel(Refine_Morse_Set_Layout);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout_2->addWidget(label_2);

        lineEdit_2 = new QLineEdit(Refine_Morse_Set_Layout);
        lineEdit_2->setObjectName(QString::fromUtf8("lineEdit_2"));
        sizePolicy6.setHeightForWidth(lineEdit_2->sizePolicy().hasHeightForWidth());
        lineEdit_2->setSizePolicy(sizePolicy6);

        horizontalLayout_2->addWidget(lineEdit_2);


        verticalLayout_2->addLayout(horizontalLayout_2);

        pushButton_2 = new QPushButton(Refine_Morse_Set_Layout);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));
        sizePolicy6.setHeightForWidth(pushButton_2->sizePolicy().hasHeightForWidth());
        pushButton_2->setSizePolicy(sizePolicy6);

        verticalLayout_2->addWidget(pushButton_2);


        verticalLayout_10->addWidget(Refine_Morse_Set_Layout);


        verticalLayout_13->addLayout(verticalLayout_10);


        horizontalLayout_24->addLayout(verticalLayout_13);

        MainWindow->setCentralWidget(centralwidget);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        groupBox_12->setTitle(QCoreApplication::translate("MainWindow", "MCG", nullptr));
        groupBox_13->setTitle(QCoreApplication::translate("MainWindow", "ECG", nullptr));
        groupBox_10->setTitle(QCoreApplication::translate("MainWindow", "Morse Triangle Observe", nullptr));
        label_17->setText(QCoreApplication::translate("MainWindow", "Morse", nullptr));
        lineEdit_15->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        pushButton_6->setText(QCoreApplication::translate("MainWindow", "Display Connection", nullptr));
        label_18->setText(QCoreApplication::translate("MainWindow", "Triangle", nullptr));
        lineEdit_16->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        pushButton_7->setText(QCoreApplication::translate("MainWindow", "Display Connection", nullptr));
        pushButton_8->setText(QCoreApplication::translate("MainWindow", "Back to SCC", nullptr));
        groupBox_9->setTitle(QCoreApplication::translate("MainWindow", "Connection Region", nullptr));
        checkBox_22->setText(QCoreApplication::translate("MainWindow", "Compute Connection Region", nullptr));
        label_14->setText(QCoreApplication::translate("MainWindow", "Tau", nullptr));
        lineEdit_12->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_15->setText(QCoreApplication::translate("MainWindow", "M1", nullptr));
        lineEdit_13->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_16->setText(QCoreApplication::translate("MainWindow", "M2", nullptr));
        lineEdit_14->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        checkBox_23->setText(QCoreApplication::translate("MainWindow", "Fast refinement (not accurate)", nullptr));
        pushButton_4->setText(QCoreApplication::translate("MainWindow", "Refine MCG edge", nullptr));
        pushButton_5->setText(QCoreApplication::translate("MainWindow", "Refine connection region", nullptr));
        groupBox_7->setTitle(QCoreApplication::translate("MainWindow", "Auto Refine", nullptr));
        label_10->setText(QCoreApplication::translate("MainWindow", "Min Priority", nullptr));
        lineEdit_8->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_11->setText(QCoreApplication::translate("MainWindow", "tau_max", nullptr));
        lineEdit_9->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_12->setText(QCoreApplication::translate("MainWindow", "Iter_max", nullptr));
        lineEdit_10->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_13->setText(QCoreApplication::translate("MainWindow", "Initial", nullptr));
        lineEdit_11->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        checkBox_18->setText(QCoreApplication::translate("MainWindow", "No t_max", nullptr));
        groupBox_8->setTitle(QCoreApplication::translate("MainWindow", "Visualize Sampling", nullptr));
        checkBox_19->setText(QCoreApplication::translate("MainWindow", "Visualize Sample Points", nullptr));
        checkBox_20->setText(QCoreApplication::translate("MainWindow", "Show triangle mapping", nullptr));
        checkBox_21->setText(QCoreApplication::translate("MainWindow", "Show Back ward", nullptr));
        pushButton_3->setText(QCoreApplication::translate("MainWindow", "Change Edge", nullptr));
        groupBox_2->setTitle(QCoreApplication::translate("MainWindow", "Visualization", nullptr));
        checkBox_11->setText(QCoreApplication::translate("MainWindow", "Animate the flow", nullptr));
        checkBox_13->setText(QCoreApplication::translate("MainWindow", "IBFV off", nullptr));
        checkBox_14->setText(QCoreApplication::translate("MainWindow", "Flip normal", nullptr));
        checkBox_12->setText(QCoreApplication::translate("MainWindow", "Disable lighting", nullptr));
        checkBox_15->setText(QCoreApplication::translate("MainWindow", "Grey texture", nullptr));
        checkBox_16->setText(QCoreApplication::translate("MainWindow", "Not Shiny", nullptr));
        checkBox_17->setText(QCoreApplication::translate("MainWindow", "Color map of vector field magnitude", nullptr));
        groupBox_6->setTitle(QCoreApplication::translate("MainWindow", "Even streamline placement", nullptr));
        label_6->setText(QCoreApplication::translate("MainWindow", "Separation", nullptr));
        lineEdit_6->setText(QCoreApplication::translate("MainWindow", "0.03", nullptr));
        label_7->setText(QCoreApplication::translate("MainWindow", "* Object radius", nullptr));
        label_8->setText(QCoreApplication::translate("MainWindow", "Shortest", nullptr));
        lineEdit_7->setText(QCoreApplication::translate("MainWindow", "0.03", nullptr));
        label_9->setText(QCoreApplication::translate("MainWindow", "* Object radius", nullptr));
        groupBox_11->setTitle(QCoreApplication::translate("MainWindow", "Debug Window", nullptr));
        groupBox->setTitle(QCoreApplication::translate("MainWindow", "Analysis", nullptr));
        checkBox->setText(QCoreApplication::translate("MainWindow", "Display fixed points", nullptr));
        checkBox_2->setText(QCoreApplication::translate("MainWindow", "Display separatrices", nullptr));
        checkBox_3->setText(QCoreApplication::translate("MainWindow", "Display periodic orbits", nullptr));
        checkBox_4->setText(QCoreApplication::translate("MainWindow", "Show Morse sets", nullptr));
        checkBox_6->setText(QCoreApplication::translate("MainWindow", "Show real ID", nullptr));
        label_3->setText(QCoreApplication::translate("MainWindow", "use a tau", nullptr));
        lineEdit_3->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_4->setText(QCoreApplication::translate("MainWindow", "Error threshold", nullptr));
        lineEdit_4->setText(QCoreApplication::translate("MainWindow", "1e-15", nullptr));
        checkBox_7->setText(QCoreApplication::translate("MainWindow", "Compute Conley index", nullptr));
        checkBox_8->setText(QCoreApplication::translate("MainWindow", "Construct minimal MCG", nullptr));
        checkBox_9->setText(QCoreApplication::translate("MainWindow", "Remove disconnected", nullptr));
        groupBox_3->setTitle(QCoreApplication::translate("MainWindow", "File Operations", nullptr));
        label_5->setText(QCoreApplication::translate("MainWindow", "Load a .ply file", nullptr));
        pushButton->setText(QCoreApplication::translate("MainWindow", "Browsers...", nullptr));
        groupBox_5->setTitle(QCoreApplication::translate("MainWindow", "Integrators", nullptr));
        radioButton->setText(QCoreApplication::translate("MainWindow", "Euler 1", nullptr));
        radioButton_2->setText(QCoreApplication::translate("MainWindow", "RK2", nullptr));
        radioButton_3->setText(QCoreApplication::translate("MainWindow", "RK4", nullptr));
        Refine_Morse_Set_Layout->setTitle(QCoreApplication::translate("MainWindow", "Refine Morse Set", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "ID", nullptr));
        lineEdit->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "tau", nullptr));
        lineEdit_2->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        pushButton_2->setText(QCoreApplication::translate("MainWindow", "Refine", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
