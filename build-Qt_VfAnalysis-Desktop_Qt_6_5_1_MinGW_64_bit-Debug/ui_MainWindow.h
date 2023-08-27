/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created by: Qt User Interface Compiler version 6.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <MCGWindow.h>
#include <QtCore/QVariant>
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
#include "ECGWindow.h"
#include "VectorFieldWindow.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout_8;
    QVBoxLayout *verticalLayout_12;
    VectorFieldWindow *VF_Window;
    QGroupBox *MCG_GB;
    QHBoxLayout *horizontalLayout_6;
    MCGWindow *MCG_Window;
    QGroupBox *ECG_GB;
    QHBoxLayout *horizontalLayout_7;
    ECGWindow *ECG_Window;
    QVBoxLayout *first_column_layout;
    QGroupBox *Analysis_GB;
    QVBoxLayout *verticalLayout_3;
    QCheckBox *display_fixed_points;
    QPushButton *Detect_FixedPts;
    QCheckBox *display_separatrices;
    QPushButton *Compute_sep;
    QCheckBox *display_periodic_orbits;
    QPushButton *Extract_POs;
    QCheckBox *show_morse_sets;
    QPushButton *Morse_Decomp;
    QCheckBox *show_real_ID;
    QCheckBox *Show_Conley_MCG;
    QHBoxLayout *horizontalLayout_3;
    QLabel *use_a_tau_label;
    QLineEdit *use_a_tau_LineEdit;
    QHBoxLayout *horizontalLayout_4;
    QLabel *error_threshold;
    QLineEdit *error_threshold_LineEdit;
    QPushButton *ComputeButton;
    QCheckBox *construct_minimal_MCG;
    QCheckBox *remove_diconnected;
    QVBoxLayout *verticalLayout_10;
    QGroupBox *File_Operations_GB;
    QHBoxLayout *horizontalLayout_5;
    QLabel *load_a_ply_file_label;
    QPushButton *Browsers_Button;
    QGroupBox *Integrators_GB;
    QVBoxLayout *verticalLayout;
    QRadioButton *Euler_1;
    QRadioButton *RK2;
    QRadioButton *RK4;
    QVBoxLayout *verticalLayout_7;
    QGroupBox *Visualize_Sampling_GB;
    QVBoxLayout *verticalLayout_13;
    QCheckBox *visualize_sample_points;
    QCheckBox *show_triangle_mapping;
    QCheckBox *show_backward;
    QPushButton *Change_edge_Button;
    QGroupBox *Visualization_GB;
    QVBoxLayout *verticalLayout_6;
    QCheckBox *Flip_normal;
    QCheckBox *Grey_texture;
    QCheckBox *Disable_lighting;
    QCheckBox *Animate;
    QCheckBox *Color_map_of_VF_magnitude;
    QCheckBox *IBFV_off;
    QGroupBox *Even_streamline_placement_GB;
    QVBoxLayout *verticalLayout_4;
    QCheckBox *display_streamlines;
    QHBoxLayout *Separation_layout;
    QLabel *separation_label;
    QLineEdit *separation_LineEdit;
    QLabel *label_7;
    QHBoxLayout *Shortest_layout;
    QLabel *shortest_label;
    QLineEdit *shortest_LineEdit;
    QLabel *label_9;
    QPushButton *Place_Streamlines_Button;
    QGroupBox *Debug_GB;
    QVBoxLayout *verticalLayout_11;
    QTextEdit *debug_Window;
    QVBoxLayout *verticalLayout_9;
    QGroupBox *Connection_Region_GB;
    QVBoxLayout *verticalLayout_8;
    QCheckBox *display_connection_region;
    QCheckBox *ComputeConnectionRegion;
    QGroupBox *Refine_Morese_set_GB;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout;
    QLabel *ID_label;
    QLineEdit *local_refine_ID_LineEdit;
    QHBoxLayout *horizontalLayout_2;
    QLabel *tau_label;
    QLineEdit *local_refine_tau_LineEdit;
    QPushButton *Refine_button;
    QGroupBox *Auto_Refine_GB;
    QVBoxLayout *verticalLayout_5;
    QHBoxLayout *min_priority;
    QLabel *label_10;
    QLineEdit *min_priority_LineEdit;
    QHBoxLayout *tau_max;
    QLabel *label_11;
    QLineEdit *tau_max_LineEdit;
    QHBoxLayout *Iter_max;
    QLabel *label_12;
    QLineEdit *Iter_max_LineEdit;
    QHBoxLayout *Initial;
    QLabel *label_13;
    QLineEdit *Inital_LineEdit;
    QCheckBox *No_T_MAX;
    QPushButton *Auto_Refine_Btn;
    QButtonGroup *Integrators_BG;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName("MainWindow");
        MainWindow->resize(1407, 789);
        MainWindow->setMinimumSize(QSize(0, 0));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName("centralwidget");
        centralwidget->setMinimumSize(QSize(0, 0));
        horizontalLayout_8 = new QHBoxLayout(centralwidget);
        horizontalLayout_8->setObjectName("horizontalLayout_8");
        verticalLayout_12 = new QVBoxLayout();
#ifndef Q_OS_MAC
        verticalLayout_12->setSpacing(-1);
#endif
        verticalLayout_12->setObjectName("verticalLayout_12");
        VF_Window = new VectorFieldWindow(centralwidget);
        VF_Window->setObjectName("VF_Window");
        QSizePolicy sizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(80);
        sizePolicy.setVerticalStretch(80);
        sizePolicy.setHeightForWidth(VF_Window->sizePolicy().hasHeightForWidth());
        VF_Window->setSizePolicy(sizePolicy);
        VF_Window->setMinimumSize(QSize(400, 400));
        VF_Window->setMaximumSize(QSize(400, 400));
        VF_Window->setMouseTracking(true);
        VF_Window->setLayoutDirection(Qt::LeftToRight);
        VF_Window->setAutoFillBackground(false);

        verticalLayout_12->addWidget(VF_Window);

        MCG_GB = new QGroupBox(centralwidget);
        MCG_GB->setObjectName("MCG_GB");
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(MCG_GB->sizePolicy().hasHeightForWidth());
        MCG_GB->setSizePolicy(sizePolicy1);
        MCG_GB->setMaximumSize(QSize(410, 300));
        horizontalLayout_6 = new QHBoxLayout(MCG_GB);
        horizontalLayout_6->setObjectName("horizontalLayout_6");
        horizontalLayout_6->setContentsMargins(0, 0, 0, 0);
        MCG_Window = new MCGWindow(MCG_GB);
        MCG_Window->setObjectName("MCG_Window");
        sizePolicy1.setHeightForWidth(MCG_Window->sizePolicy().hasHeightForWidth());
        MCG_Window->setSizePolicy(sizePolicy1);
        MCG_Window->setMinimumSize(QSize(400, 150));
        MCG_Window->setMaximumSize(QSize(400, 150));
        MCG_Window->setMouseTracking(true);

        horizontalLayout_6->addWidget(MCG_Window);


        verticalLayout_12->addWidget(MCG_GB);

        ECG_GB = new QGroupBox(centralwidget);
        ECG_GB->setObjectName("ECG_GB");
        sizePolicy1.setHeightForWidth(ECG_GB->sizePolicy().hasHeightForWidth());
        ECG_GB->setSizePolicy(sizePolicy1);
        ECG_GB->setMaximumSize(QSize(410, 300));
        horizontalLayout_7 = new QHBoxLayout(ECG_GB);
        horizontalLayout_7->setObjectName("horizontalLayout_7");
        horizontalLayout_7->setContentsMargins(0, 0, 0, 0);
        ECG_Window = new ECGWindow(ECG_GB);
        ECG_Window->setObjectName("ECG_Window");
        sizePolicy1.setHeightForWidth(ECG_Window->sizePolicy().hasHeightForWidth());
        ECG_Window->setSizePolicy(sizePolicy1);
        ECG_Window->setMinimumSize(QSize(400, 150));
        ECG_Window->setMaximumSize(QSize(400, 150));

        horizontalLayout_7->addWidget(ECG_Window);


        verticalLayout_12->addWidget(ECG_GB);


        horizontalLayout_8->addLayout(verticalLayout_12);

        first_column_layout = new QVBoxLayout();
        first_column_layout->setSpacing(20);
        first_column_layout->setObjectName("first_column_layout");
        Analysis_GB = new QGroupBox(centralwidget);
        Analysis_GB->setObjectName("Analysis_GB");
        QSizePolicy sizePolicy2(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(Analysis_GB->sizePolicy().hasHeightForWidth());
        Analysis_GB->setSizePolicy(sizePolicy2);
        Analysis_GB->setMaximumSize(QSize(200, 16777215));
        verticalLayout_3 = new QVBoxLayout(Analysis_GB);
        verticalLayout_3->setObjectName("verticalLayout_3");
        verticalLayout_3->setContentsMargins(0, 0, 0, -1);
        display_fixed_points = new QCheckBox(Analysis_GB);
        display_fixed_points->setObjectName("display_fixed_points");

        verticalLayout_3->addWidget(display_fixed_points);

        Detect_FixedPts = new QPushButton(Analysis_GB);
        Detect_FixedPts->setObjectName("Detect_FixedPts");

        verticalLayout_3->addWidget(Detect_FixedPts);

        display_separatrices = new QCheckBox(Analysis_GB);
        display_separatrices->setObjectName("display_separatrices");

        verticalLayout_3->addWidget(display_separatrices);

        Compute_sep = new QPushButton(Analysis_GB);
        Compute_sep->setObjectName("Compute_sep");

        verticalLayout_3->addWidget(Compute_sep);

        display_periodic_orbits = new QCheckBox(Analysis_GB);
        display_periodic_orbits->setObjectName("display_periodic_orbits");

        verticalLayout_3->addWidget(display_periodic_orbits);

        Extract_POs = new QPushButton(Analysis_GB);
        Extract_POs->setObjectName("Extract_POs");

        verticalLayout_3->addWidget(Extract_POs);

        show_morse_sets = new QCheckBox(Analysis_GB);
        show_morse_sets->setObjectName("show_morse_sets");

        verticalLayout_3->addWidget(show_morse_sets);

        Morse_Decomp = new QPushButton(Analysis_GB);
        Morse_Decomp->setObjectName("Morse_Decomp");
        sizePolicy2.setHeightForWidth(Morse_Decomp->sizePolicy().hasHeightForWidth());
        Morse_Decomp->setSizePolicy(sizePolicy2);

        verticalLayout_3->addWidget(Morse_Decomp);

        show_real_ID = new QCheckBox(Analysis_GB);
        show_real_ID->setObjectName("show_real_ID");

        verticalLayout_3->addWidget(show_real_ID);

        Show_Conley_MCG = new QCheckBox(Analysis_GB);
        Show_Conley_MCG->setObjectName("Show_Conley_MCG");

        verticalLayout_3->addWidget(Show_Conley_MCG);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName("horizontalLayout_3");
        use_a_tau_label = new QLabel(Analysis_GB);
        use_a_tau_label->setObjectName("use_a_tau_label");

        horizontalLayout_3->addWidget(use_a_tau_label);

        use_a_tau_LineEdit = new QLineEdit(Analysis_GB);
        use_a_tau_LineEdit->setObjectName("use_a_tau_LineEdit");
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(use_a_tau_LineEdit->sizePolicy().hasHeightForWidth());
        use_a_tau_LineEdit->setSizePolicy(sizePolicy3);

        horizontalLayout_3->addWidget(use_a_tau_LineEdit);


        verticalLayout_3->addLayout(horizontalLayout_3);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName("horizontalLayout_4");
        error_threshold = new QLabel(Analysis_GB);
        error_threshold->setObjectName("error_threshold");

        horizontalLayout_4->addWidget(error_threshold);

        error_threshold_LineEdit = new QLineEdit(Analysis_GB);
        error_threshold_LineEdit->setObjectName("error_threshold_LineEdit");
        sizePolicy3.setHeightForWidth(error_threshold_LineEdit->sizePolicy().hasHeightForWidth());
        error_threshold_LineEdit->setSizePolicy(sizePolicy3);

        horizontalLayout_4->addWidget(error_threshold_LineEdit);


        verticalLayout_3->addLayout(horizontalLayout_4);

        ComputeButton = new QPushButton(Analysis_GB);
        ComputeButton->setObjectName("ComputeButton");
        ComputeButton->setMaximumSize(QSize(1000000, 16777215));

        verticalLayout_3->addWidget(ComputeButton);

        construct_minimal_MCG = new QCheckBox(Analysis_GB);
        construct_minimal_MCG->setObjectName("construct_minimal_MCG");

        verticalLayout_3->addWidget(construct_minimal_MCG);

        remove_diconnected = new QCheckBox(Analysis_GB);
        remove_diconnected->setObjectName("remove_diconnected");

        verticalLayout_3->addWidget(remove_diconnected);


        first_column_layout->addWidget(Analysis_GB);

        verticalLayout_10 = new QVBoxLayout();
        verticalLayout_10->setObjectName("verticalLayout_10");
        File_Operations_GB = new QGroupBox(centralwidget);
        File_Operations_GB->setObjectName("File_Operations_GB");
        sizePolicy2.setHeightForWidth(File_Operations_GB->sizePolicy().hasHeightForWidth());
        File_Operations_GB->setSizePolicy(sizePolicy2);
        File_Operations_GB->setMaximumSize(QSize(200, 100));
        horizontalLayout_5 = new QHBoxLayout(File_Operations_GB);
        horizontalLayout_5->setObjectName("horizontalLayout_5");
        horizontalLayout_5->setContentsMargins(0, 0, 0, 0);
        load_a_ply_file_label = new QLabel(File_Operations_GB);
        load_a_ply_file_label->setObjectName("load_a_ply_file_label");
        QSizePolicy sizePolicy4(QSizePolicy::Preferred, QSizePolicy::Minimum);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(load_a_ply_file_label->sizePolicy().hasHeightForWidth());
        load_a_ply_file_label->setSizePolicy(sizePolicy4);

        horizontalLayout_5->addWidget(load_a_ply_file_label);

        Browsers_Button = new QPushButton(File_Operations_GB);
        Browsers_Button->setObjectName("Browsers_Button");
        sizePolicy2.setHeightForWidth(Browsers_Button->sizePolicy().hasHeightForWidth());
        Browsers_Button->setSizePolicy(sizePolicy2);

        horizontalLayout_5->addWidget(Browsers_Button);


        verticalLayout_10->addWidget(File_Operations_GB);

        Integrators_GB = new QGroupBox(centralwidget);
        Integrators_GB->setObjectName("Integrators_GB");
        sizePolicy2.setHeightForWidth(Integrators_GB->sizePolicy().hasHeightForWidth());
        Integrators_GB->setSizePolicy(sizePolicy2);
        Integrators_GB->setMaximumSize(QSize(200, 16777215));
        verticalLayout = new QVBoxLayout(Integrators_GB);
        verticalLayout->setObjectName("verticalLayout");
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        Euler_1 = new QRadioButton(Integrators_GB);
        Integrators_BG = new QButtonGroup(MainWindow);
        Integrators_BG->setObjectName("Integrators_BG");
        Integrators_BG->addButton(Euler_1);
        Euler_1->setObjectName("Euler_1");
        Euler_1->setChecked(true);

        verticalLayout->addWidget(Euler_1);

        RK2 = new QRadioButton(Integrators_GB);
        Integrators_BG->addButton(RK2);
        RK2->setObjectName("RK2");

        verticalLayout->addWidget(RK2);

        RK4 = new QRadioButton(Integrators_GB);
        Integrators_BG->addButton(RK4);
        RK4->setObjectName("RK4");

        verticalLayout->addWidget(RK4);


        verticalLayout_10->addWidget(Integrators_GB);


        first_column_layout->addLayout(verticalLayout_10);


        horizontalLayout_8->addLayout(first_column_layout);

        verticalLayout_7 = new QVBoxLayout();
        verticalLayout_7->setObjectName("verticalLayout_7");
        Visualize_Sampling_GB = new QGroupBox(centralwidget);
        Visualize_Sampling_GB->setObjectName("Visualize_Sampling_GB");
        sizePolicy1.setHeightForWidth(Visualize_Sampling_GB->sizePolicy().hasHeightForWidth());
        Visualize_Sampling_GB->setSizePolicy(sizePolicy1);
        Visualize_Sampling_GB->setMaximumSize(QSize(300, 16777215));
        verticalLayout_13 = new QVBoxLayout(Visualize_Sampling_GB);
        verticalLayout_13->setObjectName("verticalLayout_13");
        visualize_sample_points = new QCheckBox(Visualize_Sampling_GB);
        visualize_sample_points->setObjectName("visualize_sample_points");
        sizePolicy1.setHeightForWidth(visualize_sample_points->sizePolicy().hasHeightForWidth());
        visualize_sample_points->setSizePolicy(sizePolicy1);

        verticalLayout_13->addWidget(visualize_sample_points);

        show_triangle_mapping = new QCheckBox(Visualize_Sampling_GB);
        show_triangle_mapping->setObjectName("show_triangle_mapping");
        sizePolicy1.setHeightForWidth(show_triangle_mapping->sizePolicy().hasHeightForWidth());
        show_triangle_mapping->setSizePolicy(sizePolicy1);

        verticalLayout_13->addWidget(show_triangle_mapping);

        show_backward = new QCheckBox(Visualize_Sampling_GB);
        show_backward->setObjectName("show_backward");
        sizePolicy1.setHeightForWidth(show_backward->sizePolicy().hasHeightForWidth());
        show_backward->setSizePolicy(sizePolicy1);

        verticalLayout_13->addWidget(show_backward);

        Change_edge_Button = new QPushButton(Visualize_Sampling_GB);
        Change_edge_Button->setObjectName("Change_edge_Button");
        sizePolicy1.setHeightForWidth(Change_edge_Button->sizePolicy().hasHeightForWidth());
        Change_edge_Button->setSizePolicy(sizePolicy1);

        verticalLayout_13->addWidget(Change_edge_Button);


        verticalLayout_7->addWidget(Visualize_Sampling_GB);

        Visualization_GB = new QGroupBox(centralwidget);
        Visualization_GB->setObjectName("Visualization_GB");
        QSizePolicy sizePolicy5(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(Visualization_GB->sizePolicy().hasHeightForWidth());
        Visualization_GB->setSizePolicy(sizePolicy5);
        Visualization_GB->setMaximumSize(QSize(300, 16777215));
        verticalLayout_6 = new QVBoxLayout(Visualization_GB);
        verticalLayout_6->setSpacing(3);
        verticalLayout_6->setObjectName("verticalLayout_6");
        verticalLayout_6->setContentsMargins(0, 0, 0, 0);
        Flip_normal = new QCheckBox(Visualization_GB);
        Flip_normal->setObjectName("Flip_normal");
        sizePolicy1.setHeightForWidth(Flip_normal->sizePolicy().hasHeightForWidth());
        Flip_normal->setSizePolicy(sizePolicy1);

        verticalLayout_6->addWidget(Flip_normal);

        Grey_texture = new QCheckBox(Visualization_GB);
        Grey_texture->setObjectName("Grey_texture");
        sizePolicy1.setHeightForWidth(Grey_texture->sizePolicy().hasHeightForWidth());
        Grey_texture->setSizePolicy(sizePolicy1);

        verticalLayout_6->addWidget(Grey_texture);

        Disable_lighting = new QCheckBox(Visualization_GB);
        Disable_lighting->setObjectName("Disable_lighting");
        QSizePolicy sizePolicy6(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(Disable_lighting->sizePolicy().hasHeightForWidth());
        Disable_lighting->setSizePolicy(sizePolicy6);

        verticalLayout_6->addWidget(Disable_lighting);

        Animate = new QCheckBox(Visualization_GB);
        Animate->setObjectName("Animate");

        verticalLayout_6->addWidget(Animate);

        Color_map_of_VF_magnitude = new QCheckBox(Visualization_GB);
        Color_map_of_VF_magnitude->setObjectName("Color_map_of_VF_magnitude");
        sizePolicy1.setHeightForWidth(Color_map_of_VF_magnitude->sizePolicy().hasHeightForWidth());
        Color_map_of_VF_magnitude->setSizePolicy(sizePolicy1);

        verticalLayout_6->addWidget(Color_map_of_VF_magnitude);

        IBFV_off = new QCheckBox(Visualization_GB);
        IBFV_off->setObjectName("IBFV_off");
        sizePolicy4.setHeightForWidth(IBFV_off->sizePolicy().hasHeightForWidth());
        IBFV_off->setSizePolicy(sizePolicy4);

        verticalLayout_6->addWidget(IBFV_off);

        Even_streamline_placement_GB = new QGroupBox(Visualization_GB);
        Even_streamline_placement_GB->setObjectName("Even_streamline_placement_GB");
        sizePolicy1.setHeightForWidth(Even_streamline_placement_GB->sizePolicy().hasHeightForWidth());
        Even_streamline_placement_GB->setSizePolicy(sizePolicy1);
        verticalLayout_4 = new QVBoxLayout(Even_streamline_placement_GB);
        verticalLayout_4->setObjectName("verticalLayout_4");
        verticalLayout_4->setContentsMargins(0, 0, 0, 0);
        display_streamlines = new QCheckBox(Even_streamline_placement_GB);
        display_streamlines->setObjectName("display_streamlines");

        verticalLayout_4->addWidget(display_streamlines);

        Separation_layout = new QHBoxLayout();
        Separation_layout->setObjectName("Separation_layout");
        separation_label = new QLabel(Even_streamline_placement_GB);
        separation_label->setObjectName("separation_label");
        sizePolicy1.setHeightForWidth(separation_label->sizePolicy().hasHeightForWidth());
        separation_label->setSizePolicy(sizePolicy1);
        separation_label->setAlignment(Qt::AlignCenter);

        Separation_layout->addWidget(separation_label);

        separation_LineEdit = new QLineEdit(Even_streamline_placement_GB);
        separation_LineEdit->setObjectName("separation_LineEdit");
        sizePolicy1.setHeightForWidth(separation_LineEdit->sizePolicy().hasHeightForWidth());
        separation_LineEdit->setSizePolicy(sizePolicy1);
        separation_LineEdit->setMaximumSize(QSize(40, 1000));
        separation_LineEdit->setAlignment(Qt::AlignCenter);

        Separation_layout->addWidget(separation_LineEdit);

        label_7 = new QLabel(Even_streamline_placement_GB);
        label_7->setObjectName("label_7");
        sizePolicy1.setHeightForWidth(label_7->sizePolicy().hasHeightForWidth());
        label_7->setSizePolicy(sizePolicy1);
        label_7->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        Separation_layout->addWidget(label_7);


        verticalLayout_4->addLayout(Separation_layout);

        Shortest_layout = new QHBoxLayout();
        Shortest_layout->setObjectName("Shortest_layout");
        shortest_label = new QLabel(Even_streamline_placement_GB);
        shortest_label->setObjectName("shortest_label");
        sizePolicy1.setHeightForWidth(shortest_label->sizePolicy().hasHeightForWidth());
        shortest_label->setSizePolicy(sizePolicy1);
        shortest_label->setMaximumSize(QSize(100, 16777215));
        shortest_label->setAlignment(Qt::AlignCenter);

        Shortest_layout->addWidget(shortest_label);

        shortest_LineEdit = new QLineEdit(Even_streamline_placement_GB);
        shortest_LineEdit->setObjectName("shortest_LineEdit");
        sizePolicy1.setHeightForWidth(shortest_LineEdit->sizePolicy().hasHeightForWidth());
        shortest_LineEdit->setSizePolicy(sizePolicy1);
        shortest_LineEdit->setMaximumSize(QSize(40, 16777215));
        shortest_LineEdit->setAlignment(Qt::AlignCenter);

        Shortest_layout->addWidget(shortest_LineEdit);

        label_9 = new QLabel(Even_streamline_placement_GB);
        label_9->setObjectName("label_9");
        QSizePolicy sizePolicy7(QSizePolicy::Preferred, QSizePolicy::Maximum);
        sizePolicy7.setHorizontalStretch(0);
        sizePolicy7.setVerticalStretch(0);
        sizePolicy7.setHeightForWidth(label_9->sizePolicy().hasHeightForWidth());
        label_9->setSizePolicy(sizePolicy7);
        label_9->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        Shortest_layout->addWidget(label_9);


        verticalLayout_4->addLayout(Shortest_layout);

        Place_Streamlines_Button = new QPushButton(Even_streamline_placement_GB);
        Place_Streamlines_Button->setObjectName("Place_Streamlines_Button");

        verticalLayout_4->addWidget(Place_Streamlines_Button);


        verticalLayout_6->addWidget(Even_streamline_placement_GB);


        verticalLayout_7->addWidget(Visualization_GB);

        Debug_GB = new QGroupBox(centralwidget);
        Debug_GB->setObjectName("Debug_GB");
        sizePolicy1.setHeightForWidth(Debug_GB->sizePolicy().hasHeightForWidth());
        Debug_GB->setSizePolicy(sizePolicy1);
        Debug_GB->setMaximumSize(QSize(300, 10000));
        verticalLayout_11 = new QVBoxLayout(Debug_GB);
        verticalLayout_11->setObjectName("verticalLayout_11");
        verticalLayout_11->setContentsMargins(0, 0, 0, 0);
        debug_Window = new QTextEdit(Debug_GB);
        debug_Window->setObjectName("debug_Window");
        debug_Window->setEnabled(true);
        sizePolicy1.setHeightForWidth(debug_Window->sizePolicy().hasHeightForWidth());
        debug_Window->setSizePolicy(sizePolicy1);
        debug_Window->setMaximumSize(QSize(10000, 16777215));
        debug_Window->setMouseTracking(false);
        debug_Window->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
        debug_Window->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

        verticalLayout_11->addWidget(debug_Window);


        verticalLayout_7->addWidget(Debug_GB);


        horizontalLayout_8->addLayout(verticalLayout_7);

        verticalLayout_9 = new QVBoxLayout();
        verticalLayout_9->setSpacing(5);
        verticalLayout_9->setObjectName("verticalLayout_9");
        Connection_Region_GB = new QGroupBox(centralwidget);
        Connection_Region_GB->setObjectName("Connection_Region_GB");
        sizePolicy5.setHeightForWidth(Connection_Region_GB->sizePolicy().hasHeightForWidth());
        Connection_Region_GB->setSizePolicy(sizePolicy5);
        Connection_Region_GB->setMaximumSize(QSize(280, 100));
        verticalLayout_8 = new QVBoxLayout(Connection_Region_GB);
        verticalLayout_8->setSpacing(20);
        verticalLayout_8->setObjectName("verticalLayout_8");
        verticalLayout_8->setSizeConstraint(QLayout::SetMinimumSize);
        verticalLayout_8->setContentsMargins(0, 0, 0, 0);
        display_connection_region = new QCheckBox(Connection_Region_GB);
        display_connection_region->setObjectName("display_connection_region");

        verticalLayout_8->addWidget(display_connection_region);

        ComputeConnectionRegion = new QCheckBox(Connection_Region_GB);
        ComputeConnectionRegion->setObjectName("ComputeConnectionRegion");

        verticalLayout_8->addWidget(ComputeConnectionRegion);


        verticalLayout_9->addWidget(Connection_Region_GB);

        Refine_Morese_set_GB = new QGroupBox(centralwidget);
        Refine_Morese_set_GB->setObjectName("Refine_Morese_set_GB");
        sizePolicy2.setHeightForWidth(Refine_Morese_set_GB->sizePolicy().hasHeightForWidth());
        Refine_Morese_set_GB->setSizePolicy(sizePolicy2);
        Refine_Morese_set_GB->setMaximumSize(QSize(280, 16777215));
        verticalLayout_2 = new QVBoxLayout(Refine_Morese_set_GB);
        verticalLayout_2->setObjectName("verticalLayout_2");
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName("horizontalLayout");
        horizontalLayout->setSizeConstraint(QLayout::SetMinimumSize);
        ID_label = new QLabel(Refine_Morese_set_GB);
        ID_label->setObjectName("ID_label");

        horizontalLayout->addWidget(ID_label);

        local_refine_ID_LineEdit = new QLineEdit(Refine_Morese_set_GB);
        local_refine_ID_LineEdit->setObjectName("local_refine_ID_LineEdit");
        sizePolicy3.setHeightForWidth(local_refine_ID_LineEdit->sizePolicy().hasHeightForWidth());
        local_refine_ID_LineEdit->setSizePolicy(sizePolicy3);

        horizontalLayout->addWidget(local_refine_ID_LineEdit);


        verticalLayout_2->addLayout(horizontalLayout);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName("horizontalLayout_2");
        tau_label = new QLabel(Refine_Morese_set_GB);
        tau_label->setObjectName("tau_label");

        horizontalLayout_2->addWidget(tau_label);

        local_refine_tau_LineEdit = new QLineEdit(Refine_Morese_set_GB);
        local_refine_tau_LineEdit->setObjectName("local_refine_tau_LineEdit");
        sizePolicy3.setHeightForWidth(local_refine_tau_LineEdit->sizePolicy().hasHeightForWidth());
        local_refine_tau_LineEdit->setSizePolicy(sizePolicy3);

        horizontalLayout_2->addWidget(local_refine_tau_LineEdit);


        verticalLayout_2->addLayout(horizontalLayout_2);

        Refine_button = new QPushButton(Refine_Morese_set_GB);
        Refine_button->setObjectName("Refine_button");
        sizePolicy3.setHeightForWidth(Refine_button->sizePolicy().hasHeightForWidth());
        Refine_button->setSizePolicy(sizePolicy3);

        verticalLayout_2->addWidget(Refine_button);


        verticalLayout_9->addWidget(Refine_Morese_set_GB);

        Auto_Refine_GB = new QGroupBox(centralwidget);
        Auto_Refine_GB->setObjectName("Auto_Refine_GB");
        sizePolicy2.setHeightForWidth(Auto_Refine_GB->sizePolicy().hasHeightForWidth());
        Auto_Refine_GB->setSizePolicy(sizePolicy2);
        Auto_Refine_GB->setMaximumSize(QSize(280, 16777215));
        verticalLayout_5 = new QVBoxLayout(Auto_Refine_GB);
        verticalLayout_5->setObjectName("verticalLayout_5");
        verticalLayout_5->setContentsMargins(0, 0, 0, 0);
        min_priority = new QHBoxLayout();
        min_priority->setObjectName("min_priority");
        label_10 = new QLabel(Auto_Refine_GB);
        label_10->setObjectName("label_10");

        min_priority->addWidget(label_10);

        min_priority_LineEdit = new QLineEdit(Auto_Refine_GB);
        min_priority_LineEdit->setObjectName("min_priority_LineEdit");
        QSizePolicy sizePolicy8(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy8.setHorizontalStretch(0);
        sizePolicy8.setVerticalStretch(0);
        sizePolicy8.setHeightForWidth(min_priority_LineEdit->sizePolicy().hasHeightForWidth());
        min_priority_LineEdit->setSizePolicy(sizePolicy8);

        min_priority->addWidget(min_priority_LineEdit);


        verticalLayout_5->addLayout(min_priority);

        tau_max = new QHBoxLayout();
        tau_max->setObjectName("tau_max");
        label_11 = new QLabel(Auto_Refine_GB);
        label_11->setObjectName("label_11");

        tau_max->addWidget(label_11);

        tau_max_LineEdit = new QLineEdit(Auto_Refine_GB);
        tau_max_LineEdit->setObjectName("tau_max_LineEdit");
        sizePolicy8.setHeightForWidth(tau_max_LineEdit->sizePolicy().hasHeightForWidth());
        tau_max_LineEdit->setSizePolicy(sizePolicy8);

        tau_max->addWidget(tau_max_LineEdit);


        verticalLayout_5->addLayout(tau_max);

        Iter_max = new QHBoxLayout();
        Iter_max->setObjectName("Iter_max");
        label_12 = new QLabel(Auto_Refine_GB);
        label_12->setObjectName("label_12");

        Iter_max->addWidget(label_12);

        Iter_max_LineEdit = new QLineEdit(Auto_Refine_GB);
        Iter_max_LineEdit->setObjectName("Iter_max_LineEdit");
        sizePolicy8.setHeightForWidth(Iter_max_LineEdit->sizePolicy().hasHeightForWidth());
        Iter_max_LineEdit->setSizePolicy(sizePolicy8);

        Iter_max->addWidget(Iter_max_LineEdit);


        verticalLayout_5->addLayout(Iter_max);

        Initial = new QHBoxLayout();
        Initial->setObjectName("Initial");
        label_13 = new QLabel(Auto_Refine_GB);
        label_13->setObjectName("label_13");

        Initial->addWidget(label_13);

        Inital_LineEdit = new QLineEdit(Auto_Refine_GB);
        Inital_LineEdit->setObjectName("Inital_LineEdit");
        sizePolicy3.setHeightForWidth(Inital_LineEdit->sizePolicy().hasHeightForWidth());
        Inital_LineEdit->setSizePolicy(sizePolicy3);

        Initial->addWidget(Inital_LineEdit);


        verticalLayout_5->addLayout(Initial);

        No_T_MAX = new QCheckBox(Auto_Refine_GB);
        No_T_MAX->setObjectName("No_T_MAX");

        verticalLayout_5->addWidget(No_T_MAX);

        Auto_Refine_Btn = new QPushButton(Auto_Refine_GB);
        Auto_Refine_Btn->setObjectName("Auto_Refine_Btn");

        verticalLayout_5->addWidget(Auto_Refine_Btn);


        verticalLayout_9->addWidget(Auto_Refine_GB);


        horizontalLayout_8->addLayout(verticalLayout_9);

        MainWindow->setCentralWidget(centralwidget);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        MCG_GB->setTitle(QCoreApplication::translate("MainWindow", "MCG", nullptr));
        ECG_GB->setTitle(QCoreApplication::translate("MainWindow", "ECG", nullptr));
        Analysis_GB->setTitle(QCoreApplication::translate("MainWindow", "Analysis", nullptr));
        display_fixed_points->setText(QCoreApplication::translate("MainWindow", "Display fixed points", nullptr));
        Detect_FixedPts->setText(QCoreApplication::translate("MainWindow", "Detect fixed points", nullptr));
        display_separatrices->setText(QCoreApplication::translate("MainWindow", "Display separatrices", nullptr));
        Compute_sep->setText(QCoreApplication::translate("MainWindow", "Compute separatrices", nullptr));
        display_periodic_orbits->setText(QCoreApplication::translate("MainWindow", "Display periodic orbits", nullptr));
        Extract_POs->setText(QCoreApplication::translate("MainWindow", "Extract periodic orbits", nullptr));
        show_morse_sets->setText(QCoreApplication::translate("MainWindow", "Show Morse sets", nullptr));
        Morse_Decomp->setText(QCoreApplication::translate("MainWindow", "Morse Decomp", nullptr));
        show_real_ID->setText(QCoreApplication::translate("MainWindow", "Show real ID", nullptr));
        Show_Conley_MCG->setText(QCoreApplication::translate("MainWindow", "Show Conley MCG", nullptr));
        use_a_tau_label->setText(QCoreApplication::translate("MainWindow", "use a tau", nullptr));
        use_a_tau_LineEdit->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        error_threshold->setText(QCoreApplication::translate("MainWindow", "Error threshold", nullptr));
        error_threshold_LineEdit->setText(QCoreApplication::translate("MainWindow", "1e-15", nullptr));
        ComputeButton->setText(QCoreApplication::translate("MainWindow", "Compute", nullptr));
        construct_minimal_MCG->setText(QCoreApplication::translate("MainWindow", "Construct minimal MCG", nullptr));
        remove_diconnected->setText(QCoreApplication::translate("MainWindow", "Remove disconnected", nullptr));
        File_Operations_GB->setTitle(QCoreApplication::translate("MainWindow", "File Operations", nullptr));
        load_a_ply_file_label->setText(QCoreApplication::translate("MainWindow", "Load a .ply file", nullptr));
        Browsers_Button->setText(QCoreApplication::translate("MainWindow", "Browsers...", nullptr));
        Integrators_GB->setTitle(QCoreApplication::translate("MainWindow", "Integrators", nullptr));
        Euler_1->setText(QCoreApplication::translate("MainWindow", "Euler 1", nullptr));
        RK2->setText(QCoreApplication::translate("MainWindow", "RK2", nullptr));
        RK4->setText(QCoreApplication::translate("MainWindow", "RK4", nullptr));
        Visualize_Sampling_GB->setTitle(QCoreApplication::translate("MainWindow", "Visualize Sampling", nullptr));
        visualize_sample_points->setText(QCoreApplication::translate("MainWindow", "Visualize Sample Points", nullptr));
        show_triangle_mapping->setText(QCoreApplication::translate("MainWindow", "Show triangle mapping", nullptr));
        show_backward->setText(QCoreApplication::translate("MainWindow", "Show Backward", nullptr));
        Change_edge_Button->setText(QCoreApplication::translate("MainWindow", "Change Edge", nullptr));
        Visualization_GB->setTitle(QCoreApplication::translate("MainWindow", "Visualization", nullptr));
        Flip_normal->setText(QCoreApplication::translate("MainWindow", "Flip normal", nullptr));
        Grey_texture->setText(QCoreApplication::translate("MainWindow", "Grey texture", nullptr));
        Disable_lighting->setText(QCoreApplication::translate("MainWindow", "Disable lighting", nullptr));
        Animate->setText(QCoreApplication::translate("MainWindow", "Animate the flow", nullptr));
        Color_map_of_VF_magnitude->setText(QCoreApplication::translate("MainWindow", "Color map of vector field magnitude", nullptr));
        IBFV_off->setText(QCoreApplication::translate("MainWindow", "IBFV off", nullptr));
        Even_streamline_placement_GB->setTitle(QCoreApplication::translate("MainWindow", "Even streamline placement", nullptr));
        display_streamlines->setText(QCoreApplication::translate("MainWindow", "Display Streamlines", nullptr));
        separation_label->setText(QCoreApplication::translate("MainWindow", "Separation", nullptr));
        separation_LineEdit->setText(QCoreApplication::translate("MainWindow", "0.03", nullptr));
        label_7->setText(QCoreApplication::translate("MainWindow", "* Object radius", nullptr));
        shortest_label->setText(QCoreApplication::translate("MainWindow", "Length", nullptr));
        shortest_LineEdit->setText(QCoreApplication::translate("MainWindow", "0.03", nullptr));
        label_9->setText(QCoreApplication::translate("MainWindow", "* Object radius", nullptr));
        Place_Streamlines_Button->setText(QCoreApplication::translate("MainWindow", "Place Streamlines", nullptr));
        Debug_GB->setTitle(QCoreApplication::translate("MainWindow", "Debug Window", nullptr));
        Connection_Region_GB->setTitle(QCoreApplication::translate("MainWindow", "Connection Region", nullptr));
        display_connection_region->setText(QCoreApplication::translate("MainWindow", "Display Connection Region", nullptr));
        ComputeConnectionRegion->setText(QCoreApplication::translate("MainWindow", "Compute Connection Region", nullptr));
        Refine_Morese_set_GB->setTitle(QCoreApplication::translate("MainWindow", "Refine Morse Set", nullptr));
        ID_label->setText(QCoreApplication::translate("MainWindow", "ID", nullptr));
        local_refine_ID_LineEdit->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        tau_label->setText(QCoreApplication::translate("MainWindow", "tau", nullptr));
        local_refine_tau_LineEdit->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        Refine_button->setText(QCoreApplication::translate("MainWindow", "Refine", nullptr));
        Auto_Refine_GB->setTitle(QCoreApplication::translate("MainWindow", "Auto Refine", nullptr));
        label_10->setText(QCoreApplication::translate("MainWindow", "Min Priority", nullptr));
        min_priority_LineEdit->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_11->setText(QCoreApplication::translate("MainWindow", "tau_max", nullptr));
        tau_max_LineEdit->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_12->setText(QCoreApplication::translate("MainWindow", "Iter_max", nullptr));
        Iter_max_LineEdit->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_13->setText(QCoreApplication::translate("MainWindow", "Initial", nullptr));
        Inital_LineEdit->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        No_T_MAX->setText(QCoreApplication::translate("MainWindow", "No t_max", nullptr));
        Auto_Refine_Btn->setText(QCoreApplication::translate("MainWindow", "Auto Refine", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
