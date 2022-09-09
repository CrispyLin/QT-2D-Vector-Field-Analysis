#pragma once
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    // public member variables


    // public member functions
    void init();
    void print_debug_message(QString message);

    void set_MCGOn(bool);
    void set_ECGOn(bool);
    void set_ShowConleyCircle(bool);
    void set_VFDisplayOn(bool flag);
private slots:
    void on_display_fixed_points_clicked(bool checked);

    void on_display_separatrices_clicked(bool checked);

    void on_display_periodic_orbits_clicked(bool checked);

    void on_show_morse_sets_clicked(bool checked);

    void on_show_real_ID_clicked(bool checked);

    void on_Show_Conley_MCG_clicked(bool checked);

    void on_ComputeButton_clicked();

    void on_construct_minimal_MCG_clicked(bool checked);

    void on_remove_diconnected_clicked(bool checked);

    void on_Refine_button_clicked();

    void on_Euler_1_clicked(bool checked);

    void on_RK2_clicked(bool checked);

    void on_RK4_clicked(bool checked);

    void on_Flip_normal_clicked(bool checked);

    void on_Grey_texture_clicked(bool checked);

    void on_Disable_lighting_clicked(bool checked);

    void on_Color_map_of_VF_magnitude_clicked(bool checked);

    void on_IBFV_off_clicked(bool checked);

    void on_Place_Streamlines_Button_clicked();

    void on_visualize_sample_points_clicked(bool checked);

    void on_show_triangle_mapping_clicked(bool checked);

    void on_show_backward_clicked(bool checked);

    void on_Change_edge_Button_clicked(bool checked);

    void on_Animate_clicked(bool checked);

    void on_display_connection_region_clicked(bool checked);

    void on_Detect_FixedPts_clicked();

    void on_Compute_sep_clicked();

    void on_Extract_POs_clicked();

    void on_Morse_Decomp_clicked();

    void on_display_streamlines_clicked(bool checked);

    void on_No_T_MAX_clicked(bool checked);

    void on_Auto_Refine_Btn_clicked();

    void on_Browsers_Button_clicked();

    void on_ComputeConnectionRegion_clicked(bool checked);

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
