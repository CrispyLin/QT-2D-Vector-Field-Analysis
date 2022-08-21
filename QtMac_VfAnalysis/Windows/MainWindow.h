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
private slots:
    void on_display_fixed_points_clicked(bool checked);

    void on_display_separatrices_clicked(bool checked);

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
