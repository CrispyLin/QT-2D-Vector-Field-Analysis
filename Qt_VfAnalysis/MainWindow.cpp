#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "VectorFieldWindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    this->init();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::init(){
    this->ui->openGLWindow->set_up_MainWindow_ptr(this);
    qInfo()  << QString("Main Window has been initialized");
}

void MainWindow::print_debug_message(QString message)
{
    this->ui->debug_window->append(message);
}
