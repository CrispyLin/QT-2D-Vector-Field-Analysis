#include "MainWindow.h"
#include "ui_MainWindow.h"

#include "Windows/VectorFieldWindow.h"
#include "Windows/MCGWindow.h"
#include "Windows/ECGWindow.h"

// some extern varibles
extern clock_t g_start;
extern clock_t g_finish;

extern Polyhedron *object; // declared in geometry.cpp

extern MorseDecomp *morse_decomp;
extern MorseDecomp *local_decomp;
extern TrajectoryList *separatrices;
extern PeriodicOrbitList *periodic_orbits;
extern EvenStreamlinePlace *evenplace;

extern MorseDecomp *L1_morse;
extern MorseDecomp *L2_morse;
extern ECG_Graph *ecg;
extern MCG_Graph *mcg;

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
    this->ui->VF_Window->set_up_MainWindow_ptr(this);
    qInfo()  << QString("Main Window has been initialized");
}

void MainWindow::print_debug_message(QString message)
{
    this->ui->debug_Window->append(message);
}

// should only be called when you done updating MCG
void MainWindow::set_MCGOn(bool flag)
{
    this->ui->MCG_Window->ShowMCGOn = flag;
    this->ui->MCG_Window->update_scene();
}

// should only be called when you  done updating ECG
void MainWindow::set_ECGOn(bool flag)
{
    this->ui->ECG_Window->ShowECGOn = flag;
    this->ui->ECG_Window->update_scene();
}

// should only be called when you done updating MCG
// set MCG_Window's ShowConleyCircle to the corresponding flag value.
void MainWindow::set_ShowConleyCircle(bool flag)
{
    this->ui->MCG_Window->ShowConleyCircle = flag;
    this->ui->MCG_Window->update_scene();
}


void MainWindow::on_display_fixed_points_clicked(bool checked)
{
    if(!checked){
        this->ui->VF_Window->ShowFixedPtOn = 0; // set ShowFixedPtOn to false/0
        return;
    }


    QString str = "----------------\n";
    g_start = clock();
    str.append("Start detecting fixed points...\n");
    this->ui->VF_Window->detect_FixedPts();
    str.append("Finish detecting fixed points...\n");
    int num_fixedPts = this->ui->VF_Window->get_num_fixedPts();
    str.append(QString("%1 fixed points have been detected\n").arg(num_fixedPts));
    g_finish = clock();
    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;
    str.append( QString( "time for this operation is %1\n").arg( time ) );
    str.append( "----------------" );
    this->print_debug_message(str);
    this->ui->VF_Window->ShowFixedPtOn = 1; // set ShowFixedPtOn to true/1
}




void MainWindow::on_display_separatrices_clicked(bool checked)
{
    if(!checked){
        this->ui->VF_Window->ShowSeparatricesOn = 0; // set Show Separatrices to false;
        return;
    }

    // constructing fixedPts
    if(this->ui->VF_Window->get_num_fixedPts() == 0){
        this->ui->display_fixed_points->setChecked(true);
        this->on_display_fixed_points_clicked(true);
    }

    // constructing separatrices
    {
        QString str = "----------------\n";
        str.append("Start computing the separatrices...\n");

        g_start = clock();
        object->cal_separatrices();
        g_finish = clock();

        int num_sep = separatrices->ntrajs;
        str.append(QString("%1 separatrices have been detected\n").arg(num_sep));

        double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;
        str.append( QString( "time for this operation is %1\n").arg( time ) );
        str.append( "----------------" );
        this->print_debug_message(str);

        this->ui->VF_Window->ShowSeparatricesOn = 1; // set ShowSeparatricesOn to true/1
    }

    // constructing ECGs
    {
        g_start = clock();

        if(ecg != nullptr)
            delete ecg;
        if(periodic_orbits == nullptr)
            ecg = new ECG_Graph(object->slist.nsingularities);
        else
            ecg = new ECG_Graph(object->slist.nsingularities + periodic_orbits->nporbits);

        ecg->init_ECG();
        ecg->build_ecg();


        g_finish = clock();
        QString str = "----------------\n";
        str.append( QString( "%1 nodes and %2 edges of ECG have been computed.\n" ).arg( ecg->nlist->nenodes ).arg( ecg->elist->nedges ) );
        double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;
        str.append( QString("time for constructing the ECG is %1 seconds\n").arg( time ) );
        str.append( "----------------" );
        this->print_debug_message(str);

        this->set_ECGOn(true);
        qInfo()<< "ECG is done";
    }
}

