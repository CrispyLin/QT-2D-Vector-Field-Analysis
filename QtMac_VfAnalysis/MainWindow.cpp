#include <qfiledialog.h>
#include "MainWindow.h"
#include "ui_MainWindow.h"

#include "VectorFieldWindow.h"
#include "MCGWindow.h"
#include "ECGWindow.h"

// some extern varibles
extern QString dataset_folder;

extern clock_t g_start;
extern clock_t g_finish;

extern Polyhedron *object; // declared in geometry.cpp

extern MorseDecomp *morse_decomp;
extern MorseDecomp *local_decomp;
extern TrajectoryList *separatrices;
extern PeriodicOrbitList *periodic_orbits;
extern EvenStreamlinePlace *evenplace;
extern bool showRealIndex;
extern double global_tau;
extern double edge_sample_error;
extern bool RemRedundantMCGEdges;
extern bool RemoveDisconnMSOn;
extern RegionTauMap* regiontau;
extern int Integrator_opt;
extern bool NotShiny;
extern bool EnGreyTexture;
extern bool DisableLighting;
extern int ndisplay_trajs;
extern bool Cal_Regions;
bool noTmax = false;

// for displaying even streamlines
extern double StreamlinePlace_dsep;
extern double StreamlinePlace_slength;
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
}

// should only be called when you  done updating ECG
void MainWindow::set_ECGOn(bool flag)
{
    this->ui->ECG_Window->ShowECGOn = flag;
}

void MainWindow::set_VFDisplayOn(bool flag){
    this->ui->VF_Window->display = flag;
}

// should only be called when you done updating MCG
// set MCG_Window's ShowConleyCircle to the corresponding flag value.
void MainWindow::set_ShowConleyCircle(bool checked)
{
    bool flag = this->ui->Show_Conley_MCG->isChecked();
    this->ui->MCG_Window->ShowConleyCircle = flag;

}


void MainWindow::on_display_fixed_points_clicked(bool checked)
{
    if(!checked){
        this->ui->VF_Window->ShowFixedPtOn = 0; // set ShowFixedPtOn to false/0
        return;
    }

    this->ui->VF_Window->ShowFixedPtOn = 1; // set ShowFixedPtOn to true/1
}




void MainWindow::on_display_separatrices_clicked(bool checked)
{
    if(!checked){
        this->ui->VF_Window->ShowSeparatricesOn = 0; // set Show Separatrices to false;
        return;
    }

    this->ui->VF_Window->ShowSeparatricesOn = 1; // set ShowSeparatricesOn to true/1
}


void MainWindow::on_display_periodic_orbits_clicked(bool checked)
{
    if(!checked){
        this->ui->VF_Window->ShowPeriodicOrbitsOn = 0; // set Show Separatrices to false;
        return;
    }


    this->ui->VF_Window->ShowPeriodicOrbitsOn = 1;
}


void MainWindow::on_show_morse_sets_clicked(bool checked)
{
    if(!checked){
        this->ui->VF_Window->ShowSCCsOn = 0; // set Show_strongly_connected_components to flase
        return;
    }

    this->ui->VF_Window->ShowSCCsOn = 1;
}


void MainWindow::on_show_real_ID_clicked(bool checked)
{
    this->set_MCGOn(false);
    showRealIndex = checked;
    this->set_MCGOn(true);
}


void MainWindow::on_Show_Conley_MCG_clicked(bool checked)
{
    this->set_MCGOn(false);
    this->set_ShowConleyCircle(checked);
    this->set_MCGOn(true);
}


void MainWindow::on_ComputeButton_clicked()
{
    this->set_MCGOn(false);
    this->set_VFDisplayOn(false);

    g_start = clock();

    QString str = "----------------";
    print_debug_message(str);
    str = "Start performing Morse Decomposition with new tau value...";
    print_debug_message(str);

    {
        // get the new tau value from the user input
        QString new_tau_str = ui->use_a_tau_LineEdit->text();
        double new_tau = new_tau_str.toDouble();
        global_tau = new_tau; // save new_tau in global_tau

        // get the new error threshold from the user input
        QString new_error_threshold_str = ui->error_threshold_LineEdit->text();
        edge_sample_error = new_error_threshold_str.toDouble();

        if(object->slist.nsingularities == 0)
            object->capture_Singularities();

        if(morse_decomp != NULL)
            delete morse_decomp;

        morse_decomp = new MorseDecomp(); /*initialize the Morse Decomposition component*/

        /*record the performance*/
        morse_decomp->morse_decomp_tau(global_tau); // perform morse decomposition using the new tau value
    }

    g_finish = clock();
    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;

    str =  QString( "Finish computing FG and extracting SCCs. There are %1 edges in FG").arg( morse_decomp->dg->elist->nedges );
    print_debug_message(str);

    str =  QString( "time for the Morse Decomposition is %1 seconds." ).arg( time );
    print_debug_message( str );
    print_debug_message( "----------------" );

    // re-calculate MCG since the tau has been changed
    if(mcg != NULL)
        delete mcg;

    mcg = new MCG_Graph();
    mcg->init_MCG();
    mcg->build_mcg();


    this->set_MCGOn(true);
    this->set_VFDisplayOn(true);
    this->ui->VF_Window->ShowSCCsOn = 1;
    this->ui->show_morse_sets->setChecked(true);
}


void MainWindow::on_construct_minimal_MCG_clicked(bool checked)
{
    this->set_MCGOn(false);
    RemRedundantMCGEdges = checked;
    this->set_MCGOn(true);
}


void MainWindow::on_remove_diconnected_clicked(bool checked)
{
    this->set_MCGOn(false);
    RemoveDisconnMSOn = checked;
    this->set_MCGOn(true);
}


// refine rhe certain morese set with id with a specified id
void MainWindow::on_Refine_button_clicked()
{
    this->set_MCGOn(false);

    QString morse_ID_str = this->ui->local_refine_ID_LineEdit->text();
    int local_morse_id = morse_ID_str.toInt();

    QString tau_str = this->ui->local_refine_tau_LineEdit->text();
    double local_morse_tau = tau_str.toDouble();

    g_start = clock();
    print_debug_message( "----------------" );
    print_debug_message( "Start performing Local Morse Decomposition..." );

    {
        regiontau->init_local_graph_2();   // added by Guoning on 07/06/2010
        regiontau->init_tau_pts_at_tris();
        regiontau->init_tau_pts_at_verts();
        regiontau->RefineRegion(local_morse_id, local_morse_tau);
    }

    g_finish = clock();
    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;

    QString str = QString("time for this process is %1 seconds").arg(time);
    print_debug_message(str);

    this->set_MCGOn(true);
    this->set_VFDisplayOn(true);
}


void MainWindow::on_Euler_1_clicked(bool checked)
{
    Integrator_opt = 0;
}


void MainWindow::on_RK2_clicked(bool checked)
{
    Integrator_opt = 1;
}


void MainWindow::on_RK4_clicked(bool checked)
{
    Integrator_opt = 2;
}


void MainWindow::on_Flip_normal_clicked(bool checked)
{
    this->ui->VF_Window->FlipNormalOn = checked;
}


void MainWindow::on_Grey_texture_clicked(bool checked)
{
    EnGreyTexture = checked;
    this->ui->VF_Window->makePatterns();
}



void MainWindow::on_Disable_lighting_clicked(bool checked)
{
    DisableLighting = checked;
}


void MainWindow::on_Color_map_of_VF_magnitude_clicked(bool checked)
{
    this->ui->VF_Window->ShowColorVFMagOn = checked;
}


void MainWindow::on_IBFV_off_clicked(bool checked)
{
    this->ui->VF_Window->IBFVOff = checked;
}


void MainWindow::on_Place_Streamlines_Button_clicked()
{
    this->set_VFDisplayOn(false);

    print_debug_message( "----------------" );
    print_debug_message( "Start placing streamlines..." );
    g_start = clock();

    /* starting calculate streamlines */

    // get streamline separation distance
    QString sep_str = this->ui->separation_LineEdit->text();
    StreamlinePlace_dsep = sep_str.toDouble();

    // get shorest length
    QString slength_str = this->ui->shortest_LineEdit->text();
    StreamlinePlace_slength = slength_str.toDouble();

    if(evenplace == nullptr){
        evenplace = new EvenStreamlinePlace();
    }

    evenplace->init();
    evenplace->reset_placement();

    if(object->slist.nsingularities == 0){
        object->capture_Singularities();
    }
    evenplace->place_streamlines(0);
    ndisplay_trajs = evenplace->evenstreamlines->ntrajs;


    /* calculating is done */
    g_finish = clock();
    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;
    print_debug_message( QString("Done, there are total %1 streamlines").arg(ndisplay_trajs) );
    QString str = QString("time for placing streamlines is %1 seconds").arg(time);
    print_debug_message( str );
    print_debug_message( "----------------" );
    this->ui->VF_Window->EvenStreamlinePlacement = true;
    this->ui->display_streamlines->setChecked(true);
    this->set_VFDisplayOn(true);
}


void MainWindow::on_visualize_sample_points_clicked(bool checked)
{
    this->ui->VF_Window->ShowEdgeSamplesOn = checked;
}


void MainWindow::on_show_triangle_mapping_clicked(bool checked)
{
    this->ui->VF_Window->ShowTriMappingOn = checked;
}


void MainWindow::on_show_backward_clicked(bool checked)
{
    this->ui->VF_Window->ShowBackward = checked;
}


void MainWindow::on_Change_edge_Button_clicked(bool checked)
{
    if(ui->VF_Window->sampling_edge == 2)
        ui->VF_Window->sampling_edge = 0;
    else
        ui->VF_Window->sampling_edge++;
}


void MainWindow::on_Animate_clicked(bool checked)
{
    this->ui->VF_Window->MoveOrStop = checked;
}


void MainWindow::on_display_connection_region_clicked(bool checked)
{
    if(checked){
        this->ui->VF_Window->ShowConnectionRegion = 1;
    }else{
        this->ui->VF_Window->ShowConnectionRegion = 0;
    }
}


void MainWindow::on_Detect_FixedPts_clicked()
{
    this->set_VFDisplayOn(false);

    QString str = "----------------\n";
    g_start = clock();
    str.append("Start detecting fixed points...\n");
    this->ui->VF_Window->detect_FixedPts();
    str.append("Finish detecting fixed points...\n");
    int num_fixedPts = this->ui->VF_Window->get_num_fixedPts();
    str.append(QString("%1 fixed points have been detected\n").arg(num_fixedPts));
    g_finish = clock();
    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;
    str.append( QString( "time for this operation is %1 seconds\n").arg( time ) );
    str.append( "----------------" );
    this->print_debug_message(str);

    this->on_display_fixed_points_clicked(true);
    this->set_VFDisplayOn(true);
    this->ui->display_fixed_points->setChecked(true);
}


void MainWindow::on_Compute_sep_clicked()
{
    this->set_ECGOn(false);
    this->set_VFDisplayOn(false);

    g_start = clock();
    QString str = "----------------\n";
    str.append("Start computing the separatrices...\n");

    {
        // constructing fixedPts
        if( this->ui->display_fixed_points->isChecked() == false ){
            this->ui->VF_Window->detect_FixedPts();
        }


        // constructing separatrices
        object->cal_separatrices();

        // constructing ECGs
        if(ecg != nullptr)
            delete ecg;
        if(periodic_orbits == nullptr)
            ecg = new ECG_Graph(object->slist.nsingularities);
        else
            ecg = new ECG_Graph(object->slist.nsingularities + periodic_orbits->nporbits);

        ecg->init_ECG();
        ecg->build_ecg();
    }


    g_finish = clock();

    int num_sep = separatrices->ntrajs;
    str.append(QString("%1 separatrices have been detected\n").arg(num_sep));

    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;
    str.append( QString( "time for this operation is %1\n").arg( time ) );
    str.append( "----------------" );
    this->print_debug_message(str);

    this->set_ECGOn(true);

    this->on_display_separatrices_clicked(true);
    this->set_VFDisplayOn(true);
    this->ui->display_separatrices->setChecked(true);
}


void MainWindow::on_Extract_POs_clicked()
{
    // disable ECG
    this->set_ECGOn(false);
    // disable MCG
    this->set_MCGOn(false);
    this->set_VFDisplayOn(false);


    g_start = clock();
    QString str = "----------------\n";
    str.append("Start extracting periodic orbits...\n");

    {
        // if the user clicks display POs without any fixed Pts detected...
        if(object->slist.nsingularities == 0){
            this->ui->VF_Window->detect_FixedPts();
        }

        // extract POs
        object->detect_PeriodicOrbit();

        // construct MCG
        if(mcg != NULL)
            delete mcg;

        mcg = new MCG_Graph();
        mcg->init_MCG();
        mcg->build_mcg();

        // construct separatrices
        object->cal_separatrices();

        //construct ECG
        if(ecg != nullptr)
            delete ecg;

        ecg = new ECG_Graph(object->slist.nsingularities+periodic_orbits->nporbits);
        ecg->init_ECG();
        ecg->build_ecg();
    }


    g_finish = clock();
    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;

    str.append( QString("%1 periodic orbits have been detected.\n").arg(periodic_orbits->nporbits) );
    str.append( QString("time for this operation is %1 seconds.\n").arg(time) );
    str.append( "----------------" );
    this->print_debug_message(str);

    // enable MCG
    this->set_MCGOn(true);
    // enable ECG
    set_ECGOn(true);

    this->on_display_periodic_orbits_clicked(true);
    this->set_VFDisplayOn(true);
    this->ui->display_periodic_orbits->setChecked(true);
}


void MainWindow::on_Morse_Decomp_clicked()
{
    this->set_MCGOn(false);
    this->set_VFDisplayOn(false);

    if(object->slist.nsingularities == 0)
        object->capture_Singularities();

    g_start = clock();

    QString str = "----------------\n";
    str.append("Start performing Morse Decomposition...\n");

    /*perform Morse decompostion*/
    if(morse_decomp != NULL)
        delete morse_decomp;

    morse_decomp = new MorseDecomp(); /*initialize the Morse Decomposition component*/
    morse_decomp->morse_decomp();

    /*compute the MCG*/
    {
        if(mcg != NULL)
            delete mcg;
        mcg = new MCG_Graph();
        mcg->init_MCG();
        mcg->build_mcg();
    }

    g_finish = clock();
    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;

    str.append( QString("%1 Morse sets have been computed.\n").arg( mcg->nlist->nmnodes ) );
    str.append( QString("time for the Morse Decomposition is %1 seconds.\n").arg( time ) );
    str.append( "----------------" );
    this->print_debug_message(str);

    this->set_MCGOn(true);
    this->set_VFDisplayOn(true);

    this->on_show_morse_sets_clicked(true);
    this->ui->show_morse_sets->setChecked(true);
}


void MainWindow::on_display_streamlines_clicked(bool checked)
{
    this->ui->VF_Window->EvenStreamlinePlacement = checked;
}

void MainWindow::on_compute_region_clicked()
{
    if(this->ui->MCG_Window->ShowMCGOn == false) return;

    // disable mcg
    this->set_MCGOn(false);
    this->set_VFDisplayOn(false);

    if(Cal_Regions == true){
        /*compute the MCG*/
        {
            if(mcg != NULL)
                delete mcg;
            mcg = new MCG_Graph();
            mcg->init_MCG();
            mcg->build_mcg();
        }
    }

    this->set_MCGOn(true);

    this->on_display_connection_region_clicked(true);
    this->ui->display_connection_region->setChecked(true);
    this->set_VFDisplayOn(true);
}


void MainWindow::on_No_T_MAX_clicked(bool checked)
{
    noTmax = checked;
}


void MainWindow::on_Auto_Refine_Btn_clicked()
{
    this->set_MCGOn(false);
    this->set_VFDisplayOn(false);

    QString min_pri_str = this->ui->min_priority_LineEdit->text();
    this->ui->VF_Window->min_pri = min_pri_str.toInt();

    QString tau_max_str = this->ui->tau_max_LineEdit->text();
    this->ui->VF_Window->tau_max = tau_max_str.toDouble();

    QString iter_max_str = this->ui->Iter_max_LineEdit->text();
    this->ui->VF_Window->iter_max = iter_max_str.toDouble();

    QString initial_tau = this->ui->Inital_LineEdit->text();
    this->ui->VF_Window->init_tau = initial_tau.toDouble();


    g_start = clock();

    QString str = "----------------\n";
    str.append("Start performing Auto Morse Decomposition...\n");

    for(int i = 0; i<object->tlist.ntris; i++){
        object->tlist.tris[i]->used_tau = 0;
    }

    if(!noTmax){
        regiontau->init_local_graph_2();
        regiontau->init_tau_pts_at_tris();
        regiontau->init_tau_pts_at_verts();
        regiontau->Auto_Refine(
                    this->ui->VF_Window->tau_max,
                    this->ui->VF_Window->init_tau,
                    this->ui->VF_Window->min_pri,
                    this->ui->VF_Window->iter_max);
    }
    else{
        regiontau->init_local_graph_2();
        regiontau->init_tau_pts_at_tris();
        regiontau->init_tau_pts_at_verts();
        regiontau->Auto_Refine(
                    this->ui->VF_Window->init_tau,
                    this->ui->VF_Window->min_pri,
                    4);
    }

    double tau_min = 1e10;
    double tau_max = -1e10;

    for(int i = 0; i < object->tlist.ntris; i++){
        if(object->tlist.tris[i]->used_tau > tau_max)
            tau_max = object->tlist.tris[i]->used_tau;

        if(object->tlist.tris[i]->used_tau < tau_min)
            tau_min = object->tlist.tris[i]->used_tau;
    }

    g_finish = clock();
    double time = (double)(g_finish-g_start)/CLOCKS_PER_SEC;

    str.append( QString("time for this Auto Refinemet is %1 seconds.\n").arg( time ) );

    int n_mcg_nodes = 0;
    for(int i = 0; i<mcg->nlist->nmnodes; i++){
        if(mcg->nlist->mnodes[i]->cancelled) continue;
        n_mcg_nodes ++;
    }

    int n_mcg_edges = 0;
    for(int i = 0; i<mcg->elist->nedges; i++){
        if(mcg->elist->edges[i]->cancel) continue;
        n_mcg_edges ++;
    }

    str.append( QString("Finish computing the MCG. It contains %1 nodes and %2 edges\n").arg(n_mcg_nodes).arg(n_mcg_edges) );
    str.append( QString("The maximum tau value is %1\n").arg(tau_max) );
    str.append( QString("The minimum tau value is %1\n").arg(tau_min) );

    str.append( "----------------" );
    this->print_debug_message(str);


    // enable MCG and show the morse sets
    this->set_MCGOn(true);
    this->set_VFDisplayOn(true);
    this->on_show_morse_sets_clicked(true);
    this->ui->show_morse_sets->setChecked(true);
}


void MainWindow::on_Browsers_Button_clicked()
{
    QString dataset_filePath = QString(dataset_folder);


    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    dataset_filePath,
                                                    tr("models (*.ply)"));
    // if the user doesn't select any file...
    if(QString::compare(fileName, "", Qt::CaseInsensitive) == 0)
        return;

    this->set_MCGOn(false);
    this->set_ECGOn(false);
    this->set_VFDisplayOn(false);

    // restore ui and flags
    this->on_display_connection_region_clicked(false);
    this->ui->display_connection_region->setChecked(false);

    this->on_No_T_MAX_clicked(false);
    this->ui->No_T_MAX->setChecked(false);

    this->on_visualize_sample_points_clicked(false);
    this->ui->visualize_sample_points->setChecked(false);

    this->on_show_triangle_mapping_clicked(false);
    this->ui->show_triangle_mapping->setChecked(false);

    this->on_show_backward_clicked(false);
    this->ui->show_backward->setChecked(false);

    this->on_Flip_normal_clicked(false);
    this->ui->Flip_normal->setChecked(false);

    this->on_Grey_texture_clicked(false);
    this->ui->Grey_texture->setChecked(false);

    this->on_Disable_lighting_clicked(false);
    this->ui->Disable_lighting->setChecked(false);

    this->on_Animate_clicked(false);
    this->ui->Animate->setChecked(false);

    this->on_Color_map_of_VF_magnitude_clicked(false);
    this->ui->Color_map_of_VF_magnitude->setChecked(false);

    this->on_IBFV_off_clicked(false);
    this->ui->IBFV_off->setChecked(false);

    this->on_display_streamlines_clicked(false);
    this->ui->display_streamlines->setChecked(false);

    this->on_display_fixed_points_clicked(false);
    this->ui->display_fixed_points->setChecked(false);

    this->on_display_separatrices_clicked(false);
    this->ui->display_separatrices->setChecked(false);

    this->on_display_periodic_orbits_clicked(false);
    this->ui->display_periodic_orbits->setChecked(false);

    this->on_show_morse_sets_clicked(false);
    this->ui->show_morse_sets->setChecked(false);

    this->on_show_real_ID_clicked(false);
    this->ui->show_real_ID->setChecked(false);

    this->on_Show_Conley_MCG_clicked(false);
    this->ui->Show_Conley_MCG->setChecked(false);

    this->on_construct_minimal_MCG_clicked(false);
    this->ui->construct_minimal_MCG->setChecked(false);

    this->on_remove_diconnected_clicked(false);
    this->ui->remove_diconnected->setChecked(false);

    this->on_Euler_1_clicked(true);
    this->ui->Euler_1->setChecked(true);

    ui->VF_Window->initializeGL2(fileName);
}
