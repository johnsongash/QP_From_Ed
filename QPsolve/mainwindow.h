#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <iostream>
#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <qcustomplot.h>
#include "qptrajectory.h"

//#ifdef CGAL_USE_GMP
//#include <CGAL/Gmpz.h>
//typedef CGAL::Gmpz ET;
//#else
//#include <CGAL/MP_Float.h>
//typedef CGAL::MP_Float ET;
//#endif
//#include <CGAL/MP_Float.h>
//typedef CGAL::MP_Float ET;
#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    QCPCurve *newCurve;// = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
    QVector<QCPCurveData> data;
    QCPCurve *newCurve2;
    QVector<QCPCurveData> data2;
    QVector<QCPCurveData> data3;
    QCPCurve *newCurve3;

    QCPCurve *path_curve ;
    QVector<QCPCurveData> path_data;
    QCPCurve *payload_curve ;
    QVector<QCPCurveData> payload_data;


    double ax,ay;
    double vx,vy;
    double px,py;
    qptrajectory *plan;
    std::vector<trajectory_profile> profile;
    double t;
    QTimer *time;
    double loop_size;
    double loop_count;
    QCPCurve *circle ;//= new QCPCurve(ui->customplot->xAxis, ui->customplot->yAxis);
    QVector< QCPCurveData > circle_data;
    QCPCurve *follower ;//= new QCPCurve(ui->customplot->xAxis, ui->customplot->yAxis);
    QVector< QCPCurveData > follower_data;
    QCPItemLine *item;
    QCPItemLine *link;



    double ax_2;
    double ay_2;
    double vx_2 , vy_2;
    double px_2 , py_2;

    //impedance coefficient
    double mass ;
    double damping;
    double spring;

    //link coefficient

    double l_mass;
    double l_damping;
    double l_spring;

    double last_p , last_v;
    double theta;
    double forcex , forcey ;
    double force_gain;
    double interaction_time;
    double mouse_corx ;
    double mouse_cory ;


    //link
    double fint_x , fint_y ;
    double u_lx , u_ly;
    double length ;
    //follower
    double lm;
    double lpx;
    double lpy;
    double lvx;
    double lvy;
    double lax;
    double lay;
private:
    Ui::MainWindow *ui;

private slots:
    void mouseClick(QMouseEvent* event);

    void update();

};

#endif // MAINWINDOW_H
