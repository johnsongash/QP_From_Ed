#include <qptrajectory.h>
qptrajectory::qptrajectory(){




}
qptrajectory::~qptrajectory(){

}



std::vector<double>
qptrajectory::qpsovle(profile begin, profile end, double time_interval, int c){
Program solver(CGAL::EQUAL, false, 0, false, 0);
//CGAL::Quadratic_program_options options;
//options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule



std::vector<double> polynomial;
polynomial.clear();
double t =time_interval ;
double b0 = 24,
       b1 = 120,
       b2 = 360,
       b3 = 840;

Eigen::MatrixXd D(8,8);
Eigen::MatrixXd d(4,4);
Eigen::MatrixXd A(8,8);
Eigen::MatrixXd Al(8,8);
Eigen::MatrixXd As(8,8);
Eigen::MatrixXd B(8,1);
//std::cout << "here" <<std::endl;
D.setZero();
d.setZero();
A.setZero();
Al.setZero();
As.setZero();
B.setZero();

double d11 = b0*b0*t        , d12 = b0*b1*t*t        , d13 = b0*b2*t*t*t         , d14 = b0*b3*t*t*t*t ;
double d21 = b0*b1*t*t      , d22 = b1*b1*t*t*t      , d23 = b1*b2*t*t*t*t       , d24 = b1*b3*t*t*t*t*t ;
double d31 = b0*b2*t*t*t    , d32 = b2*b1*t*t*t*t    , d33 = b2*b2*t*t*t*t*t     , d34 = b2*b3*t*t*t*t*t*t;
double d41 = b3*b0*t*t*t*t  , d42 = b3*b1*t*t*t*t*t  , d43 = b3*b2*t*t*t*t*t*t   , d44 = b3*b3*t*t*t*t*t*t*t ;


d << (1/1.0)*d11*1.0 , (1/2.0)*d21*1.0 ,(1/3.0)*d31*1.0  ,      (1/4.0)*d41*1.0     ,
     (1/2.0)*d21*1.0 , (1/3.0)*d22*1.0 ,(1/4.0)*d32*1.0  ,      (1/5.0)*d42*1.0     ,
     (1/3.0)*d31*1.0 , (1/4.0)*d32*1.0 , (1/5.0)*d33*1.0 ,     (1/6.0)*d43*1.0    ,
     (1/4.0)*d41*1.0 , (1/5.0)*d42*1.0 , (1/6.0)*d43*1.0 , (1/7.0)*d44*1.0 ;

    D.block<4,4>(4,4) = d;

    //std::cout << "constant" << d <<std::endl;


     Al<< 1 ,  0  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  1  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   2  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   0  ,  6  , 0 , 0 , 0  , 0 ,
         1 , 1*t , 1*t*t , 1*t*t*t , 1*t*t*t*t , 1*t*t*t*t*t , 1*t*t*t*t*t*t , 1*t*t*t*t*t*t*t,
         0 ,  1  ,   2*t , 3*t*t   , 4*t*t*t   , 5*t*t*t*t   , 6*t*t*t*t*t   , 7*t*t*t*t*t*t  ,
         0 ,  0  ,    2 ,    6*t    ,12*t*t    , 20*t*t*t     ,30*t*t*t*t    , 42*t*t*t*t*t   ,
         0 ,  0  ,    0  ,    6     ,24*t      , 60*t*t       ,120*t*t*t     , 210*t*t*t*t    ;

     As<< 1 ,  0  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  1  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   2  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   0  ,  6  , 0 , 0 , 0  , 0 ,
         1 , 1*t , 1*t*t , 1*t*t*t , 1*t*t*t*t , 1*t*t*t*t*t , 1*t*t*t*t*t*t , 1*t*t*t*t*t*t*t,
         0 ,  0  ,   0 ,   0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   0 ,   0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   0 ,   0  , 0 , 0 , 0  , 0 ;

     A<< 1 ,  0  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  1  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   2  ,  0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   0  ,  6  , 0 , 0 , 0  , 0 ,
         1 , 1*t , 1*t*t , 1*t*t*t , 1*t*t*t*t , 1*t*t*t*t*t , 1*t*t*t*t*t*t , 1*t*t*t*t*t*t*t,
         0 ,  0  ,   0 ,   0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   0 ,   0  , 0 , 0 , 0  , 0 ,
         0 ,  0  ,   0 ,   0  , 0 , 0 , 0  , 0 ;

//    Al<< 1 ,  0  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  1  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  0  ,   2  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  0  ,   0  ,  6  , 0 , 0 , 0  , 0 ,
//         1 , 1*t , 1*t*t , 1*t*t*t , 1*t*t*t*t , 1*t*t*t*t*t , 1*t*t*t*t*t*t , 1*t*t*t*t*t*t*t,
//         0 ,  1  ,   2*t , 3*t*t   , 4*t*t*t   , 5*t*t*t*t   , 6*t*t*t*t*t   , 7*t*t*t*t*t*t  ,
//         0 ,  0  ,    2 ,    6*t    ,12*t*t    , 20*t*t*t     ,30*t*t*t*t    , 42*t*t*t*t*t   ,
//         0 ,  0  ,    0  ,    6     ,24*t      , 60*t*t       ,120*t*t*t     , 210*t*t*t*t    ;

//    Al<< 1 ,  0  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  1  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  0  ,   2  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  0  ,   0  ,  6  , 0 , 0 , 0  , 0 ,
//         1 , 1*t , 1*t*t , 1*t*t*t , 1*t*t*t*t , 1*t*t*t*t*t , 1*t*t*t*t*t*t , 1*t*t*t*t*t*t*t,
//         0 ,  1  ,   2*t , 3*t*t   , 4*t*t*t   , 5*t*t*t*t   , 6*t*t*t*t*t   , 7*t*t*t*t*t*t  ,
//         0 ,  0  ,    2 ,    6*t    ,12*t*t    , 20*t*t*t     ,30*t*t*t*t    , 42*t*t*t*t*t   ,
//         0 ,  0  ,    0  ,    6     ,24*t      , 60*t*t       ,120*t*t*t

//    Al<< 1 ,  0  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  1  ,   0  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  0  ,   2  ,  0  , 0 , 0 , 0  , 0 ,
//         0 ,  0  ,   0  ,  6  , 0 , 0 , 0  , 0 ,
//         1 , 1*t , 1*t*t , 1*t*t*t , 1*t*t*t*t , 1*t*t*t*t*t , 1*t*t*t*t*t*t , 1*t*t*t*t*t*t*t,
//         0 ,  1  ,   2*t , 3*t*t   , 4*t*t*t   , 5*t*t*t*t   , 6*t*t*t*t*t   , 7*t*t*t*t*t*t  ,
//         0 ,  0  ,    2 ,    6*t    ,12*t*t    , 20*t*t*t     ,30*t*t*t*t    , 42*t*t*t*t*t   ,
//         0 ,  0  ,    0  ,    6     ,24*t      , 60*t*t       ,120*t*t*t

//         B<< begin.V , end.V;



//     B <<  begin.pos, begin.vel , begin.acc  , begin.jerk ,
//             end.pos,  end.vel  , end.acc    , end.jerk ;
        B<< begin.V , end.V;

     for(int i=0 ; i<8; i++){
         for(int j=0;j<(i+1);j++){
             solver.set_d(i, j, D(i,j));
         }
     }

     if(c==0){
         for(int i=0;i<8;i++){
             for(int j=0;j<8;j++){
                 solver.set_a(j, i,  As(i,j));
             }
         }
     }
     else if(c==2){
         //std::cout << "end" <<std::endl;
         for(int i=0;i<8;i++){
             for(int j=0;j<8;j++){
                 solver.set_a(j, i,  Al(i,j));
             }
         }

     }
     else{
         for(int i=0;i<8;i++){
             for(int j=0;j<8;j++){
                 solver.set_a(j, i,  A(i,j));
             }
         }

     }
     for(int i=0;i<8;i++){
            solver.set_b(i,B(i,0));
     }
    Solution s =CGAL::solve_quadratic_program(solver, ET() );
    assert (s.is_valid());


    for (Solution::Variable_value_iterator it_ =  s.variable_values_begin();
         it_ != s.variable_values_end();
         ++it_) {
        CGAL::Quotient<ET> data = *it_;

      polynomial.push_back(data.numerator().to_double()/data.denominator().to_double());
    }

    return polynomial;
}

std::vector<trajectory_profile> qptrajectory::get_profile(std::vector<segments> seg , double time_interval , double dt ){
    double t=0.0 ;
    std::vector<double> polyx , polyy;
    std::vector<trajectory_profile> tprofile;
    std::vector<trajectory_profile> pprofile;
    tprofile.clear();
    Eigen::Vector3d d(0,0,0);
    trajectory_profile data(d ,d,d ,d,0.01);
    int c;

    for(int i=0 ; i < seg.size() ; i++){

        profile begin ,end;

        if(i==0){
        c=0;
        }
        else if(i==(seg.size()-1)){
        c=2;
        }
        else{
        c=1;
        }
        std::cout << "c=" << c <<std::endl;

        if(c==0){
            begin.V<< seg[i].b_c.pos[0] , seg[i].b_c.vel[0] , seg[i].b_c.acc[0] , seg[i].b_c.jerk[0];
            end.V  << seg[i].t_c.pos[0] , seg[i].t_c.vel[0] , seg[i].t_c.acc[0] , seg[i].t_c.jerk[0];
        }
        //std::cout << "x=" << seg[i].t_c.pos[0] <<std::endl;

        else{
            begin.V<< seg[i].b_c.pos[0] , pprofile[pprofile.size()-1].vel[0] , pprofile[pprofile.size()-1].acc[0] , pprofile[pprofile.size()-1].jerk[0];
            end.V  << seg[i].t_c.pos[0] , seg[i].t_c.vel[0] , seg[i].t_c.acc[0] , seg[i].t_c.jerk[0];
        }
        std::cout <<  "pva:"  << begin.V   <<std::endl;
//        begin.V<< seg[i].b_c.pos[0] , seg[i].b_c.vel[0] , seg[i].b_c.acc[0] , 0;
//        end.V  << seg[i].t_c.pos[0] , seg[i].t_c.vel[0] , seg[i].t_c.acc[0] , 0;

        polyx  = qpsovle(begin , end , seg[i].time_interval, c);

        begin.V.setZero();
        end.V.setZero();
        if(c==0){
            begin.V<< seg[i].b_c.pos[1] , seg[i].b_c.vel[1] , seg[i].b_c.acc[1] , seg[i].b_c.jerk[1];
            end.V  << seg[i].t_c.pos[1] , seg[i].t_c.vel[1] , seg[i].t_c.acc[1] , seg[i].t_c.jerk[1];
        }
        else{
            begin.V<< seg[i].b_c.pos[1] , pprofile[pprofile.size()-1].vel[1] , pprofile[pprofile.size()-1].acc[1] , pprofile[pprofile.size()-1].jerk[1];
            end.V  << seg[i].t_c.pos[1] , seg[i].t_c.vel[1] , seg[i].t_c.acc[1] , seg[i].t_c.jerk[1];
        }
//        begin.V<< seg[i].b_c.pos[1] , seg[i].b_c.vel[1] , seg[i].b_c.acc[1] , 0;
//        end.V<< seg[i].t_c.pos[1] , seg[i].t_c.vel[1] , seg[i].t_c.acc[1] , 0;
        polyy =   qpsovle(begin , end , seg[i].time_interval, c);
        t=0;

        for(int j=0;j<(seg[i].time_interval/dt);j++){
            t =(double) dt*j;
            data.pos << polynomial(polyx , t) ,polynomial(polyy , t) , 0;
            data.vel << polynomial_d1(polyx,t),polynomial_d1(polyy , t) ,0;
            data.acc << polynomial_d2(polyx,t) , polynomial_d2(polyy,t) , 0;
            data.jerk << polynomial_d3(polyx,t) , polynomial_d3(polyy,t) , 0;
            tprofile.push_back(data);
        }
        pprofile.clear();
        pprofile=tprofile;
           for(int i=0 ;i< tprofile.size();i++){
               //std::cout << "step : " <<i<<std::endl;
               //std::cout <<     tprofile[i].pos[0]    <<std::endl;
           }
        //std::cout <<     pprofile[pprofile.size()-1].vel[0]     <<std::endl;
    }

return tprofile;
}

/* Position, A0+A1*t+A2*t^2...... */
double qptrajectory::polynomial(std::vector<double> data ,double t){
    double sum =0.0 , var =1;
    for(int i =0 ; i<data.size();i++){
        sum+= data[i]*var;
        var*=t;
    }
    return sum;
}

double cpow(double t,int times){
    double var =1.0;
    for(int i=0;i<times ;i++){
        var*=t;
    }
    return var;
}

double qptrajectory::polynomial_d1(std::vector<double> data ,double t){
    double sum=0.0;
    for(int i=1;i<data.size();i++){
        sum+=  (double) data[i] * i * cpow( t , i-1);
    }
    return sum;
}
double qptrajectory::polynomial_d2(std::vector<double> data ,double t){
    double sum =0.0 , var =1.0;
    for(int i =2 ; i<data.size();i++){
        sum += data[i] * i *(i-1)* cpow( t , i-2);
    }
    return sum;
}
double qptrajectory::polynomial_d3(std::vector<double> data ,double t){
    double sum =0.0 , var =1.0;
    for(int i = 3 ; i<data.size();i++){
        sum += data[i] * i *(i-1)* (i-2) * cpow( t , i-3);
    }
    return sum;
}
