#ifndef DISTANCE_H
#define DISTANCE_H

#include "utils.h"



namespace DIS {



    void d_PP(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        double& d);

    void g_PP(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        Vector6d& g);

    void H_PP(Matrix6d& H);

    void derivTest_PP(void);
    

    void d_PE(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        double& d);

    void g_PE(double v01, double v02, double v03, double v11, double v12, double v13, double v21, double v22, double v23, double g[9]);
    

    void g_PE(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        Vector9d& g);


    void H_PE(double v01, double v02, double v03, double v11, double v12, double v13,
        double v21, double v22, double v23, double H[81]);
    

    void H_PE(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        Matrix9d& H);


    void derivTest_PE(void);
    

    void d_PT(const Eigen::Vector3d& v0,
            const Eigen::Vector3d& v1,
            const Eigen::Vector3d& v2,
            const Eigen::Vector3d& v3,
            double& d);
    

    void g_PT(double v01, double v02, double v03, double v11, double v12, double v13,
            double v21, double v22, double v23, double v31, double v32, double v33,
            double g[12]);
    

    void g_PT(const Eigen::Vector3d& v0,
            const Eigen::Vector3d& v1,
            const Eigen::Vector3d& v2,
            const Eigen::Vector3d& v3,
            Vector12d& g);
   

    void H_PT(double v01, double v02, double v03, double v11, double v12, double v13,
            double v21, double v22, double v23, double v31, double v32, double v33,
            double H[144]);
    

    void H_PT(const Eigen::Vector3d& v0,
            const Eigen::Vector3d& v1,
            const Eigen::Vector3d& v2,
            const Eigen::Vector3d& v3,
            Matrix12d& H);
    

    void derivTest_PT(void);
    

    void d_EE(const Eigen::Vector3d& v0,
            const Eigen::Vector3d& v1,
            const Eigen::Vector3d& v2,
            const Eigen::Vector3d& v3,
            double& d);


    void g_EE(double v01, double v02, double v03, double v11, double v12, double v13,
                double v21, double v22, double v23, double v31, double v32, double v33,
                double g[12]);
    

    void g_EE(const Eigen::Vector3d& v0,
                const Eigen::Vector3d& v1,
                const Eigen::Vector3d& v2,
                const Eigen::Vector3d& v3,
                Vector12d& g);


    void H_EE(double v01, double v02, double v03, double v11, double v12, double v13,
                    double v21, double v22, double v23, double v31, double v32, double v33,
                    double H[144]);
    

    void H_EE(const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    Matrix12d& H);
    

     void derivTest_EE(void);
    

    // http://geomalgorithms.com/a07-_distance.html
     int dType_EE(const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3);
    

     int dType_PT(const Eigen::Vector3d& v0_,
                    const Eigen::Vector3d& v1_,
                    const Eigen::Vector3d& v2_,
                    const Eigen::Vector3d& v3_);
    

     void checkDType(void);
    

    // http://geomalgorithms.com/a02-_lines.html
     void computePointEdgeD(const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    double& d);
    

     void computePointTriD(const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    double& d);
    

     void computeEdgeEdgeD(const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    double& d);
    

     void checkEdgeEdgeD(void);
    




    // for fixing tiny step size issue
     void computeEECrossSqNorm(
                    const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    double& result);
    

    // cal_EEM_c_wrt_x
     void computeEECrossSqNormGradient(double v01, double v02, double v03, double v11,
                    double v12, double v13, double v21, double v22, double v23, double v31, double v32, double v33, double g[12]);
    

     void computeEECrossSqNormGradient(
                    const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    Vector12d& grad);
    

    // cal_EEM_2c_wrt_x2
     void computeEECrossSqNormHessian(double v01, double v02, double v03, double v11,
                    double v12, double v13, double v21, double v22, double v23, double v31, double v32, double v33, double H[144]);
    

     void computeEECrossSqNormHessian(
                    const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    Matrix12d& Hessian);
    

     void derivTest_EECross(void);
    

     void compute_q(double input, double eps_x, double& e);
    

     void compute_q_g(double input, double eps_x, double& g);
    

     void compute_q_H(double input, double eps_x, double& H);
    

     void compute_e(
                    const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    double eps_x, double& e);
    

     void compute_e_g(
                    const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    double eps_x, Vector12d& g);
    

     void compute_e_H(
                    const Eigen::Vector3d& v0,
                    const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2,
                    const Eigen::Vector3d& v3,
                    double eps_x, Matrix12d& H);
    

     void derivTest_e(double eps_x = 10.0);
    


     double cal_EEM_eps_x(Eigen::Vector3d& P1_Rest, Eigen::Vector3d& P2_Rest, Eigen::Vector3d& Q1_Rest, Eigen::Vector3d& Q2_Rest);





     /////////////////////////////////////////////////////////////////////////
     // The following codes calculate the edge-edge mollifier
     /////////////////////////////////////////////////////////////////////////

     // calculate the value of eps_x
     double cal_EEM_eps_x(Eigen::Vector3d& P1_Rest, Eigen::Vector3d& P2_Rest, Eigen::Vector3d& Q1_Rest, Eigen::Vector3d& Q2_Rest);


     // calculate the mollifier c
     double cal_EEM_c(Eigen::Vector3d& P1, Eigen::Vector3d& P2, Eigen::Vector3d& Q1, Eigen::Vector3d& Q2);


     // calculate the mollifier value wrt c
     double cal_EEM_ek_wrt_c(double eps_x, double c);


     // calculate the second-order derivative of mollifier value wrt c
     double cal_EEM_2ek_wrt_c2(double eps_x, double c);


     // calculate the mollifier c wrt x
     Vector12d cal_EEM_c_wrt_x(Eigen::Vector3d& P1, Eigen::Vector3d& P2, Eigen::Vector3d& Q1, Eigen::Vector3d& Q2);
     

     // calculate the second-order derivative of mollifier c wrt x
     Matrix12d cal_EEM_2c_wrt_x2(Eigen::Vector3d& P1, Eigen::Vector3d& P2, Eigen::Vector3d& Q1, Eigen::Vector3d& Q2);
     

 

} // namespace IPC






#endif