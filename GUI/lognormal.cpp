#include "mainwindow.h"
#include "log.h"
#include <math.h>
#include <random>
#include <iostream>

LOG_USE();

std::uniform_real_distribution<double> dist(0.0, 1.0);
std::mt19937 mt(12345);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
double random_uni()
{
//    std::random_device rd;
//    std::mt19937 mt(rd());
//    std::uniform_real_distribution<double> dist(0.0, 1.0);

//    std::cout << dist(mt) << "\n";
//    for (int i=0; i<16; ++i)
//        std::cout << dist(mt) << "\n";
//    }

    return dist(mt);
//    return rand()/(RAND_MAX + 1.0);
}

//-----------------------------------------------------------------------------------------
// from: http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
// Implements the Gauss error function.
//   erf(z) = 2 / sqrt(pi) * integral(exp(-t*t), t = 0..z)
//
// fractional error in math formula less than 1.2 * 10 ^ -7.
// although subject to catastrophic cancellation when z in very close to 0
// from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
//-----------------------------------------------------------------------------------------
double erf(double z)
{
   double t = 1.0 / (1.0 + 0.5 * fabs(z));
   // use Horner's method
   double ans = 1 - t * exp( -z*z -  1.26551223 +
         t * ( 1.00002368 +
         t * ( 0.37409196 +
         t * ( 0.09678418 +
         t * (-0.18628806 +
         t * ( 0.27886807 +
         t * (-1.13520398 +
         t * ( 1.48851587 +
         t * (-0.82215223 +
         t * ( 0.17087277))))))))));
   if (z >= 0.0)
       return ans;
   else
       return -ans;
}

//-----------------------------------------------------------------------------------------
// When X is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double pnorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;

    z1 = (x1-mu)/sig;
    e1 = erf(z1/sqrt(2.0))/2;
    z2 = (x2-mu)/sig;
    e2 = erf(z2/sqrt(2.0))/2;
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// When log(X) is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double plognorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;

    z1 = 0;
    z2 = 0;
    if (x1 == 0)
        e1 = -0.5;
    else {
        z1 = (log(x1)-mu)/sig;
        e1 = erf(z1/sqrt(2.0))/2;
    }
    if (x2 == 0)
        e2 = -0.5;
    else {
        z2 = (log(x2)-mu)/sig;
        e2 = erf(z2/sqrt(2.0))/2;
    }
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// Create the lognormal distribution with median = p1, shape = p2
// at n points stored in x[], probability values stored in prob[].
// Note that x[0] = 0.
// The range of x is currently just less than 4*median.  This should be
// OK for values of shape < 2.
// Convert probability into probability density
//-----------------------------------------------------------------------------------------
void MainWindow::create_lognorm_dist(double p1, double p2, int n, double *x, double *prob)
{
//    bool use_exp = false;

//    if (!use_exp) {
    double xmax, dx, mu_l, sig_l, x1, x2;

    if (p1 >= 0.5)
        xmax = p1*4;
    else
        xmax = p1*8;

    dx = xmax/n;
    mu_l = log(p1);
    sig_l = log(p2);
    for (int ix=0; ix<n; ix++) {
        x1 = (ix - 0.5)*dx;
        x2 = x1 + dx;
        x1 = max(x1,0.0);
        x[ix] = (x1+x2)/2;
        prob[ix] = plognorm(x1,x2,mu_l,sig_l)/(x2-x1);
    }
}

//------------------------------------------------------------------------------------------------------
// Distribution of sum of nexp exponentially distributed RVs, with rate constants lambda[].
// Note: lambda = 1/mean
//------------------------------------------------------------------------------------------------------
void MainWindow::create_expon_dist(double tbase, double mean[], int n, double *x, double *prob)
{
    int i, j, nexp;
    double lambda[3], a1, a2, a3, z, dR;

    nexp = 0;
    for (i=0; i<3; i++) {
        if (mean[i] > 0) {
            lambda[nexp] = 1/mean[i];
            nexp++;
        }
    }
    if (nexp == 1) {
        dR = (1/lambda[0])/20;
    } else if (nexp == 2) {
        dR = (1/lambda[0] + 1/lambda[1])/40;
    } else if (nexp == 3) {
        dR = (1/lambda[0] + 1/lambda[1] + 1/lambda[2])/60;
    }
    for (j=0; j<n; j++) {
        prob[j] = 0;
        x[j] = tbase + j*dR;
    }
    prob[0] = 0;
    // Analytical
    if (nexp == 1) {
        a1 = lambda[0];
        for (j=1; j<n; j++) {
            z = x[j] - tbase;
            prob[j] = a1*exp(-a1*z);
        }
    } else if (nexp == 2) {
        // two cases
        if (lambda[0] == lambda[1]) {
            a1 = lambda[0];
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = a1*a1*z*exp(-a1*z);
            }
        } else {
            a1 = lambda[0];
            a2 = lambda[1];
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = a1*a2*(exp(-a1*z) - exp(-a2*z))/(a2-a1);
            }
        }
    } else if (nexp == 3) {
        // three cases: all different, two the same, all the same
        if (lambda[0] == lambda[1] && lambda[0] == lambda[2]) {
            a1 = lambda[0];
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = a1*a1*a1*z*z*exp(-a1*z)/2;
            }
        } else if (lambda[0] != lambda[1] && lambda[0] != lambda[2] && lambda[1] != lambda[2]) {
            a1 = lambda[0];
            a2 = lambda[1];
            a3 = lambda[2];
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = a1*a2*a3*(exp(-a1*z)/((a2-a1)*(a3-a1)) +
                                    exp(-a2*z)/((a1-a2)*(a3-a2)) +
                                    exp(-a3*z)/((a1-a3)*(a2-a3)));
            }
        } else {
            if (lambda[0] == lambda[1]) {
                a1 = lambda[0];
                a3 = lambda[2];
            } else if (lambda[0] == lambda[2]) {
                a1 = lambda[0];
                a3 = lambda[1];
            } else {
                a1 = lambda[1];
                a3 = lambda[0];
            }
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = (a1*a1*a3/((a3-a1)*(a3-a1)))*((z*(a3-a1)-1)*exp(-a1*z) + exp(-a3*z));
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: drawDistPlots(bool dummy)
{
    double *x, *prob;
    x = new double[nDistPts];
    prob = new double[nDistPts];
    QString median_qstr, shape_qstr;
    double median, shape;
    double expMean[3];
    double expTbase;
    bool use_lognormal = cbox_USE_LOGNORMAL_DIST->isChecked();

    groupBox_divisiondistributions->setEnabled(use_lognormal);
    for (int j=0; j<ndistplots; j++) {
        QwtPlot *qp = distplot_list[j];
        if (j == 0) {
            qp->setTitle("Type 1 division time (hrs)");
            if (use_lognormal) {
                median_qstr = line_DIVIDE_TIME_1_MEDIAN->text();
                shape_qstr = line_DIVIDE_TIME_1_SHAPE->text();
            } else {
                expTbase = line_T_G1_1->text().toDouble() + line_T_S_1->text().toDouble() + line_T_G2_1->text().toDouble();
                expMean[0] = line_G1_MEAN_DELAY_1->text().toDouble();
                expMean[1] = line_S_MEAN_DELAY_1->text().toDouble();
//                expMean[2] = line_G2_MEAN_DELAY_1->text().toDouble();
                expMean[2] = line_G2_DELAY_FACTOR_1->text().toDouble();
            }
        } else if (j == 1) {
            qp->setTitle("Type 2 division time (hrs)");
            if (use_lognormal) {
                median_qstr = line_DIVIDE_TIME_2_MEDIAN->text();
                shape_qstr = line_DIVIDE_TIME_2_SHAPE->text();
            } else {
                expTbase = line_T_G1_2->text().toDouble() + line_T_S_2->text().toDouble() + line_T_G2_2->text().toDouble();
                expMean[0] = line_G1_MEAN_DELAY_2->text().toDouble();
                expMean[1] = line_S_MEAN_DELAY_2->text().toDouble();
//                expMean[2] = line_G2_MEAN_DELAY_2->text().toDouble();
                expMean[2] = line_G2_DELAY_FACTOR_2->text().toDouble();
            }
        } else if (j == 2) {
            qp->setTitle("Apoptosis delay time 1 (hrs)");
            median_qstr = line_APOPTOSIS_MEDIAN_1->text();
            shape_qstr = line_APOPTOSIS_SHAPE_1->text();
        } else if (j == 3) {
            qp->setTitle("Apoptosis delay time 2 (hrs)");
            median_qstr = line_APOPTOSIS_MEDIAN_2->text();
            shape_qstr = line_APOPTOSIS_SHAPE_2->text();
        }
        if (use_lognormal || j >= 2) {
            median = median_qstr.toDouble();
            shape = shape_qstr.toDouble();
            create_lognorm_dist(median, shape, nDistPts, x, prob);
        } else {
            create_expon_dist(expTbase, expMean, nDistPts, x, prob);
        }

        int n = dist_limit(prob,nDistPts);
        if (j >= 2) {
            for (int i=0; i<=n; i++)
                x[i] += 12;
        }
        double xmax = x[n];
        qp->setAxisScale(QwtPlot::xBottom, 0.0, xmax, 0.0);

        if (first_plot) {
            curve_list[j] = new QwtPlotCurve("title");
            curve_list[j]->setSamples(x, prob, n);
            curve_list[j]->attach(qp);
        } else {
            curve_list[j]->setSamples(x, prob, n);
        }
        qp->replot();
    }
    delete [] x;
    x = NULL;
    delete [] prob;
    prob = NULL;
    first_plot = false;
}
