// FACS and histogram functions

#include <QtGui>

#include "mainwindow.h"
#include "log.h"
#include "params.h"
#include "global.h"
#include "malloc.h"
#include <QPainter>
#include <qwt_scale_engine.h>

#ifdef linux
#include <QTcpServer>
#else
#include <QTcpServer.h>
#endif

LOG_USE();

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
/*
void MainWindow::processGroupBoxClick(QString text)
{
    LOG_QMSG("processGroupBoxClick: " + text);
    QwtPlot *plot;

    if (text.compare("Histo") == 0) {
        LOG_MSG("save Histo plot");
        bool use_HistoBar = (buttonGroup_histotype->checkedId() == 1);
        if (use_HistoBar) {
            plot = qpHistoBar;
            qpHistoLine->hide();
        } else {
            plot = qpHistoLine;
            qpHistoBar->hide();
        }
    } else if (text.compare("FACS") == 0) {
        LOG_MSG("save FACS plot");
        plot = qpFACS;
    } else {
        return;
    }

    QString fileName = QFileDialog::getSaveFileName(0,"Select image file", ".",
        "Image files (*.png *.jpg *.tif *.bmp)");
    if (fileName.isEmpty()) {
        return;
    }
    QSizeF size(120,120);
    QwtPlotRenderer renderer;
    renderer.renderDocument(plot,fileName,size,85);
}
*/

void MainWindow::createFACSPage()
{
    LOG_MSG("createFACSPage");
    groupBox_FACS = new QGroupBox;
    QHBoxLayout *layout1 = new QHBoxLayout;
    layout1->setMargin(1);
    qpFACS = new QwtPlot;
    layout1->addWidget(qpFACS);
    QGroupBox *vbox1 = new QGroupBox;
    QVBoxLayout *layout1v = new QVBoxLayout;
    layout1v->setMargin(1);
    vbox1->setLayout(layout1v);
    vbox1->setMinimumWidth(103);
    vbox1->setMaximumWidth(103);
    layout1v->setStretch(0,0);
    QLabel *xAxisLabel = new QLabel(" X axis");
    checkBox_FACS_log_x = new QCheckBox;
    checkBox_FACS_log_x->setText("Log scale?");
    groupBox_FACS_x_vars = new QGroupBox;
    groupBox_FACS_x_vars->setMinimumHeight(101);
    layout1v->addWidget(xAxisLabel);
    layout1v->addWidget(checkBox_FACS_log_x);
    layout1v->addWidget(groupBox_FACS_x_vars);
    layout1v->addSpacing(20);
    QLabel *yAxisLabel = new QLabel(" Y axis");
    checkBox_FACS_log_y = new QCheckBox;
    checkBox_FACS_log_y->setText("Log scale?");
    groupBox_FACS_y_vars = new QGroupBox;
    groupBox_FACS_y_vars->setMinimumHeight(101);
    layout1v->addWidget(yAxisLabel);
    layout1v->addWidget(checkBox_FACS_log_y);
    layout1v->addWidget(groupBox_FACS_y_vars);
    layout1v->addStretch(1);
    layout1->addWidget(vbox1);
    groupBox_FACS->setLayout(layout1);

    groupBox_Histo = new QGroupBox;
    QGridLayout *layout2a = new QGridLayout;
    layout2a->setMargin(1);
    qpHistoBar = new QwtPlot;
    qpHistoLine = new QwtPlot;
    layout2a->addWidget(qpHistoBar,0,0);
    layout2a->addWidget(qpHistoLine,0,0);
    QGroupBox *vBox2 = new QGroupBox;
    QVBoxLayout *layout2v = new QVBoxLayout;
    layout2v->setMargin(1);
    vBox2->setLayout(layout2v);
    vBox2->setMinimumWidth(103);
    vBox2->setMaximumWidth(103);
    QGroupBox *groupBox_celltype = new QGroupBox;
    groupBox_celltype->setMinimumHeight(80);
    groupBox_celltype->setMaximumHeight(80);
    groupBox_Histo_x_vars = new QGroupBox;
    groupBox_Histo_x_vars->setMinimumHeight(101);
    checkBox_histo_logscale = new QCheckBox;
    checkBox_histo_logscale->setText("Log scale?");
    QGroupBox *groupBox_histotype = new QGroupBox;
    groupBox_histotype->setMinimumHeight(60);
    groupBox_histotype->setMaximumHeight(60);
    layout2v->addWidget(groupBox_celltype);
    layout2v->addSpacing(20);
    layout2v->addWidget(groupBox_Histo_x_vars);
    layout2v->addStretch(1);
    layout2v->addWidget(checkBox_histo_logscale);
    layout2v->addWidget(groupBox_histotype);

    layout2a->addWidget(vBox2,0,1);
    groupBox_Histo->setLayout(layout2a);

    QRadioButton *radioButton_celltype_1 = new QRadioButton("Cell type 1");
    QRadioButton *radioButton_celltype_2 = new QRadioButton("Cell type 2");
    QRadioButton *radioButton_celltype_3 = new QRadioButton("Both types");
    buttonGroup_celltype = new QButtonGroup;
    buttonGroup_celltype->addButton(radioButton_celltype_1);
    buttonGroup_celltype->addButton(radioButton_celltype_2);
    buttonGroup_celltype->addButton(radioButton_celltype_3);
    QVBoxLayout *vbox_rb1 = new QVBoxLayout;
    vbox_rb1->addWidget(radioButton_celltype_1);
    vbox_rb1->addWidget(radioButton_celltype_2);
    vbox_rb1->addWidget(radioButton_celltype_3);
    radioButton_celltype_1->setChecked(true);
    groupBox_celltype->setLayout(vbox_rb1);

    QRadioButton *radioButton_histotype_1 = new QRadioButton("Bar plot");
    QRadioButton *radioButton_histotype_2 = new QRadioButton("Line plot");
    buttonGroup_histotype = new QButtonGroup;
    buttonGroup_histotype->addButton(radioButton_histotype_1);
    buttonGroup_histotype->setId(radioButton_histotype_1,1);
    buttonGroup_histotype->addButton(radioButton_histotype_2);
    buttonGroup_histotype->setId(radioButton_histotype_2,2);
    QVBoxLayout *vbox_rb2 = new QVBoxLayout;
    vbox_rb2->addWidget(radioButton_histotype_1);
    vbox_rb2->addWidget(radioButton_histotype_2);
    radioButton_histotype_1->setChecked(true);
    groupBox_histotype->setLayout(vbox_rb2);

    QHBoxLayout *biglayout = new QHBoxLayout;
    biglayout->addWidget(groupBox_FACS,1);
    biglayout->addWidget(groupBox_Histo,1);
    biglayout->setMargin(25);
    page_FACS->setLayout(biglayout);

    qpFACS->setTitle("FACS");
    qpFACS->replot();
    qpHistoBar->setTitle("Histogram");
    qpHistoBar->replot();
    qpHistoLine->hide();

    connect(checkBox_histo_logscale,SIGNAL(stateChanged(int)),this,SLOT(checkBox_histo_logscale_stateChanged(int)));
    connect(buttonGroup_celltype,SIGNAL(buttonClicked(QAbstractButton*)),this,SLOT(buttonGroup_celltype_buttonClicked(QAbstractButton*)));
    connect(buttonGroup_histotype,SIGNAL(buttonClicked(QAbstractButton*)),this,SLOT(buttonGroup_histotype_buttonClicked(QAbstractButton*)));
//    connect((QObject *)groupBox_FACS,SIGNAL(groupBoxClicked(QString)),this,SLOT(processGroupBoxClick(QString)));
//    connect((QObject *)groupBox_Histo,SIGNAL(groupBoxClicked(QString)),this,SLOT(processGroupBoxClick(QString)));
}


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::buttonGroup_celltype_buttonClicked(QAbstractButton* button)
{
    LOG_MSG("on_buttonGroup_celltype_buttonClicked");
    if (button->text() == "Cell type 1") {
        Global::histo_celltype = 1;
    } else if (button->text() == "Cell type 2") {
        Global::histo_celltype = 2;
    } else {
        Global::histo_celltype = 0; // both cell types
    }
    showHisto();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::checkBox_histo_logscale_stateChanged(int state)
{
    showHisto();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::buttonGroup_histotype_buttonClicked(QAbstractButton* button)
{
    showHisto();
}

//--------------------------------------------------------------------------------------------------------
// Possible variables to plot are Global::vars_used[]
// Use this to trigger FACS plot.
//--------------------------------------------------------------------------------------------------------
void MainWindow::showFACS()
{
    double xmin, xmax, ymin, ymax, x, y, xscale, yscale;
    int i, ivar, kvar_x, kvar_y, ichemox, ichemoy;
    bool x_logscale, y_logscale;
    QString xlabel, ylabel;
    QRadioButton *rb;
    QTime t;

//    if (!videoFACS->record)
//        if (!paused && !exthread->stopped) return;
    LOG_MSG("showFACS");

    qpFACS->size();
    qpFACS->setTitle("FACS");

    // Determine which x button is checked:
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        rb = FACS_x_vars_rb_list[ivar];
        if (rb->isChecked()) {
            ichemox = Global::GUI_to_DLL_index[ivar];
            break;
        }
    }
    kvar_x = ivar;
    xlabel = Global::var_string[ivar];
    xmin = Global::FACS_vmin[ivar];
    xmax = Global::FACS_vmax[ivar];
    xscale = 1;
    switch(ichemox) {
    case CFSE:
        xscale = 1000;
        xlabel = "CFSE";
        xmin = 0.1;
        xmax = 1500;
        break;
    case OXYGEN:
        xscale = 1;
        xlabel = "Oxygen";
        xmin = 1.0e-4;
        xmax = 1.0;
        break;
    case GLUCOSE:
        xscale = 1;
        xlabel = "Glucose";
        xmin = 1.0e-3;
        xmax = 10.0;
        break;
    case LACTATE:
        xscale = 1;
        xlabel = "Lactate";
        xmin = 1.0e-3;
        xmax = 10.0;
        break;
    case TRACER:
        break;
    case DRUG_A_PARENT:
//        xlabel = "Drug A";
        break;
    case DRUG_A_METAB_1:
        xscale = 1;
//        xlabel = "Drug A metabolite1";
        break;
    case DRUG_A_METAB_2:
//        xlabel = "Drug A metabolite2";
        break;
    case DRUG_B_PARENT:
//        xlabel = "Drug B";
        break;
    case DRUG_B_METAB_1:
//        xlabel = "Drug B metabolite1";
        break;
    case DRUG_B_METAB_2:
//        xlabel = "Drug B metabolite2";
        break;
    case GROWTH_RATE:
        xscale = 1;
        xlabel = "Growth rate";
        xmin = 1.0e-8;
        xmax = 2.0e-5;
        break;
    case CELL_VOLUME:
        xscale = 1;
        xlabel = "Cell volume";
        xmin = 0.5;
        xmax = 2;
        break;
    case O2_BY_VOL:
        xscale = 1;
        xlabel = "O2 x volume";
        xmin = 0.05;
        xmax = 2;
        break;
    case CYCLE_PHASE:
        xscale = 1;
        xlabel = "Cell cycle phase";
        xmin = 0;
        xmax = 8;
        break;
    }
    // Determine which y button is checked:
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        rb = FACS_y_vars_rb_list[ivar];
        if (rb->isChecked()) {
            ichemoy = Global::GUI_to_DLL_index[ivar];
            break;
        }
    }
    kvar_y = ivar;
    ylabel = Global::var_string[ivar];
    ymin = Global::FACS_vmin[ivar];
    ymax = Global::FACS_vmax[ivar];
    yscale = 1;
    switch(ichemoy) {
    case CFSE:
        yscale = 1000;
        ylabel = "CFSE";
        ymin = 0.1;
        ymax = 1500;
        break;
    case OXYGEN:
        yscale = 1;
        ylabel = "Oxygen";
        ymin = 1.0e-4;
        ymax = 1.0;
        break;
    case GLUCOSE:
        yscale = 1;
        ylabel = "Glucose";
        ymin = 1.0e-3;
        ymax = 10.0;
        break;
    case LACTATE:
        yscale = 1;
        ylabel = "Lactate";
        ymin = 1.0e-3;
        ymax = 10.0;
        LOG_MSG("Lactate selected");
        break;
    case TRACER:
        break;
    case DRUG_A_PARENT:
        break;
    case DRUG_A_METAB_1:
        break;
    case DRUG_A_METAB_2:
        break;
    case DRUG_B_PARENT:
        break;
    case DRUG_B_METAB_1:
        break;
    case DRUG_B_METAB_2:
        break;
    case GROWTH_RATE:
        yscale = 1;
        ylabel = "Growth rate";
        ymin = 1.0e-8;
        ymax = 2.0e-5;
        break;
    case CELL_VOLUME:
        yscale = 1;
        ylabel = "Cell volume";
        ymin = 0.5;
        ymax = 2;
        break;
    case O2_BY_VOL:
        yscale = 1;
        ylabel = "O2 x volume";
        ymin = 0.05;
        ymax = 2;
        break;
    case CYCLE_PHASE:
        yscale = 1;
        ylabel = "Cell cycle phase";
        ymin = 0;
        ymax = 8;
        break;
    }

    x_logscale = checkBox_FACS_log_x->isChecked();
    y_logscale = checkBox_FACS_log_y->isChecked();

    double cfse_max = 0;
    if (ichemox == CFSE || ichemoy == CFSE) {
        for (i=0; i<Global::nFACS_cells; i++) {
            if (ichemox == CFSE) {
                x = Global::FACS_data[Global::nvars_used*i+kvar_x];
                cfse_max = max(cfse_max,xscale*x);
            } else {
                y = Global::FACS_data[Global::nvars_used*i+kvar_y];
                cfse_max = max(cfse_max,yscale*y);
            }
        }
    }
    if (ichemox == CFSE) xmax = cfse_max;
    if (ichemoy == CFSE) ymax = cfse_max;

    if (xQpval) free(xQpval);
    if (yQpval) free(yQpval);
    double *xQpval = (double *)malloc((Global::nFACS_cells)*sizeof(double));
    double *yQpval = (double *)malloc((Global::nFACS_cells)*sizeof(double));
    QwtSymbol *sym = new QwtSymbol();
    QwtPlotCurve *crv = &FACS_crv;
    sym->setStyle(QwtSymbol::Ellipse);
    sym->setPen(QColor(Qt::blue));
    sym->setBrush(QColor(Qt::blue));
    int dotSize = lineEdit_dotsize->text().toInt();
    sym->setSize(dotSize);
    crv->setSymbol(sym);
    crv->setPen(QColor(Qt::red));
    crv->setStyle(QwtPlotCurve::NoCurve);
    crv->attach(qpFACS);
    QTime tt;
    tt.start();
    for (i=0; i<Global::nFACS_cells; i++) {
        x = Global::FACS_data[Global::nvars_used*i+kvar_x];
        y = Global::FACS_data[Global::nvars_used*i+kvar_y];
        x = xscale*x;
        y = yscale*y;
        y = max(y,1.01*ymin);
        if (x >= xmin) {
            xQpval[i] = x;
            yQpval[i] = y;
        }
    }

//    crv->setRawData(xQpval,yQpval,Global::nFACS_cells);
    crv->setSamples(xQpval,yQpval,Global::nFACS_cells);

    qpFACS->setAxisScale(QwtPlot::yLeft, ymin, max(1.5*ymin,ymax), 0);
    qpFACS->setAxisTitle(QwtPlot::yLeft, ylabel);
    if (y_logscale) {
        qpFACS->setAxisScaleEngine(QwtPlot::yLeft, new QwtLogScaleEngine(10));
    } else {
        qpFACS->setAxisScaleEngine(QwtPlot::yLeft, new QwtLinearScaleEngine);
    }
    qpFACS->setAxisMaxMinor(QwtPlot::yLeft, 10);
    qpFACS->setAxisMaxMajor(QwtPlot::yLeft, 5);

    qpFACS->setAxisScale(QwtPlot::xBottom, xmin, xmax, 0);
    qpFACS->setAxisTitle(QwtPlot::xBottom, xlabel);
    if (x_logscale) {
        qpFACS->setAxisScaleEngine(QwtPlot::xBottom, new QwtLogScaleEngine(10));
    } else {
        qpFACS->setAxisScaleEngine(QwtPlot::xBottom, new QwtLinearScaleEngine);
    }
    qpFACS->setAxisMaxMinor(QwtPlot::xBottom, 10);
    qpFACS->setAxisMaxMajor(QwtPlot::xBottom, 5);

    qpFACS->replot();
    free(xQpval);
    free(yQpval);
    if (videoFACS->record) {
        videoFACS->recorder();
        exthread->mutex1.unlock();
        LOG_MSG("exthread->mutex1.unlock");
    }
    LOG_QMSG("showFACS display time (ms): " + QString::number(tt.elapsed()));
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::saveFACSImage()
{
    QString path = QFileDialog::getSaveFileName(this, tr("Save as image"), "", tr("PNG file (*.png)"));
    if (path.isEmpty())
        return;

    QSizeF size(120,120);
    QwtPlotRenderer renderer;
    renderer.renderDocument(qpFACS,path,size,85);
}


//-----------------------------------------------------------------------------------------
// Get FACS data TESTING WITH IT HERE
//-----------------------------------------------------------------------------------------
void MainWindow::zzzgetFACS()
{
    get_nfacs(&Global::nFACS_cells);
    if (!Global::FACS_data || Global::nFACS_cells*Global::nvars_used > Global::nFACS_dim) {
        if (Global::FACS_data) free(Global::FACS_data);
        Global::nFACS_dim = 3*Global::nFACS_cells*Global::nvars_used;   // 3* to avoid excessive malloc/free
        Global::FACS_data = (double *)malloc(Global::nFACS_dim*sizeof(double));
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::test_histo()
{
    QString testlabel = "test";
    int numValues = 20;
    double width = 10, xmin = 0;
    QVector<double> values(numValues);
    for (int i=0; i<numValues; i++) {
        values[i] = rand() %100;
    }
    makeHistoPlot(numValues,xmin,width,values,testlabel);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: showHisto()
{
    int ivar, k, k0, numValues;
    QRadioButton *rb;
    QString xlabel;
    double width, xmin;
    bool log_scale;

    exthread->getHisto();
    log_scale = checkBox_histo_logscale->isChecked();
    numValues = Global::nhisto_bins;
    QVector<double> values(numValues);
    // Determine which button is checked:
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        rb = histo_rb_list[ivar];
        if (rb->isChecked()) {
            break;
        }
    }
    xlabel = Global::var_string[ivar];
    k0 = Global::histo_celltype*numValues*Global::nvars_used;
    if (!Global::histo_data) {
        LOG_MSG("No histo_data");
        return;
    }
    for (int i=0; i<numValues; i++) {
        k = k0 + ivar*numValues + i;
        if (log_scale)
            values[i] = Global::histo_data_log[k];
        else
            values[i] = Global::histo_data[k];
    }
    if (log_scale) {
        xmin = Global::histo_vmin_log[ivar];
        width = (Global::histo_vmax_log[ivar] - Global::histo_vmin_log[ivar])/numValues;
    } else {
        xmin = Global::histo_vmin[ivar];
        width = (Global::histo_vmax[ivar] - Global::histo_vmin[ivar])/numValues;
    }
    if (xlabel.compare("Cycle phase") == 0) {
        LOG_QMSG("xlabel: " + xlabel)
        xmin = 1;
        double xmax = 8;
        numValues = 7;
        width = (xmax - xmin)/numValues;
    }
    makeHistoPlot(numValues,xmin,width,values,xlabel);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeHistoPlot(int numValues, double xmin, double width,
                               QVector<double> values, QString xlabel)
{
    QwtPlot *plot;
    double pos;

    bool use_HistoBar = (buttonGroup_histotype->checkedId() == 1);
    if (use_HistoBar) {
        plot = qpHistoBar;
        qpHistoLine->hide();
    } else {
        plot = qpHistoLine;
        qpHistoBar->hide();
    }
    plot->setCanvasBackground(QColor(Qt::white));
    plot->setTitle("Histogram");

    QwtPlotGrid *grid = new QwtPlotGrid;
    grid->enableXMin(true);
    grid->enableYMin(true);
    grid->setMajorPen(QPen(Qt::black, 0, Qt::DotLine));
    grid->setMinorPen(QPen(Qt::gray, 0 , Qt::DotLine));
    grid->attach(plot);

    if (use_HistoBar) {
        if (histogram) {
            histogram->detach();
        } else {
            histogram = new QwtPlotHistogram();
        }

        QColor c = Qt::darkCyan;
        c.setAlpha( 180 );
        histogram->setBrush( QBrush( c ) );
        QVector<QwtIntervalSample> samples( numValues );
        pos = xmin;
        for ( uint i = 0; i < numValues; i++ )
        {
            QwtInterval interval( pos, pos+width );
            interval.setBorderFlags( QwtInterval::ExcludeMaximum );
            samples[i] = QwtIntervalSample( values[i], interval );
            pos += width;
        }
        histogram->setData( new QwtIntervalSeriesData( samples ) );
        histogram->attach(plot);
    } else {
        plot->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
        double x[100], y[100];
        for ( int i = 0; i < numValues; i++ ) {
            x[i] = xmin + (i + 0.5)*width;
            y[i] = values[i];
        }
        pos = x[numValues-1] + width/2;
        QwtPlotCurve *curve = new QwtPlotCurve("");
        QPen *pen = new QPen();
        pen->setColor(Qt::black);
        pen->setWidth(3);
        curve->attach(plot);
        curve->setPen(*pen);
        curve->setSamples(x, y, numValues);
    }

    plot->setAxisTitle(QwtPlot::xBottom, xlabel);
    plot->setAxisScale(QwtPlot::yLeft, 0.0, 100.0);
    plot->setAxisScale(QwtPlot::xBottom, xmin, pos);
    plot->replot();
    plot->show();

}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::saveHistoImage(bool save_data)
{
    QwtPlot *qplot;

    bool use_HistoBar = (buttonGroup_histotype->checkedId() == 1);
    if (use_HistoBar) {
        qplot = qpHistoBar;
    } else {
        qplot = qpHistoLine;
    }

    QString path = QFileDialog::getSaveFileName(this, tr("Save as image"), "", tr("PNG file (*.png)"));
    if (path.isEmpty())
        return;

    QSizeF size(120,120);
    QwtPlotRenderer renderer;
    renderer.renderDocument(qplot,path,size,85);

    if (save_data) {
        QString data_path = QFileDialog::getSaveFileName(this, tr("Save Histo data as text"), "", tr("TXT file (*.txt)"));
        if (data_path.isEmpty())
            return;
        saveHistoData(data_path);
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: saveHistoData(QString dataFile)
{
    int ivar, k, k0, numValues;
    QRadioButton *rb;
    QString xlabel;
    double width, xmin;
    bool log_scale;

    QFile file(dataFile);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(dataFile)
                             .arg(file.errorString()));
        LOG_MSG("File open failed");
        return;
    }
    QTextStream out(&file);

    log_scale = checkBox_histo_logscale->isChecked();
    numValues = Global::nhisto_bins;
    // Determine which button is checked:
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        rb = histo_rb_list[ivar];
        if (rb->isChecked()) {
            break;
        }
    }
    xlabel = Global::var_string[ivar];
    out << xlabel << "\n";
    k0 = Global::histo_celltype*numValues*Global::nvars_used;
    if (!Global::histo_data) {
        LOG_MSG("No histo_data");
        return;
    }
    xmin = Global::histo_vmin[ivar];
    width = (Global::histo_vmax[ivar] - Global::histo_vmin[ivar])/numValues;
    out << "width: " + QString::number(width) + "\n";
    for (int i=0; i<numValues; i++) {
        k = k0 + ivar*numValues + i;
        sprintf(msg,"%f %f",xmin+(i+0.5)*width,Global::histo_data[k]);
        out << msg << "\n";
    }
    file.close();
}
