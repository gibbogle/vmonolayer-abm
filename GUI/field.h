#ifndef FIELD_H
#define FIELD_H

#include <QDialog>
#include <QMessageBox>
#include <QtGui>
#include <QMouseEvent>
#include "plot.h"
#include "global.h"
#include "myqgraphicsview.h"

#define CANVAS_WIDTH 696
#define MAX_CONC 9 // must = MAX_CHEMO in DLL for conc[] in FIELD_DATA to be the right size
#define NEXTRA 4    // must = N_EXTRA in DLL

struct field_data {
    int site[3];
    int state;
    double volume;
    double conc[MAX_CONC+NEXTRA+1];    // added CFSE, dVdt, volume, O2byVol
};

typedef field_data FIELD_DATA;

struct new_field_data {
    int NX, NY, NZ;
    int NCONST;
    double DX;
    double *Cave;   // Cslice(NX,NY,NZ,NCONST) where NCONST = MAX_CHEMO
    int ncells;
    CELL_DATA *cell_data;
};
typedef new_field_data NEW_FIELD_DATA;


#define X_AXIS 1
#define Y_AXIS 2
#define Z_AXIS 3

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern "C" {
    void get_fieldinfo(int *, int *, double *, int *, int *, int *, int *);
    void get_fielddata(int *, double *, int *, int *, FIELD_DATA *, int *);
    void new_get_fielddata(int *, double *, new_field_data *, int *, int *);
}

class Field : public QWidget
{
public:
    Field(QWidget *);
    ~Field();
    void displayField(int, int *);
    void displayField1();
    void setSliceChanged();
    void chooseFieldColor(double c, double cmin, double cmax, bool use_log, int rgbcol[]);
    void chooseRateColor(double fr, int rgbcol[]);
    void getTitle(int iconst, QString *title);
    bool isConcPlot();
    void setConcPlot(bool);
    bool isVolPlot();
    void setVolPlot(bool);
    bool isOxyPlot();
    void setOxyPlot(bool);
    void selectCellConstituent();
    void selectFieldConstituent();
    void setExecuting(bool);
    void setSaveImages(bool);
    void setUseLogScale(bool);
    void setCellConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QList<QRadioButton *> *rb_list, QString tag);
    void setFieldConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QList<QRadioButton *> *rb_list, QString tag);

    QWidget *field_page;
    bool save_images;
    bool use_log;
    MyQGraphicsView* view;
    QGraphicsScene* scene;
    int NX;
    int axis;
    double fraction;
    int hour;
    int ifield;
    int nsites, nconst, const_used[MAX_CONC+NEXTRA+1];
    int nvars_used;
    int cvar_index[32];
    QString const_name[16];
    QString constituentText;
    int cell_constituent;
    int field_constituent;
    bool slice_changed;
    bool show_cells;
    bool useConcPlot;
    bool useVolPlot;
    bool useOxyPlot;
    FIELD_DATA *data;
    Plot *pGconc, *pGvol, *pGoxy;
    bool executing;
    char msg[1024];

    QButtonGroup *buttonGroup_cell_constituent;
    QButtonGroup *buttonGroup_field_constituent;
    QVBoxLayout *vbox_cell_constituent;
    QVBoxLayout *vbox_field_constituent;
    QList<QRadioButton *> cell_constituent_rb_list;
    QList<QRadioButton *> field_constituent_rb_list;
    QList<QLineEdit *> line_maxConc_list;
    QVBoxLayout *vbox_cell_max_concentration;

    void setPlane(QAbstractButton* button);
	void setFraction(QString text);
    void setCellConstituent(QAbstractButton* button);
    void setFieldConstituent(QAbstractButton* button);
    void setMaxConcentrations(QGroupBox *gbox);
};

#endif // FIELD_H
