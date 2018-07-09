#include "qdebug.h"
#include "field.h"
#include "log.h"

#include "global.h"

LOG_USE();

//NEW_FIELD_DATA fdata;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::Field(QWidget *aParent) : QWidget(aParent)
{
    field_page = aParent;
    axis = Z_AXIS;
    fraction = 0;
    cell_constituent = OXYGEN;
    field_constituent = OXYGEN;
    slice_changed = true;
    setConcPlot(false);
    setVolPlot(false);
    setOxyPlot(false);
    pGconc = NULL;
    pGvol = NULL;
    pGoxy = NULL;
    ifield = 0;
    view = new MyQGraphicsView(field_page);
    scene = new QGraphicsScene(QRect(0, 0, CANVAS_WIDTH, CANVAS_WIDTH));
    vbox_cell_constituent = NULL;
    vbox_field_constituent = NULL;
    vbox_cell_max_concentration = NULL;
    buttonGroup_cell_constituent = new QButtonGroup;
    buttonGroup_field_constituent = new QButtonGroup;
    line_maxConc_list.clear();
    cell_constituent_rb_list.clear();
    field_constituent_rb_list.clear();
    show_cells = true;
    data = NULL;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::~Field()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isConcPlot()
{
    return useConcPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setConcPlot(bool status)
{
    useConcPlot = status;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isVolPlot()
{
    return useVolPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setVolPlot(bool status)
{
    useVolPlot = status;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isOxyPlot()
{
    return useOxyPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setOxyPlot(bool status)
{
    useOxyPlot = status;
}

//------------------------------------------------------------------------------------------------
// To create the group of radiobuttons for cell constituent selection.
// This uses information about active constituents fetched from the DLL.
// Need to distinguish field constituent from cell constituent, because for the cell we are
// interested in CFSE, volume, O2byvol in addition to the dissolved constituents:
// oxygen, glucose, drugA, drugB, metabolites...
//------------------------------------------------------------------------------------------------
void Field::setCellConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QList<QRadioButton *> *rb_list, QString tag)
{
    int ivar;
    QString name, str, ivar_str;
    QRadioButton *rb;

//    LOG_QMSG("setCellConstituentButtons: " + tag);
    if (rb_list->length() != 0) {
        LOG_MSG("rb_list not NULL, delete it");
        for (ivar=0; ivar<rb_list->length(); ivar++) {
            rb = (*rb_list)[ivar];
            bg->removeButton(rb);
            delete rb;
        }
        rb_list->clear();
    }
    if (!*vbox) {
//        LOG_MSG("vbox = NULL, create it");
        *vbox = new QVBoxLayout;
        gbox->setLayout(*vbox);
    }
    name = "rb_cell_constituent_"+tag;
//    LOG_QMSG(name);
//    sprintf(msg,"rb_list: %p vbox: %p bg: %p nvars_used: %d",rb_list,*vbox,bg,Global::nvars_used);
//    LOG_MSG(msg);
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        ivar_str = QString::number(ivar);
        str = Global::var_string[ivar];
        rb = new QRadioButton;
        rb->setText(str);
        rb->setObjectName(name+ivar_str);
        (*vbox)->addWidget(rb);
        rb->setEnabled(true);
        bg->addButton(rb,ivar);
        rb_list->append(rb);
//        LOG_QMSG(rb->objectName());
    }
//    LOG_MSG("added buttons");
    if (tag.contains("FACS")) {
        (*rb_list)[0]->setChecked(true);   // CFSE
    } else {
        (*rb_list)[1]->setChecked(true);   // Oxygen
    }
    QRect rect = gbox->geometry();
    rect.setHeight(25*Global::nvars_used);
//    LOG_MSG("gbox->setGeometry");
    gbox->setGeometry(rect);
//    LOG_MSG("gbox->show");
    gbox->show();
//    LOG_QMSG("did setCellConstituentButtons: " + tag);
}


//------------------------------------------------------------------------------------------------
// To create the group of radiobuttons for field constituent selection.
// This uses information about active constituents fetched from the DLL.
// Need to distinguish field constituent from cell constituent, because for the
// field we have only the dissolved constituents:
// oxygen, glucose, drugA, drugB, metabolites...
//------------------------------------------------------------------------------------------------
void Field::setFieldConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QList<QRadioButton *> *rb_list, QString tag)
{
    int ivar;
    QString name, str, ivar_str;
    QRadioButton *rb;

    Global::nfieldvars_used = Global::nvars_used - Global::N_EXTRA;
//    LOG_QMSG("setFieldConstituentButtons: " + tag + " nfieldvars_used: "
//             + QString::number(Global::nfieldvars_used));
    if (rb_list->length() != 0) {
        LOG_MSG("rb_list not NULL, delete it");
        for (ivar=0; ivar<rb_list->length(); ivar++) {
            rb = (*rb_list)[ivar];
            bg->removeButton(rb);
            delete rb;
        }
        rb_list->clear();
    }
    if (!*vbox) {
//        LOG_MSG("vbox = NULL, create it");
        *vbox = new QVBoxLayout;
        gbox->setLayout(*vbox);
    }
    name = "rb_field_constituent_"+tag;
//    LOG_QMSG(name);
    for (ivar=0; ivar<Global::nfieldvars_used; ivar++) {
        ivar_str = QString::number(ivar);
        str = Global::var_string[ivar+1];
        rb = new QRadioButton;
        rb->setText(str);
        QString objname = name+ivar_str;
        rb->setObjectName(objname);
        (*vbox)->addWidget(rb);
        rb->setEnabled(true);
        bg->addButton(rb,ivar);
        rb_list->append(rb);
//        QRadioButton *chkrb = gbox->findChild<QRadioButton *>(objname);
//        if (chkrb) {
//            QString chkstr = chkrb->objectName();
//            LOG_QMSG(chkstr);
//        } else {
//            chkrb = (*rb_list)[ivar];
//            LOG_QMSG("findChild failed, but: " + chkrb->objectName());
//        }
    }
    (*rb_list)[0]->setChecked(true);   // Oxygen
    QRect rect = gbox->geometry();
    rect.setHeight(25*(Global::nfieldvars_used + 1));
    gbox->setGeometry(rect);
    gbox->show();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setMaxConcentrations(QGroupBox *gbox)
{
    int ivar;
    QLineEdit *line_maxConc;

    if (!vbox_cell_max_concentration) {
        LOG_MSG("vbox_cell_max_concentration = NULL, create it");
        vbox_cell_max_concentration = new QVBoxLayout;
        gbox->setLayout(vbox_cell_max_concentration);
    }

    if (line_maxConc_list.length() != 0) {
        for (ivar=0; ivar < line_maxConc_list.length(); ivar++) {
            line_maxConc = line_maxConc_list[ivar];
            vbox_cell_max_concentration->removeWidget(line_maxConc);
            line_maxConc->setVisible(false);
        }
    }
    line_maxConc_list.clear();

    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        line_maxConc = new QLineEdit;
        line_maxConc->setObjectName("maxConc"+ivar);
        line_maxConc->setText("1.0");
        line_maxConc->setMaximumWidth(50);
        line_maxConc->setEnabled(true);
        vbox_cell_max_concentration->addWidget(line_maxConc);
        line_maxConc_list.append(line_maxConc);
    }

    line_maxConc_list[OXYGEN]->setText("0.18");
    line_maxConc_list[GLUCOSE]->setText("5.5");
    QRect rect = gbox->geometry();
    rect.setHeight(25*Global::nvars_used);
    gbox->setGeometry(rect);
    gbox->show();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::selectCellConstituent()
{
    int iconst;
    QStringList items;

//    LOG_MSG("selectCellConstituent");
    for (iconst=0; iconst<Global::nvars_used; iconst++) {
        if (iconst == cell_constituent) continue;
        items << Global::var_string[iconst];
    }
    bool ok;
    QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                          tr("Constituent:"), items, 0, false, &ok);
    if (ok && !item.isEmpty()) {
        for (iconst=0; iconst<Global::nvars_used; iconst++) {
            if (item == Global::var_string[iconst]) {
                cell_constituent = iconst;
//                LOG_QMSG("selectCellConstituent: " + QString::number(cell_constituent));
            }
        }
    }
    slice_changed = true;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setCellConstituent(QAbstractButton *button)
{
    int res;
//    LOG_MSG("setCellConstituent");
    int prev_constituent = cell_constituent;
    QString text = button->text();
    for (int ivar=0; ivar<Global::nvars_used; ivar++) {
        if (text == Global::var_string[ivar]) {
            cell_constituent = ivar;
        }
    }

    if (cell_constituent != prev_constituent) {
//        LOG_MSG("setCellConstituent");
        slice_changed = true;
        displayField(hour,&res);
    }
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::selectFieldConstituent()
{
    int iconst;
    QStringList items;

//    LOG_MSG("selectFieldConstituent");
    for (iconst=1; iconst<Global::nfieldvars_used+1; iconst++) {
        if (iconst == field_constituent) continue;
        items << Global::var_string[iconst];
    }
    bool ok;
    QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                          tr("Constituent:"), items, 0, false, &ok);
    if (ok && !item.isEmpty()) {
        for (iconst=1; iconst<Global::nfieldvars_used+1; iconst++) {
            if (item == Global::var_string[iconst]) {
                field_constituent = iconst;
            }
        }
    }
    slice_changed = true;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setFieldConstituent(QAbstractButton *button)
{
    int res;
//    LOG_MSG("setFieldConstituent");
    int prev_constituent = field_constituent;
    QString text = button->text();
    for (int ivar=0; ivar<Global::nfieldvars_used; ivar++) {
        if (text == Global::var_string[ivar+1]) {
            field_constituent = ivar+1;
        }
    }

    if (field_constituent != prev_constituent) {
//        LOG_MSG("setFieldConstituent");
        slice_changed = true;
        displayField(hour,&res);
    }
}


//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setPlane(QAbstractButton *button)
{
    int res;
    LOG_MSG("setPlane");
    QString text = button->text();
	int prev_axis = axis;
    if (text.compare("X-Y") == 0)
        axis = Z_AXIS;
    else if (text.compare("Y-Z") == 0)
        axis = X_AXIS;
    else if (text.compare("X-Z") == 0)
        axis = Y_AXIS;
	if (axis != prev_axis) {
        slice_changed = true;
        LOG_MSG("setPlane");
        displayField(hour,&res);
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setFraction(QString text)
{
    int res;
	double prev_fraction = fraction;
	fraction = text.toDouble();
	if (fraction != prev_fraction) {
        displayField(hour,&res);
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setSliceChanged()
{
    slice_changed = true;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setSaveImages(bool save)
{
    save_images = save;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setUseLogScale(bool use_logscale)
{
    use_log = use_logscale;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseFieldColor(double c, double cmin, double cmax, bool use_logscale, int rgbcol[])
{
    double f, denom;
    int rgb_lo[3], rgb_hi[3], i;
    bool copy_scell = false;

    if (copy_scell) {
        rgbcol[0] = 0;
        if (cmax > 0) {
            rgbcol[1] = 255*min(c,cmax)/cmax;
            rgbcol[2] = 255*min(c,cmax)/cmax;
        } else {
            rgbcol[1] = 0;
            rgbcol[2] = 0;
        }
        return;
    }
    if (use_logscale) {
        if (cmin == cmax) {
            f = 1;
        } else {
//            if (cmin > 0.0001)
//                logcmin = log(cmin);
//            else
//                logcmin = 1.0e6;
//            if (c > 0.0001)
//                logc = log(c);
//            else
//                logc = 1.0e6;
            cmin = max(cmin, 0.00001);
            c = max(c, 0.00001);
            denom = (log(cmax) - log(cmin));
            if (denom < 0.001)
                f = 1;
            else
                f = (log(c) - log(cmin))/denom;
        }
    } else {
        f = c/cmax;
    }
//    if (cell_constituent == OXYGEN) {
    if (field_constituent == OXYGEN) {
        rgb_hi[0] =   0; rgb_hi[1] =   0; rgb_hi[2] = 0;
        rgb_lo[0] = 255; rgb_lo[1] =   0; rgb_lo[2] = 0;
        for (i=0; i<3; i++) {
            rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
            if (rgbcol[i] < 0 || rgbcol[i] > 255) {
                sprintf(msg,"chooseFieldColor: %f %f %f %f %d %d",c,cmin,cmax,f,i,rgbcol[i]);
                LOG_MSG(msg);
                exit(1);
            }
        }
    } else {
        rgb_hi[0] =   0; rgb_hi[1] =   255; rgb_hi[2] = 255;
        rgb_lo[0] =   0; rgb_lo[1] =   0; rgb_lo[2] = 0;
        for (i=0; i<3; i++) {
            rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
            if (rgbcol[i] < 0 || rgbcol[i] > 255) {
                sprintf(msg,"chooseFieldColor: %f %f %f %f %d %d",c,cmin,cmax,f,i,rgbcol[i]);
                LOG_MSG(msg);
                exit(1);
            }
        }
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseRateColor(double f, int rgbcol[])
{
    int rgb_lo[3], rgb_hi[3], i;

    rgb_hi[0] = 0; rgb_hi[1] = 255; rgb_hi[2] = 0;
    rgb_lo[0] = 0; rgb_lo[1] =  64; rgb_lo[2] = 0;
    for (i=0; i<3; i++) {
        rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
    }
}

//-----------------------------------------------------------------------------------------
// Now 'constituent' is an index of the active constituents: 0 - nvars_used
//-----------------------------------------------------------------------------------------
void Field::getTitle(int iconst, QString *title)
{
    QString name = Global::var_string[iconst];
    if (Global::GUI_to_DLL_index[iconst] <= Global::MAX_CHEMO) {
        *title = name + " Concentration";
    } else {
        *title = name;
    }
}


//-----------------------------------------------------------------------------------------
// New version, site/cell size is fixed, the blob grows
// Removed for monolayer
//-----------------------------------------------------------------------------------------
void Field::displayField(int hr, int *res)
{
    *res = 0;
}

