#include "mainwindow.h"
#include "QMessageBox"
#include "QFile"
#include <QDebug>

void MainWindow::SetupProtocol()
{
    QStringList tableHeader;

    tableWidget->setRowCount(64);
    tableWidget->setColumnCount(12);
    tableHeader<<"Hour"<<"Drug"<<"Duration"<<"Volume"<<"O2 conc"<<"O2 flush"<<"Conc"<<"Radiation"<<"Medium vol"<<"Full? (Y/N)"<<"Medium O2"<<"Medium glu";
    tableWidget->setHorizontalHeaderLabels(tableHeader);

    connect(tableWidget, SIGNAL(cellChanged(int, int)  ),this, SLOT(ProtocolChanged(int, int)) );

}

void MainWindow::LoadProtocol(QString fileName)
{
    int nTimes;
    QTableWidget *table = tableWidget;
    table->clearContents();
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }
    QTextStream in(&file);
    QString line;
    // Skip lines
    for (;;) {
        line = in.readLine();
        if (line.compare("PROTOCOL")==0) break;
        if (in.atEnd()) {
            QMessageBox::warning(this, "LoadProtocol", "No protocol lines in the input file");
            return;
        }
    }
    line = in.readLine();
    nTimes = line.toInt();
    for (int row=0; row<nTimes; row++) {
        QString mode = in.readLine();
        if (mode.compare("DRUG") == 0) {
            QString drug = in.readLine();
            setField(table, row, 1, drug);
            QString hour = in.readLine();
            setField(table, row, 0, hour);
            QString duration = in.readLine();
            setField(table, row, 2, duration);
            QString volume = in.readLine();
            setField(table, row, 3, volume);
            QString O2conc = in.readLine();
            setField(table, row, 4, O2conc);
            QString O2flush = in.readLine();
            setField(table, row, 5, O2flush);
            QString conc = in.readLine();
            setField(table, row, 6, conc);
        } else if (mode.compare("RADIATION") == 0) {
            QString hour = in.readLine();
            setField(table, row, 0, hour);
            QString dose = in.readLine();
            setField(table, row, 7, dose);
        } else if (mode.compare("MEDIUM") == 0) {
            QString hour = in.readLine();
            setField(table, row, 0, hour);
            QString volume = in.readLine();
            setField(table, row, 8, volume);
            QString full = in.readLine();
            setField(table, row, 9, full);
            QString O2level = in.readLine();
            setField(table, row, 10, O2level);
            QString glulevel = in.readLine();
            setField(table, row, 11, glulevel);
        }
    }
    paramSaved = false;
}

void MainWindow::SaveProtocol(QTextStream *out, int ndrugs)
{
    int nTimes, eventType, kevents;
    int err;
    QString entry, hour;
    QString drugEntry, radiationEntry, mediumEntry;
    QTableWidgetItem *item;
    QTableWidget *table = tableWidget;
    QMessageBox msgBox;
    int idrug=1, iradiation=2,imedium=3;

    int row = tableWidget->currentRow();
    int col = tableWidget->currentColumn();
    tableWidget->setCurrentCell(row+1,col);
    *out << "PROTOCOL\n";
    nTimes = 0;
    if (ndrugs == 0) {
        for (int row=0; row<tableWidget->rowCount(); row++) {
            item = tableWidget->item(row,7);
            if (item != 0) {
                entry = item->text();
                if (entry.compare("")) {    // true if <>
                    nTimes++;
                }
            }
            item = tableWidget->item(row,8);
            if (item != 0) {
                entry = item->text();
                if (entry.compare("")) {    // true if <>
                    nTimes++;
                }
            }
        }
        if (nTimes == 0) {      // nothing in the protocol
            *out << 0 << "\n";
            return;
        }
    }
    nTimes = 0;
    for (int row=0; row<tableWidget->rowCount(); row++) {
        item = tableWidget->item(row,0);
        if (item != 0) {
            entry = item->text();
            if (entry.compare("")) {    // true if <>
                nTimes++;
            }
        }
    }
    *out << nTimes << "\n";
    for (int row=0; row<tableWidget->rowCount(); row++) {
        item = tableWidget->item(row,0);
        if (item != 0) {
            entry = item->text();
            if (entry.compare("")) {    // true if <>
                hour = entry;
                kevents = 0;
                err = getField(table,row,1,&drugEntry);
                if (drugEntry.compare("")) {   // Entry in DRUG column
                    eventType = idrug;
                    kevents++;
                }
                err = getField(table,row,7,&radiationEntry);
                if (radiationEntry.compare("")) {   // Entry in RADIATION column
                    eventType = iradiation;
                    kevents++;
                }
                err = getField(table,row,8,&mediumEntry);
                if (mediumEntry.compare("")) {   // Entry in MEDIUM column
                    eventType = imedium;
                    kevents++;
                }
                if (kevents == 0) {
                    QString msg = "No DRUG, RADIATION or MEDIUM data for event at hour: " + hour;
                    msgBox.setText(msg);
                    msgBox.exec();
                    return;
                }
                if (kevents > 1) {
                    QString msg = "More than one event at hour: " + hour;
                    msgBox.setText(msg);
                    msgBox.exec();
                    return;
                }
                if (eventType == idrug) {
                    *out << "DRUG" << "\n";
                    *out << drugEntry << "\n";
                    *out << hour << "\n";
                    err = getField(table,row,2,&entry);
                    if (entry.compare("")) {   // Entry in Duration column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in Duration column");
                        msgBox.exec();
                        qDebug() << "Missing entry in Duration column";
                        return;
                    }
                    err = getField(table,row,3,&entry);
                    if (entry.compare("")) {   // Entry in Volume column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in Volume column");
                        msgBox.exec();
                        qDebug() << "Missing entry in Volume column";
                        return;
                    }
                    err = getField(table,row,4,&entry);
                    if (entry.compare("")) {   // Entry in O2conc column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in O2conc column");
                        msgBox.exec();
                        qDebug() << "Missing entry in O2conc column";
                        return;
                    }
                    err = getField(table,row,5,&entry);
                    if (entry.compare("")) {   // Entry in O2flush column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in O2flush column");
                        msgBox.exec();
                        qDebug() << "Missing entry in O2flush column";
                        return;
                    }
                    err = getField(table,row,6,&entry);
                    if (entry.compare("")) {   // Entry in Conc column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in Conc column");
                        msgBox.exec();
                        qDebug() << "Missing entry in Conc column";
                        return;
                    }
                    continue;
                }
                else if (eventType == iradiation) {
                    *out << "RADIATION" << "\n";
                    *out << hour << "\n";
                    *out << radiationEntry << "\n"; // dose
                    continue;
                }
                else if (eventType == imedium) {
                    *out << "MEDIUM" << "\n";
                    *out << hour << "\n";
                    *out << mediumEntry << "\n";    // volume
                    err = getField(table,row,9,&entry);
                    if (!entry.compare("")) entry = "N";
                    *out << entry << "\n";          // full
                    err = getField(table,row,10,&entry);
                    if (entry.compare("")) {   // Entry in Medium O2conc column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in Medium O2conc column");
                        msgBox.exec();
                        qDebug() << "Missing entry in Medium O2conc column";
                        return;
                    }
                    err = getField(table,row,11,&entry);
                    if (entry.compare("")) {   // Entry in Medium gluconc column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in Medium gluconc column");
                        msgBox.exec();
                        qDebug() << "Missing entry in Medium gluconc column";
                        return;
                    }
                    continue;
                }
            }
        }
    }
}

bool MainWindow::ProtocolUsesDrug()
{
    QTableWidgetItem *item;

    for (int row=0; row<tableWidget->rowCount(); row++) {
        item = tableWidget->item(row,1);
        if (item != 0) {
            if (item->text().compare("")) {    // true if <>
                return true;
            }
        }
    }
    return false;
}

int MainWindow::getField(QTableWidget *table, int row, int col, QString *entry)
{
    int err=0;
    QTableWidgetItem *item;
    item = table->item(row,col);
    if (item != 0) {
        *entry = item->text();
    } else {
        *entry = "";
    }
    return err;
}

void MainWindow::setField(QTableWidget *table, int row, int col, QString entry)
{
    QTableWidgetItem *item = new QTableWidgetItem;
    item->setText(entry);
    table->setItem(row,col,item);
}

void MainWindow::ProtocolChanged(int row, int col)
{
    LOG_MSG("ProtocolChanged");
//    QMessageBox::information(this,"Changed","Changed");
    paramSaved = false;
    if (col == 1) {
        QTableWidgetItem *item;
        item = tableWidget->item(row,col);
        QString drugname = item->text();
        LOG_QMSG(drugname);
        if (drugname == "") return;
        if (drugname != text_DRUG_A_NAME->text() && drugname != text_DRUG_B_NAME->text()) {
            LOG_QMSG("drug " + drugname + " is not selected");
            QMessageBox::warning(this, tr("Protocol"),
                                 tr("This drug has not been selected: %1\n %2")
                                 .arg(drugname)
                                 .arg("Select on the Treatment screen, then reload the input file."));
            item->setText("");
        }
    }
}
