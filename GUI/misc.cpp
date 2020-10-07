#include <string>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#endif
#include <QTcpServer>
#include <QTcpSocket>
#include <QtGui>
#include <QTcpServer>
#include <QMessageBox>

#include "misc.h"
#include "log.h"
#include "transfer.h"

#include "libvmonolayer.h"

#include "global.h"

LOG_USE();
char msg[2048];

class SleeperThread : public QThread
{
public:
    static void msleep(unsigned long msecs)
    {
        QThread::msleep(msecs);
    }
};

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
SocketHandler::SocketHandler(int newport, QObject *parent)
	: QThread(parent)
{
    exiting = false;
    port = newport;
	sprintf(msg,"SocketHandler: port: %d",port);
	LOG_MSG(msg);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
SocketHandler::~SocketHandler() // make sure the worker object is destroyed
{
    exiting = true;
    wait();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::stop()
{
	LOG_MSG("SocketHandler::stop: set stopped");
	stopped = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::run()
{
//	QObject::moveToThread(this);
	sprintf(msg,"run: port: %d", port);
	LOG_MSG(msg);
	quint16 qport = port;
	QString addressStr = "127.0.0.1";
	QHostAddress hostAddress;
	hostAddress.setAddress(addressStr);
    tcpServer = new QTcpServer(this);
	stopped = false;
	connect(tcpServer, SIGNAL(newConnection()), this, SLOT(processor()), Qt::DirectConnection);
    if (!tcpServer->listen(hostAddress,qport)) {
 //       QMessageBox::critical(this, tr("Fortune Server"),
 //                              tr("Unable to start the server: %1.")
 //                              .arg(tcpServer->errorString()));
		sprintf(msg,"Unable to start the server: port: %d", port);
		LOG_MSG(msg);
        return;
    }
	sprintf(msg,"Listening on port: %d",tcpServer->serverPort());
	LOG_MSG(msg);
	LOG_MSG("serverAddress:");
	LOG_QMSG((tcpServer->serverAddress()).toString());
	bool timedOut = false;
	tcpServer->waitForNewConnection(-1,&timedOut);
	exec();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::processor()
{
//	LOG_MSG("In processor");
    socket = tcpServer->nextPendingConnection();
	sprintf(msg,"got server connection: %p",socket);	
	LOG_MSG(msg);
    emit sh_connected();
	QString qdata;
	QByteArray ba;
	ba.resize(1024);
	while (true) {
		if (stopped) {
			LOG_MSG("SocketHandler::processor: stopped!");
			break;
		}
		socket->waitForReadyRead(100);
		int nb = socket->bytesAvailable();
		if (nb > 0) {
			ba = socket->readLine(1024);
			qdata = QString(ba);
			QStringList s = qdata.split("^",QString::SkipEmptyParts);
			for (int k=0; k<s.length(); k++) {
				emit sh_output(s[k]); // Emit signal to update GUI
				if (port == CPORT0) {
					LOG_QMSG(s[k]);
				}
			}
			if (quitMessage(qdata)) {
				sprintf(msg,"Closing connection: port: %d", port);
				LOG_MSG(msg);
		        break;
			} else {
//				LOG_MSG("No bytes yet");
			}
		}
	}
	socket->close();
    LOG_MSG("closed socket");
	tcpServer->close();
    LOG_MSG("closed tcpServer");
	if (port == CPORT0) {
		emit sh_disconnected();		// Is it right that both threads emit this?
    }
    LOG_MSG("processor ends");
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
ExecThread::ExecThread(QString infile)
{
	inputFile = infile;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::run()
{
	LOG_MSG("Invoking DLL...");
	int res=0;
    double hour;
//	const char *infile, *outfile;
    char *infile, *outfile;
    char version[12];
	QString infile_path, outfile_path;
	int len_infile, len_outfile;
    int len_version;
    QString dll_version;
    bool cused[32];
    bool USING_PI = false;

    get_dll_build_version(version,&len_version);
    dll_version = version;
    // We can compare the version of the actual DLL that is linked
    // with the GUI's idea of the DLL version (i.e. the version of the library the GUI was built with)
    if (dll_version != Global::DLL_build_version) {
        LOG_QMSG("Bad DLL version: " + dll_version);
        emit(badDLL(dll_version));
        return;
    }
    paused = false;

    infile_path = inputFile;
	QString casename = QFileInfo(inputFile).baseName();
    outfile_path = casename.append(".res");
    len_outfile = outfile_path.length();
    len_infile = infile_path.length();

//	std::string std_infile = infile_path.toStdString();
//	infile = std_infile.c_str();
//	std::string std_outfile = outfile_path.toStdString();
//    outfile = std_outfile.c_str();
//    execute(&ncpu,const_cast<char *>(infile),&len_infile,const_cast<char *>(outfile),&len_outfile,&res);

    QByteArray ba = infile_path.toLocal8Bit();
    infile = ba.data();
    ba = outfile_path.toLocal8Bit();
    outfile = ba.data();

    LOG_MSG("call execute");
    execute(&ncpu,infile,&len_infile,outfile,&len_outfile,&res);
    LOG_MSG("did execute");
    if (res) {
        terminate_run(&res);
        return;
    }

//    char *b;
//    get_string(&b);
//    LOG_MSG(b);

    get_dimensions(&nsteps, &Global::DELTA_T, &Global::NT_DISPLAY, &Global::MAX_CHEMO, &Global::N_EXTRA, cused);
//    summary_interval = int(3600./Global::DELTA_T);
    summary_interval = Global::NT_DISPLAY;
    sprintf(msg,"exthread: nsteps: %d summary_interval: %d",nsteps,summary_interval);
    LOG_MSG(msg);
    Global::conc_nc_ic = 0;
    Global::conc_nc_ex = 0;
    hour = 0;

    int ndisplay = nsteps/Global::NT_DISPLAY + 1;
    mutex1.lock();
    emit setupC(ndisplay);
    mutex1.unlock();

//    LOG_MSG("call tester");
//    tester();
//    LOG_MSG("did tester");
//    emit run_tester();

    get_summary(Global::summaryData, &Global::i_hypoxia_cutoff, &Global::i_growth_cutoff);
//    getProfiles();
    sprintf(msg,"did get_summary: start hour: %f glucose: %f",hour,Global::summaryData[27]);
    LOG_MSG(msg);
    mutex1.lock();
    LOG_MSG("locked: emit summary\n");
    emit summary(hour);		// Emit signal to initialise summary plots

    for (int i=1; i <= nsteps+1; i++) {
        sleep(10);
        bool updated = false;
		if (paused && !updated) {
//			snapshot();
//            sprintf(msg,"got snapshot: i: %d",i);
//            LOG_MSG(msg);
            updated = true;
		}
        while(paused || Global::leftb) {
            sleep(100);
		}
        if (stopped) {
            LOG_MSG("stopped");
            res = -1;
            break;
        }

        simulate_step(&res);
//        LOG_MSG("did simulate_step");
        if (res != 0) {
            LOG_MSG("simulate_step: res != 0");
            break;
        }

        if (i%summary_interval == 0) {
            get_summary(Global::summaryData, &Global::i_hypoxia_cutoff, &Global::i_growth_cutoff);
//            getProfiles();
//            get_volprob(&Global::vol_nv, &Global::vol_v0, &Global::vol_dv, Global::volProb);
//            get_oxyprob(&Global::oxy_nv, &Global::oxy_v0, &Global::oxy_dv, Global::oxyProb);
//            get_distdata(&Global::dist_nv, Global::distParams, Global::distData);
            get_concdata(&Global::conc_nvars, &Global::conc_nc_ex, &Global::conc_dx_ex, Global::concData);
            if (USING_PI) get_pi_dist(PI_NBINS, Global::PI_fract, &Global::PI_max_fract, &Global::PI_max_fluor);
//            get_ic_concdata(&Global::conc_nvars, &Global::conc_nc_ic, &Global::conc_dx_ic, Global::IC_concData);

// The problem is that getFACS takes a long time when the simulation is running and the number of cells is large.  FIXED
            if (Global::showingFACS || Global::recordingFACS) {
                getFACS();
            }
            if (Global::showingFACS || Global::recordingFACS) {
                LOG_MSG("emit facs_update");
                mutex1.lock();
                emit facs_update();
                mutex1.lock();
                emit histo_update();
                mutex1.unlock();
            }
            hour = hour + (Global::DELTA_T*Global::NT_DISPLAY)/3600.;
            mutex1.lock();
            emit summary(hour);		// Emit signal to update summary plots, at hourly intervals
            mutex1.lock();
            mutex1.unlock();
        }

        if (stopped) {
            res = -1;
            break;
        }
//        if (i%Global::nt_vtk == 0) {
//            if (Global::showingVTK || Global::recordingVTK) {
//				snapshot();
//                Global::istep = i;
//                sleep(10);
//			}
//		}
        if (stopped) {
            res = -1;
            break;
        }
        if (i == nsteps+1) {
//            mutex1.lock();
            stopped = true;
//            mutex1.unlock();
            break;
        }
    }
    LOG_MSG("ExecThread::run: stopped or completed");
//    snapshot();
//    LOG_MSG("got snapshot:");
    sleep(100);

    if (Global::simulate_colony) {
        make_colony_distribution(&Global::colony_days, Global::dist, &Global::ddist, &Global::ndist);
    }

	LOG_MSG("ExecThread::run: call terminate_run");
	terminate_run(&res);

	return;
}

/*
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::wait_to_go()
{
    for (;;) {
        sprintf(msg,"waiting");
        if (goflag || stopped) break;
    }
}
*/

//-----------------------------------------------------------------------------------------
// Note that storage for BC_list, DC_list, bond_list is provided in the GUI code
// (see mainwindow.cpp and transfer.h).  The DLL fills in the data, and the number of elements
// in the lists is returned in nBC_list, nDC_list, nbond_list.
//-----------------------------------------------------------------------------------------
void ExecThread::snapshot()
{
    return;
    mutex2.lock();
//    get_scene(&Global::ncell_list,Global::cell_list);
    if (Global::ncell_list > MAX_CELLS) {
        LOG_MSG("Error: MAX_CELLS exceeded");
        exit(1);
    }
    mutex2.unlock();
    emit display(); // Emit signal to update VTK display
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::getProfiles()
{
    int k;

//    k = PROFILE_OXYGEN;
//    get_profile_oxygen(Global::profile_x[k],Global::profile_y[k],&Global::profile_n[k]);

//    k = PROFILE_S1PR1;
//    get_profile_s1pr1(profile_x[k],profile_y[k],&profile_n[k]);
//    k = PROFILE_CFSE;
//    get_profile_cfse(profile_x[k],profile_y[k],&profile_n[k]);
//    k = PROFILE_STIM;
//    get_profile_stim(profile_x[k],profile_y[k],&profile_n[k]);
//    k = PROFILE_STIMRATE;
//    get_profile_stimrate(profile_x[k],profile_y[k],&profile_n[k]);
//    k = PROFILE_AVIDITY_LN;
//    get_profile_avidity_ln(profile_x[k],profile_y[k],&profile_n[k]);
//    k = PROFILE_AVIDITY_PER;
//    get_profile_avidity_per(profile_x[k],profile_y[k],&profile_n[k]);
//    k = PROFILE_GENERATION_LN;
//    get_profile_generation_ln(profile_x[k],profile_y[k],&profile_n[k]);
//    k = PROFILE_FIRSTDCCONTACTTIME;
//    get_profile_firstdccontacttime(profile_x[k],profile_y[k],&profile_n[k]);
//    k = PROFILE_DCBINDTIME;
//    get_profile_dcbindtime(profile_x[k],profile_y[k],&profile_n[k]);
}


//-----------------------------------------------------------------------------------------
// Get FACS data
//-----------------------------------------------------------------------------------------
void ExecThread::getFACS()
{
    LOG_MSG("get_nfacs");
    get_nfacs(&Global::nFACS_cells);
    sprintf(msg,"nFACS_cells: %d",Global::nFACS_cells);
    LOG_MSG(msg);
    if (!Global::FACS_data || Global::nFACS_cells*Global::nvars_used > Global::nFACS_dim) {
        if (Global::FACS_data) {
            LOG_MSG("Free FACS_data");
            free(Global::FACS_data);
        }
        Global::nFACS_dim = Global::nFACS_cells*Global::nvars_used;   // was 3* to avoid excessive malloc/free
        Global::FACS_data = (double *)malloc(Global::nFACS_dim*sizeof(double));
        sprintf(msg,"Allocated FACS_data: %d",Global::nFACS_dim);
        LOG_MSG(msg);
    }
    LOG_MSG("get_facs");
    QTime t;
    t.start();
    get_facs(Global::FACS_data, Global::FACS_vmin, Global::FACS_vmax, Global::FACS_vmin_log, Global::FACS_vmax_log, Global::volume_scaling);
    LOG_MSG("did get_facs");
    LOG_QMSG("get_facs time (ms): " + QString::number(t.elapsed()));

}

//-----------------------------------------------------------------------------------------
// Get histogram data
//-----------------------------------------------------------------------------------------
void ExecThread::getHisto()
{
    LOG_MSG("getHisto");
    if (!Global::histo_data || Global::nhisto_bins*Global::nvars_used > Global::nhisto_dim) {
        if (Global::histo_data) free(Global::histo_data);
        if (Global::histo_data_log) free(Global::histo_data_log);
        Global::nhisto_dim = 6*Global::nhisto_bins*Global::nvars_used;   // 2*3 to avoid excessive malloc/free (only 3* used)
        Global::histo_data = (double *)malloc(Global::nhisto_dim*sizeof(double));
        Global::histo_data_log = (double *)malloc(Global::nhisto_dim*sizeof(double));
        sprintf(msg,"Allocated histo_data: %d",Global::nhisto_dim);
        LOG_MSG(msg);
    }
    LOG_MSG("get_histo");
    get_histo(Global::nhisto_bins, Global::histo_data, Global::histo_vmin, Global::histo_vmax,
              Global::histo_data_log, Global::histo_vmin_log, Global::histo_vmax_log, Global::volume_scaling);
    LOG_MSG("did get_histo");
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::stop()
{
	stopped = true;
	LOG_MSG("ExecThread::stop: stopped");
}

	//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::pause()
{
	paused = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::unpause()
{
	paused = false;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool quitMessage(QString msg)
{
	if (msg.contains("__EXIT__",Qt::CaseSensitive))
		return true;
	else
		return false;
}
