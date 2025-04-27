//-----------------------------------------------------------------------------
// sdrmain.c : SDR main functions
//
// Copyright (C) 2014 Taro Suzuki <gnsssdrlib@gmail.com>
// Edits from Don Kelly, don.kelly@mac.com, 2025
//-----------------------------------------------------------------------------
#include "sdr.h"

// Thread handles and mutex
thread_t hmainthread;
thread_t hsyncthread;
thread_t hkeythread;
thread_t hdatathread;
thread_t hguithread;

mlock_t hbuffmtx;
mlock_t hreadmtx;
mlock_t hfftmtx;
mlock_t hobsmtx;
mlock_t hresetmtx;
mlock_t hobsvecmtx;

// SDR structs
sdrini_t sdrini={0};
sdrstat_t sdrstat={0};
sdrch_t sdrch[MAXSAT]={{0}};

// Keyboard thread ------------------------------------------------------------
// keyboard thread for program termination
// args   : void   *arg      I   not used
// return : none
//-----------------------------------------------------------------------------
extern void *keythread(void * arg)
{
    do {
        switch(getchar()) {
        case 'q':
        case 'Q':
            sdrstat.stopflag=1;
            break;
        default:
            SDRPRINTF("press 'q' to exit...\n");
            break;
        }
    } while (!sdrstat.stopflag);

    return THRETVAL;
}

// main function --------------------------------------------------------------
// main entry point in CLI application
// args   : none
// return : none
//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Set processor to Performance mode
    // Need to add this as system command

    /* read ini file */
    if (readinifile(&sdrini)<0) {
        return -1;
    }

    // Declare CPU affinity variables
    int num_cpus;
    cpu_set_t cpu_set;

    // Find the number of available processor CPUs with sysconf call.
    num_cpus = sysconf(_SC_NPROCESSORS_ONLN);

    // We will leave the last four CPUs for OS tasks, then enable the rest for
    // real-time use by main. Select CPUs to use with CPU_SET.
    for (int n=0;n<(num_cpus-4);n++) {
        CPU_SET(n, &cpu_set);
    }

    // Schedule the cores for use by main using sched_setaffinity
    //if (sched_setaffinity(0, sizeof(cpu_set_t), &cpu_set) == -1) {
    if (sched_setaffinity(getpid(), sizeof(cpu_set_t), &cpu_set) == -1) {
        perror("error: sched_setaffinity\n");
    }
    //printf(BGRN "CPU cores scheduled!\n\n" reset);

    // Start SDR and threads
    startsdr();

    return 0;
}

// sdr start ------------------------------------------------------------------
// start sdr function
// args   : void   *arg      I   not used
// return : none
// note : This is called as a function in CLI application
//-----------------------------------------------------------------------------
extern void startsdr(void)
{
    // Timer struct for terminal printf
    //struct timespec start, end, orig;

    int i;
    SDRPRINTF("GNSS-SDRLIB start!\n");

    // Establish attr and param for real time threads
    struct sched_param param1;
    struct sched_param param2;
    pthread_attr_t attr1;
    pthread_attr_t attr2;
    int ret; // thread return state

    // Initialize EKF
    ret = ekfInit();

    /*
    // Test print
    printf("xkk_v (satartsdr): \n");
    for (int i=0; i<8; i++) {
      printf("%.1f ", sdrekf.xkk_v[i]);
    }
    printf("\n\n");
    */

    // Declare thread attribute attr1 and attr2
    ret = pthread_attr_init(&attr1);
    if (ret) {
        printf("Init for thread attr1 failed: %s\n", strerror(ret));
    }
    ret = pthread_attr_init(&attr2);
    if (ret) {
        printf("Init for thread attr2 failed: %s\n", strerror(ret));
    }

    // Set stack size
    ret = pthread_attr_setstacksize(&attr1, PTHREAD_STACK_MIN);
    if (ret) {
        printf("Stack size for thread attr1 failed: %s\n", strerror(ret));
    }
    ret = pthread_attr_setstacksize(&attr2, PTHREAD_STACK_MIN);
    if (ret) {
        printf("Stack size for thread attr2 failed: %s\n", strerror(ret));
    }

    // Set scheduler policy and priority
    ret = pthread_attr_setschedpolicy(&attr1, SCHED_FIFO);
    if (ret) {
        printf("Policy for thread attr1 failed: %s\n", strerror(ret));
    }
    ret = pthread_attr_setschedpolicy(&attr2, SCHED_OTHER);
    if (ret) {
        printf("Policy for thread attr2 failed: %s\n", strerror(ret));
    }

    // Set thread priority in attr and param
    param1.sched_priority = 40;
    ret = pthread_attr_setschedparam(&attr1, &param1);
    if (ret) {
        printf("Priority for thread attr1 failed: %s\n", strerror(ret));
    }
    param2.sched_priority = 0;
    ret = pthread_attr_setschedparam(&attr2, &param2);
    if (ret) {
        printf("Priority for thread attr2 failed: %s\n", strerror(ret));
    }

    // Direct attr to use the scheduling parameters we set
    ret = pthread_attr_setinheritsched(&attr1, PTHREAD_EXPLICIT_SCHED);
    if (ret) {
        printf("Inherit for thread attr1 failed: %s\n", strerror(ret));
    }
    ret = pthread_attr_setinheritsched(&attr2, PTHREAD_EXPLICIT_SCHED);
    if (ret) {
        printf("Inherit for thread attr2 failed: %s\n", strerror(ret));
    }

    // check initial value
    if (chk_initvalue(&sdrini)<0) {
        SDRPRINTF("error: chk_initvalue\n");
        quitsdr(&sdrini,1);
        return;
    }

    // receiver initialization
    if (rcvinit(&sdrini)<0) {
        SDRPRINTF("error: rcvinit\n");
        quitsdr(&sdrini,1);
        return;
    }
    // initialize sdr channel struct
    for (i=0;i<sdrini.nch;i++) {
        if (initsdrch(i+1,sdrini.sys[i],sdrini.prn[i],sdrini.ctype[i],
            sdrini.dtype[sdrini.ftype[i]-1],sdrini.ftype[i],
            sdrini.f_gain[sdrini.ftype[i]-1],sdrini.f_bias[sdrini.ftype[i]-1],
            sdrini.f_clock[sdrini.ftype[i]-1],sdrini.f_cf[sdrini.ftype[i]-1],
            sdrini.f_sf[sdrini.ftype[i]-1],sdrini.f_if[sdrini.ftype[i]-1],
            &sdrch[i])<0) {

            SDRPRINTF("error: initsdrch\n");
            quitsdr(&sdrini,2);
            return;
        }
    }

    // mutexes and events
    openhandles();

    // Create threads ---------------------------------------------------------
    // Keyboard thread
    //ret = pthread_create(&hkeythread,&attr2,keythread,NULL);
    ret = pthread_create(&hkeythread,NULL,keythread,NULL);
    if (ret) {
        printf(BRED "Create for keyboard thread failed: %s\n" reset,
             strerror(ret));
    }

    // Sync thread
    //ret = pthread_create(&hsyncthread,&attr1,syncthread,NULL);
    ret = pthread_create(&hsyncthread,NULL,syncthread,NULL);
    if (ret) {
        printf(BRED "Create for sync thread failed: %s\n" reset,
             strerror(ret));
    }

    // GUI thread
    //ret = pthread_create(&hguithread,&attr2,guithread,NULL);
    ret = pthread_create(&hguithread,NULL,guithread,NULL);
    if (ret) {
        printf(BRED "Create for GUI thread failed: %s\n" reset,
             strerror(ret));
    }

    // SDR channel threads
    for (i=0;i<sdrini.nch;i++) {
        // GPS/QZS/GLO/GAL/CMP L1
        if (sdrch[i].ctype==CTYPE_L1CA  || sdrch[i].ctype==CTYPE_L1SBAS){
            //ret=pthread_create(&sdrch[i].hsdr,&attr1,sdrthread,&sdrch[i]);
            ret=pthread_create(&sdrch[i].hsdr,NULL,sdrthread,&sdrch[i]);
            if (ret) {
                printf(BRED "Create for sdr thread failed: %s\n" reset,
                     strerror(ret));
            } // if
        } // if
    } // for (sdrch threads)

    // Data grabber thread
    //ret=pthread_create(&hdatathread,&attr1,datathread,NULL);
    ret=pthread_create(&hdatathread,NULL,datathread,NULL);
    if (ret) {
        printf(BRED "Create for data thread failed: %s\n" reset, strerror(ret));
    }

    // Create flag that is set if a channel is reset. If this flag
    // is set, then sync thread will be reset.
    int tempFlag;

    // Get a start time
    time_t start_time_reset, current_time_reset;
    int initflag = 0;

    // Main while loop
    while (!sdrstat.stopflag) {

      // Check all channels for reset flag
      for (i=0;i<sdrini.nch;i++) {

        // Reset a particular channel if directed
        mlock(hresetmtx);
        tempFlag = sdrch[i].resetflag;
        unmlock(hresetmtx);

        if (tempFlag) {

          // Start timer for when create thread will occur
          if (!initflag) {
            start_time_reset = time(NULL);
            initflag = 1;
          }

          // Give sdrthread time to close (5s good, but maybe overkill)
          sleepms(500);

          // Free up channel memory
          freesdrch(&sdrch[i]);

          // Zero out struct parms
          memset(&sdrch[i], 0, sizeof(sdrch[i]));

          // Re-init the channel
          current_time_reset = time(NULL);
          if ( (current_time_reset - start_time_reset) > 300) {
            initsdrch(i+1,sdrini.sys[i],sdrini.prn[i],sdrini.ctype[i],
              sdrini.dtype[sdrini.ftype[i]-1],sdrini.ftype[i],
              sdrini.f_gain[sdrini.ftype[i]-1],sdrini.f_bias[sdrini.ftype[i]-1],
              sdrini.f_clock[sdrini.ftype[i]-1],sdrini.f_cf[sdrini.ftype[i]-1],
              sdrini.f_sf[sdrini.ftype[i]-1],sdrini.f_if[sdrini.ftype[i]-1],
              &sdrch[i]);

            // Re-start the channel
            cratethread(sdrch[i].hsdr, sdrthread,&sdrch[i]);

            // Reset the resetflags
            mlock(hresetmtx);
            sdrch[i].resetflag = 0;
            unmlock(hresetmtx);
          } // end if
        } // end if
      } // end for

        // Sleep before checking again
        sleepms(1000);

    } // end main while

    // Wait (pthreads join) threads
    waitthread(hsyncthread);
    for (i=0;i<sdrini.nch;i++)
        waitthread(sdrch[i].hsdr);
    waitthread(hdatathread);
    waitthread(hguithread);

    // SDR termination
    quitsdr(&sdrini,0);

    SDRPRINTF("GNSS-SDRLIB is finished!\n");
}

// sdr termination ------------------------------------------------------------
// sdr termination process
// args   : sdrini_t *ini    I   sdr init struct
// args   : int    stop      I   stop position in function 0: run all
// return : none
//-----------------------------------------------------------------------------
extern void quitsdr(sdrini_t *ini, int stop)
{
    int i;

    if (stop==1) return;

    // SDR termination
    rcvquit(ini);
    if (stop==2) return;

    // Free memory
    for (i=0;i<ini->nch;i++) freesdrch(&sdrch[i]);
    if (stop==3) return;

    // Mutexes and events
    closehandles();
    if (stop==4) return;
}

// sdr channel thread ---------------------------------------------------------
// sdr channel thread for signal acquisition and tracking
// args   : void   *arg      I   sdr channel struct
// return : none
// note : This thread handles the acquisition and tracking of one of the signals.
//        The thread is created at startsdr function.
//-----------------------------------------------------------------------------
extern void *sdrthread(void *arg)
{
    sdrch_t *sdr=(sdrch_t*)arg;
    uint64_t buffloc=0,bufflocnow=0,cnt=0,loopcnt=0;
    double *acqpower=NULL;

    int tempFlag;
    double snr;

    // Establish timer parameters
    time_t start_time_nav, start_time_snr, current_time;
    start_time_nav = time(NULL);
    start_time_snr = time(NULL);
    double elapsed_time_snr, elapsed_time_nav;

    // Slightly delay the start of each thread so they don't all start at
    // same time
    sleepms(sdr->no*500);
    //SDRPRINTF("**** %s sdr thread %d start! ****\n",sdr->satstr,sdr->no);

    mlock(hresetmtx);
    tempFlag = sdr->resetflag;
    unmlock(hresetmtx);

    //-------------------------------------------------------------------------
    // While loop for sdrch thread
    //-------------------------------------------------------------------------
    while (!sdrstat.stopflag && !tempFlag) {
    //while (!sdrstat.stopflag) {

        // SDR Channel Reset Checks -------------------------------------------
        // Calculate elapsed time since flagacq was set.
        current_time = time(NULL);
        elapsed_time_snr = current_time - start_time_snr;
        elapsed_time_nav = current_time - start_time_nav;

        // Every 30s check SNR to make sure it is not too low
        if (sdr->flagacq && elapsed_time_snr>30) {
            start_time_snr = time(NULL); // restart snr timer
            mlock(hobsmtx);
            snr = sdr->trk.S[0];
            unmlock(hobsmtx);

            // If SNR is low, set resetflag for this channel
            //if (snr<30.0) {
            if (snr<SNR_RESET_THRES) {
                mlock(hresetmtx);
                sdr->resetflag = 1;
                unmlock(hresetmtx);

                // Break out of while loop
                break;
            } // end if
        } // end if

        // Check every 90s to see if tracking and nav decode is successful.
        // Reset channel if not.
        // FIX: Need to not hard-code check for GPS week !!!!!!!!!!!!!!!!
        if (sdr->flagacq && elapsed_time_nav>90) {
            // Reset clock to check flagdec
            start_time_nav = time(NULL);

            //if (!sdr->flagtrk ||
            if (!sdr->nav.flagdec ||
                !sdr->nav.flagsync ||
                (sdr->nav.sdreph.week_gpst<GPS_WEEK) ) {

                mlock(hresetmtx);
                sdr->resetflag = 1;
                unmlock(hresetmtx);

                //SDRPRINTF(BRED);
                //SDRPRINTF("flagtrk: %d, flagdec: %d, swloop: %d, flagsync: %d, week: %d:\n",
                //          sdr->flagtrk, sdr->nav.flagdec, sdr->nav.swloop,
                //          sdr->nav.flagsync, sdr->nav.sdreph.week_gpst);
                //SDRPRINTF("Resetflag set for channel %d due to trk or nav\n",
                //          sdr->prn);
                //SDRPRINTF(reset);

                // Break out of while loop
                break;
            }  // end if
        } // end if

        // Acquisition --------------------------------------------------------
        if (!sdr->flagacq) {
            /* memory allocation */
            if (acqpower!=NULL) free(acqpower);
            acqpower=(double*)calloc(sizeof(double),sdr->nsamp*sdr->acq.nfreq);

            /* fft correlation */
            buffloc=sdraqcuisition(sdr,acqpower);

            // Start snr and timers
            start_time_snr = time(NULL);
            start_time_nav = time(NULL);
        }

        // Tracking -----------------------------------------------------------
        if (sdr->flagacq) {
            bufflocnow=sdrtracking(sdr,buffloc,cnt);
            if (sdr->flagtrk) {

                // correlation output accumulation
                cumsumcorr(&sdr->trk,sdr->nav.ocode[sdr->nav.ocodei]);

                sdr->trk.flagloopfilter=0;
                if (!sdr->nav.flagsync) {
                    pll(sdr,&sdr->trk.prm1,sdr->ctime);
                    dll(sdr,&sdr->trk.prm1,sdr->ctime);
                    sdr->trk.flagloopfilter=1;
                }
                else if (sdr->nav.swloop) {
                    pll(sdr,&sdr->trk.prm2,(double)sdr->trk.loopms/1000);
                    dll(sdr,&sdr->trk.prm2,(double)sdr->trk.loopms/1000);
                    sdr->trk.flagloopfilter=2;

                    mlock(hobsmtx);

                    // calculate observation data
                    if (loopcnt%(SNSMOOTHMS/sdr->trk.loopms)==0) {
                        setobsdata(sdr,buffloc,cnt,&sdr->trk,1);
                    } else {
                        setobsdata(sdr,buffloc,cnt,&sdr->trk,0);
                    }
                    unmlock(hobsmtx);

                    // Increment loop counter
                    loopcnt++;
                }

                if (sdr->trk.flagloopfilter) clearcumsumcorr(&sdr->trk);
                cnt++;
                buffloc+=sdr->currnsamp;

            } // end tracking if (flagtrk)
        } // end tracking if (flagacq)

        sdr->trk.buffloc=buffloc;
    } // end while

    // Thread finished
    if (sdr->flagacq) {
        SDRPRINTF("SDR channel %s thread finished! Delay=%d [ms]\n",
            sdr->satstr,(int)(bufflocnow-buffloc)/sdr->nsamp);
    } else {
        SDRPRINTF("SDR channel %s thread finished!\n",sdr->satstr);
    }

    return THRETVAL;
}

//-----------------------------------------------------------------------------
// Data grabber thread
//-----------------------------------------------------------------------------
extern void *datathread(void *arg)
{
    // start grabber
    if (rcvgrabstart(&sdrini)<0) {
        quitsdr(&sdrini,4);
        //return;  // is this needed?
    }

    // data grabber loop
    while (!sdrstat.stopflag) {
        if (rcvgrabdata(&sdrini)<0) {
            sdrstat.stopflag=ON;
            //break; // is this needed?
        }
    } // end while loop

    return THRETVAL;
}
