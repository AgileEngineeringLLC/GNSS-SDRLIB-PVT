//-----------------------------------------------------------------------------
// sdrsync.c : SDR sync thread functions
//
// Edits from Don Kelly, don.kelly@mac.com, 2025
//-----------------------------------------------------------------------------
#include "sdr.h"

/* synchronization thread ------------------------------------------------------
* synchronization thread for pseudo range computation
* args   : void   *arg      I   not used
* return : none
* note : this thread collects all data of sdr channel thread and compute pseudo
*        range at every output timing.
*-----------------------------------------------------------------------------*/
extern void *syncthread(void * arg)
{
    int i,j,nsat,isat[MAXOBS],ind[MAXSAT]={0},refi;
    uint64_t sampref,sampbase,codei[MAXSAT],diffcnt,mincodei;
    double codeid[OBSINTERPN],remcode[MAXSAT],samprefd,reftow=0,oldreftow;
    sdrobs_t obs[MAXSAT];
    sdrtrk_t trk[MAXSAT]={{0}};
    int ret=0; // used for function output

    // Timer for nav messages
    //time_t start_time, current_time;
    //start_time = time(NULL);

    while (!sdrstat.stopflag) {

        mlock(hobsmtx);
        /* copy all tracking data */
        for (i=nsat=0;i<sdrini.nch;i++) {
            if (sdrch[i].nav.flagdec&&sdrch[i].nav.sdreph.eph.week!=0) {
                memcpy(&trk[nsat],&sdrch[i].trk,sizeof(sdrch[i].trk));
                isat[nsat]=i;
                nsat++;
            }
        }
        unmlock(hobsmtx);

        /* find minimum tow channel (most distant satellite) */
        oldreftow=reftow;
        reftow=3600*24*7;
        for (i=0;i<nsat;i++) {
            if (trk[i].tow[0]<reftow)
                reftow=trk[i].tow[0];
        }
        /* output timing check */
        if (nsat==0||oldreftow==reftow||((int)(reftow*1000)%sdrini.outms)!=0) {
            continue;
        }
        /* select same timing index */
        for (i=0;i<nsat;i++) {
            for (j=0;j<OBSINTERPN;j++) {
                if (fabs(trk[i].tow[j]-reftow)<1E-4) {
                    ind[i]=j;
                    break;
                }
            }
            if (j==OBSINTERPN-1&&ind[i]==0)
                SDRPRINTF("error:%s reftow=%.1f tow=%.1f\n",
                    sdrch[isat[i]].satstr,trk[i].tow[OBSINTERPN-1],reftow);
        }

        /* decide reference satellite (nearest satellite) */
        mincodei=UINT64_MAX;
        refi=0;
        for (i=0;i<nsat;i++) {
            codei[i]=trk[i].codei[ind[i]];
            remcode[i]=trk[i].remcout[ind[i]];
            if (trk[i].codei[ind[i]]<mincodei) {
                refi=i;
                mincodei=trk[i].codei[ind[i]];
            }
        }
        /* reference satellite */
        diffcnt=trk[refi].cntout[ind[refi]]-sdrch[isat[refi]].nav.firstsfcnt;
        sampref=sdrch[isat[refi]].nav.firstsf+
            (uint64_t)(sdrch[isat[refi]].nsamp*
            (-PTIMING/(1000*sdrch[isat[refi]].ctime)+diffcnt));
        sampbase=trk[refi].codei[OBSINTERPN-1]-10*sdrch[isat[refi]].nsamp;
        samprefd=(double)(sampref-sampbase);

        /* computation observation data */
        for (i=0;i<nsat;i++) {
            obs[i].sys=sdrch[isat[i]].sys;
            obs[i].prn=sdrch[isat[i]].prn;
            obs[i].week=sdrch[isat[i]].nav.sdreph.week_gpst;
            obs[i].tow=reftow+(double)(PTIMING)/1000;
            obs[i].P=CLIGHT*sdrch[isat[i]].ti*
                ((double)(codei[i]-sampref)-remcode[i]); /* pseudo range */

            /* uint64 to double for interp1 */
            uint64todouble(trk[i].codei,sampbase,OBSINTERPN,codeid);
            obs[i].L=interp1(codeid,trk[i].L,OBSINTERPN,samprefd);
            obs[i].D=interp1(codeid,trk[i].D,OBSINTERPN,samprefd);
            obs[i].S=trk[i].S[0];
        }

        // Populate the obs_v vector with the initial obs data and nsat.
        // First, zero out obs_v, altho this is not really necessary, but
        // is a good precaution.
        mlock(hobsvecmtx);
        //memset(sdrstat.obs_v, 0, sizeof(sdrstat.obs_v));
        sdrstat.nsat = nsat;
        //printf("Raw obs: \n");
        //printf("nsat: %d\n", (int)sdrstat.nsat);
        for (i=0;i<nsat;i++) {
            sdrstat.obs_v[i*10+0] = obs[i].prn;  // PRN
            //sdrstat.obs_v[i*10+1] = 0.0;         // SV X posn
            //sdrstat.obs_v[i*10+2] = 0.0;         // SV Y posn
            //sdrstat.obs_v[i*10+3] = 0.0;         // SV Z posn
            sdrstat.obs_v[i*10+4] = obs[i].P;    // PR
            sdrstat.obs_v[i*10+5] = obs[i].tow;  // TOW
            sdrstat.obs_v[i*10+6] = obs[i].week; // GPS week
            sdrstat.obs_v[i*10+7] = obs[i].S;    // SNR
            //sdrstat.obs_v[i*10+8] = 0.0;         // SV Az
            //sdrstat.obs_v[i*10+9] = 0.0;         // SV El
            /*
            printf("%.0f %.1f %.1f %.1f %.1f %.1f %.0f %.1f\n",
            sdrstat.obs_v[i*8+0],
            sdrstat.obs_v[i*8+1],
            sdrstat.obs_v[i*8+2],
            sdrstat.obs_v[i*8+3],
            sdrstat.obs_v[i*8+4],
            sdrstat.obs_v[i*8+5],
            sdrstat.obs_v[i*8+6],
            sdrstat.obs_v[i*8+7]);
            */
        }
        unmlock(hobsvecmtx);

        // Perform prechecks on obs and eph data
        // Artificially mess up one PR to check precheckObs
        //sdrstat.obs_v[4] = 16000000.0;
        precheckObs();
        precheckEph();

        /*
        printf("Culled obs: \n");
        printf("nsat: %d\n", (int)sdrstat.nsat);
        for (int k=0; k<(int)sdrstat.nsat; k++) {
            printf("%.0f %.1f %.1f %.1f %.1f %.1f %.0f %.1f\n",
               sdrstat.obs_v[k*8+0],
               sdrstat.obs_v[k*8+1],
               sdrstat.obs_v[k*8+2],
               sdrstat.obs_v[k*8+3],
               sdrstat.obs_v[k*8+4],
               sdrstat.obs_v[k*8+5],
               sdrstat.obs_v[k*8+6],
               sdrstat.obs_v[k*8+7]);
        }
        */

        // Call pvtProcessor if there are at least four observations and
        // that the eph appears valid.
        mlock(hobsvecmtx);
        nsat = sdrstat.nsat;
        unmlock(hobsvecmtx);

        //printf("nsat: %d\n",nsat);
        if (nsat >= 4) {
        //if (nsat >= 5) {
            ret = pvtProcessor();
            if (ret != 0) {
              printf("errorDetected: exiting pvtProcessor\n");
            }
        }
        else {
          //printf("Less than 4 SVs, pvtProcessor bypassed (nsat: %d)\n", nsat);
        }

        // Print obs and nav to file if printflag selected (DK added)
        if (sdrstat.printflag) {
          FILE *fptr;
          char *fname = "obsData.txt";
          fptr = fopen(fname,"w");
          if (fptr == NULL) {
            perror("Error opening file\n");
          }

          for (int k=0;k<nsat;k++) {
            fprintf(fptr, "%d\n", obs[k].prn);
            fprintf(fptr, "%.14e\n", obs[k].tow);
            fprintf(fptr, "%.14e\n", obs[k].P);
          }
          for (int k=0;k<nsat;k++) {
             int ind = obs[k].prn-1;
             fprintf(fptr, "%d\n", sdrch[ind].nav.sdreph.eph.sat);
             fprintf(fptr, "%.14e\n", 0.0);
             //fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.toc);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.toes);
             fprintf(fptr, "%.14e\n", 0.0);
             //fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.ttr);
             fprintf(fptr, "%.14e\n", sqrt(sdrch[ind].nav.sdreph.eph.A));
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.e);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.M0);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.omg);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.i0);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.OMG0);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.deln);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.idot);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.OMGd);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.cuc);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.cus);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.crc);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.crs);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.cic);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.cis);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.f0);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.f1);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.f2);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.tgd[0]);
             fprintf(fptr, "%d\n", sdrch[ind].nav.sdreph.eph.week);
             fprintf(fptr, "%d\n", sdrch[ind].nav.sdreph.eph.iodc);
             fprintf(fptr, "%d\n", sdrch[ind].nav.sdreph.eph.iode);
             fprintf(fptr, "%d\n", sdrch[ind].nav.sdreph.eph.sat);
             fprintf(fptr, "%d\n", sdrch[ind].nav.sdreph.eph.svh);
             fprintf(fptr, "%d\n", sdrch[ind].nav.sdreph.eph.sva);
             fprintf(fptr, "%.14e\n", sdrch[ind].nav.sdreph.eph.fit);
          }

          SDRPRINTF("sdrData file written.\n");
          //fprintf(fptr, " Test write to file.\n");
          fclose(fptr);
          sdrstat.printflag = 0;
        }

    }

    SDRPRINTF("SDR syncthread finished!\n");

    return THRETVAL;
}
