//-----------------------------------------------------------------------------
//
// sdrpvt.c.c
//
// Function for computing PVT from a set of observations
// Copyright 2025, Don Kelly, don.kelly@mac.com
//-----------------------------------------------------------------------------

// Include libraries
#include "nml.h"
#include "sdr.h"

//-----------------------------------------------------------------------------
// PVT functions
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function that drives PVT calculations
//-----------------------------------------------------------------------------
extern int pvtProcessor(void) {

  mlock(hobsvecmtx);
  int numSat = sdrstat.nsat;
  unmlock(hobsvecmtx);

  double *pr_v = (double *)malloc(numSat * sizeof(double));
  double *Xs_v = (double *)malloc(numSat * 3 * sizeof(double));
  //double *G_v = (double *)malloc(numSat * 3 * sizeof(double));

  // Create dynamic variables to store obs input data
  int *satList_v = (int *)malloc(numSat * sizeof(int));
  double *rcvr_tow_v = (double *)malloc(numSat * sizeof(double));
  double *prRaw_v = (double *)malloc(numSat * sizeof(double));
  double *prSvClkCorr_v = (double *)malloc(numSat * sizeof(double));
  double *snr_v = (double *)malloc(numSat * sizeof(double));

  double tau;
  double lambda;
  double phi;
  double height;
  double rcvr_tow;
  double lat, lon;

  // Initialize obs input data
  mlock(hobsvecmtx);
  for (int i=0; i<numSat; i++) {
    satList_v[i]       = (int)sdrstat.obs_v[i*10];
    sdrstat.satList[i] = (int)sdrstat.obs_v[i*10];
    prRaw_v[i]          = sdrstat.obs_v[i*10+4];
    rcvr_tow_v[i]      = sdrstat.obs_v[i*10+5];
    snr_v[i]           = sdrstat.obs_v[i*10+7];
  }
  unmlock(hobsvecmtx);

  // Set rcvr time
  rcvr_tow = rcvr_tow_v[0];

  // Initialize parameters
  double xyzdt_v[] = {0,0,0,0};
  double xs_v[3];
  double svClkCorr;
  double transmitTime;
  double gdop = 0.0;
  int ret = 0;

  //---------------------------------------------------------------------------
  // Calculate SV positions
  //---------------------------------------------------------------------------
  for (int i=0; i<numSat; i++) {
    // Correct PR for user clock bias

    mlock(hobsvecmtx);
    pr_v[i] = prRaw_v[i] - sdrstat.xyzdt[3];
    unmlock(hobsvecmtx);

    // Calculate transit time
    tau = pr_v[i] / CTIME;

    // Update transmit time and get satellite position
    transmitTime = rcvr_tow - tau;
    ret = satPos(&sdrch[satList_v[i]-1].nav.sdreph, transmitTime, xs_v,
             &svClkCorr);
    if (ret != 0) {
      printf("Function satPos has xs NaN for G%02d, exiting pvtProcessor\n",
             satList_v[i]);
      goto errorDetected;
    }

    // If xs calculated as NaN, exit pvtProcessor
    if (isnan(xs_v[0]) || isnan(xs_v[0]) || isnan(xs_v[0])) {
      //printf("Function satPos has xs NaN for i=%02d, exiting pvtProcessor\n",
      //       i);
      goto errorDetected;
    }

    // Correct raw PR for satellite clock bias
    prSvClkCorr_v[i] = prRaw_v[i] + (CTIME * svClkCorr);

    // Build Xs_v (vector of sat positions) from each individual sat
    // posn xs_v
    Xs_v[(i*3)+0] = xs_v[0];
    Xs_v[(i*3)+1] = xs_v[1];
    Xs_v[(i*3)+2] = xs_v[2];

    // Load SV positions into obs_v
    mlock(hobsvecmtx);
    sdrstat.obs_v[i*10+1] = xs_v[0];
    sdrstat.obs_v[i*10+2] = xs_v[1];
    sdrstat.obs_v[i*10+3] = xs_v[2];
    unmlock(hobsvecmtx);

  } // end for numSat

  /*
  double pos_v[] = {0,0,0,0};
  mlock(hobsvecmtx);
  if (sdrekf.startEKF) {
    ret = bancroft(sdrstat.obs_v, pos_v, numSat);
    //printf("Bancroft initialized xkk for EKF\n");
    printf("numSat: %d, pos_v: %.1f, %.1f, %.1f, %.1f\n",
            numSat, pos_v[0],pos_v[1],pos_v[2],pos_v[3]);
    //sdrekf.xkk_v[0] = pos_v[0];
    //sdrekf.xkk_v[1] = pos_v[1];
    //sdrekf.xkk_v[2] = pos_v[2];
    //sdrekf.xkk_v[6] = pos_v[3];
    //sdrekf.startEKF = 0;
  }
  unmlock(hobsvecmtx);
  */

  // Check SV elevations if we have a receiver position solution to
  // start with
  if (sdrstat.pvtflag) {
    checkSvEl();
  }

  // Continue if we still have at least 4 SVs
  if (sdrstat.nsat<4) {
    goto errorDetected;
  }

  // Run least squares to calculate new receiver position and receiver
  // clock bias
  if (sdrini.ekfFilterOn==0) {
    ret = blsFilter(Xs_v, prSvClkCorr_v, numSat, xyzdt_v, &gdop);
  } else {
    ret = ekfFilter(Xs_v, prSvClkCorr_v, numSat, xyzdt_v, &gdop);
  }

  // If x calculated as NaN, exit pvtProcessor
  if (isnan(xyzdt_v[0]) || isnan(xyzdt_v[1]) || isnan(xyzdt_v[2])) {
    printf("Function estRcvrPosn gets NaN for xu, exiting pvtProcessor\n");
    goto errorDetected;
  }

  // Convert ECEF to LLA
  ecef2lla(xyzdt_v[0], xyzdt_v[1], xyzdt_v[2], &lambda, &phi, &height);

  // Convert radians to degrees
  lat = phi * 180.0/M_PI;
  lon = lambda * 180.0/M_PI;

  // Save LLA to sdrstat struct
  mlock(hobsvecmtx);
  sdrstat.lat = lat;
  sdrstat.lon = lon;
  sdrstat.hgt = height;
  sdrstat.gdop = gdop;
  sdrstat.xyzdt[0] = xyzdt_v[0];
  sdrstat.xyzdt[1] = xyzdt_v[1];
  sdrstat.xyzdt[2] = xyzdt_v[2];
  sdrstat.xyzdt[3] = xyzdt_v[3];
  unmlock(hobsvecmtx);

  // Free memory
  free(Xs_v);
  free(pr_v);

  free(satList_v);
  free(prRaw_v);
  free(prSvClkCorr_v);
  free(snr_v);
  free(rcvr_tow_v);

  return 0;

errorDetected:
  // Save LLA to sdrstat struct
  mlock(hobsvecmtx);
  sdrstat.lat = 0.0;
  sdrstat.lon = 0.0;
  sdrstat.hgt = 0.0;
  sdrstat.gdop = 0.0;
  unmlock(hobsvecmtx);

  // Free memory
  free(Xs_v);
  free(pr_v);

  free(satList_v);
  free(prRaw_v);
  free(prSvClkCorr_v);
  free(snr_v);
  free(rcvr_tow_v);

  sdrstat.pvtflag = 0; // Not good solution
  return -1;
} // end function

//-----------------------------------------------------------------------------
// Estimate receiver position function with BLS filter
//-----------------------------------------------------------------------------
extern int blsFilter(double *X_v, double *pr_v, int numSat,
                 double xyzdt_v[], double *gdop) {

  // Declare matrix types
  nml_mat *x, *A, *pos, *omc;
  nml_mat *tempMat1, *tempMat2, *tempMat3, *tempMat5, *tempMat6, *tempPos;
  nml_mat_lup *tempMat4;

  // Set up empty matrices
  x = nml_mat_new(4,1);
  pos = nml_mat_new(4,1);
  A = nml_mat_new(numSat,4);
  omc = nml_mat_new(numSat,1);
  tempMat2 = nml_mat_new(4,numSat);
  tempMat3 = nml_mat_new(4,4);
  tempMat5 = nml_mat_new(4,4);
  tempMat6 = nml_mat_new(4,numSat);
  tempPos = nml_mat_new(4,1);

  // Dynamic and static variables
  double Rot_X_v[] = {0,0,0};
  double pos_v[] = {0,0,0,0};
  double rho2 = 0.0;
  double travelTime = 0.0;
  double omegatau = 0.0;
  double rhoSq = 0.0;
  double *A_v = (double *)calloc(numSat * 4, sizeof(double));
  double det = 0.0;
  double detTol = 1e-12;
  double trop = 0.0;
  double *omc_v = (double *)calloc(numSat, sizeof(double));
  int ret = 0;
  double *az = (double *)calloc(numSat, sizeof(double));
  double *el = (double *)calloc(numSat, sizeof(double));
  double D = 0.0;
  double X[] = {0,0,0};
  double dx[] = {0,0,0};
  double normX = 100.0;
  int iter = 0;

  // Try initializing pos_v
  //pos_v[0] = (double)sdrini.xu0_v[0];
  //pos_v[1] = (double)sdrini.xu0_v[1];
  //pos_v[2] = (double)sdrini.xu0_v[2];
  //pos_v[0] = sdrstat.xyzdt[0];
  //pos_v[1] = sdrstat.xyzdt[1];
  //pos_v[2] = sdrstat.xyzdt[2];
  //pos_v[3] = sdrstat.xyzdt[3]; // leave out?

  // Initialization
  //int nmbOfIterations = 10;

  // Loop through BLS equation 10 times
  // Note: Should revise to exit on tolerance (norm of x small)
  for (int j=0; j<10; j++) {
  //while (normX > 1.0e-10) {

    // Break out of BLS loop if desired accuracy is achieved
    if (normX < 1.0e-10) {
      //printf("breaking from BLS loop early, iteration %d\n",j);
      break;
    }

    // Loop through obs and rotate xs, solve for el angle, and solve for tropo
    // correction.
    for (int i=0; i<numSat; i++) {
      // For first iteration, set Rot_X and trop
      if (iter==0) {
        for (int k=0; k<3; k++) {
          Rot_X_v[k] = X_v[i*3+k];
        }
        trop = 2.0;
      } else {

      // Solve for estimated range squared
      rho2 = ( ( (X_v[i*3+0] - pos_v[0]) * (X_v[i*3+0] - pos_v[0]) ) +
               ( (X_v[i*3+1] - pos_v[1]) * (X_v[i*3+1] - pos_v[1]) ) +
               ( (X_v[i*3+2] - pos_v[2]) * (X_v[i*3+2] - pos_v[2]) ) );

      // Solve for travelTime
      travelTime = sqrt(rho2) / CTIME;

      // Rotate SV position to account for Earth rotation during
      // travel time
      omegatau =  OMEGAEDOT * travelTime;
      Rot_X_v[0] =  (cos(omegatau) * X_v[i*3+0])  + (sin(omegatau) * X_v[i*3+1]);
      Rot_X_v[1] = (-sin(omegatau) * X_v[i*3+0]) + (cos(omegatau) * X_v[i*3+1]);
      Rot_X_v[2] = X_v[i*3+2];

      // Update SV el angles
      for (int i=0; i<3; i++) {
        X[i] = pos_v[i];
        dx[i] = Rot_X_v[i] - pos_v[i];
      }
      ret = topocent(X, dx, &az[i], &el[i], &D);
      if (ret!=0) {
        printf("topocent function: error\n");
      }

      // Load SV positions and angles for current SV into obs
      mlock(hobsvecmtx);
      sdrstat.obs_v[i*10+1] = Rot_X_v[0];
      sdrstat.obs_v[i*10+2] = Rot_X_v[1];
      sdrstat.obs_v[i*10+3] = Rot_X_v[2];
      sdrstat.obs_v[i*10+8] = az[i];
      sdrstat.obs_v[i*10+9] = el[i];
      unmlock(hobsvecmtx);

      //printf("iter: %d, az: %.1f, el: %.1f, D: %.1f\n",
      //  iter, az[i], el[i], D);

      // Calculate tropo correction (do later)
      ret = tropo(sin(el[i]*D2R), 0.0, 1013.0, 293.0, 50.0,
                                  0.0, 0.0, 0.0, &trop);
      if (ret!=0) {
        printf("tropo function: error\n");
      }

      } // end of if else

      // Apply the corrections to the observations to get estimated
      // range squared
      rhoSq = ( ( (Rot_X_v[0] - pos_v[0]) * (Rot_X_v[0] - pos_v[0]) ) +
                ( (Rot_X_v[1] - pos_v[1]) * (Rot_X_v[1] - pos_v[1]) ) +
                ( (Rot_X_v[2] - pos_v[2]) * (Rot_X_v[2] - pos_v[2]) ) );

      // Calculate the measurement residuals, PR observed minus PR calc
      omc_v[i] = pr_v[i] - sqrt(rhoSq) - pos_v[3] - trop;

      // Construct the A matrix
      A_v[i*4+0] = (-(Rot_X_v[0] - pos_v[0])) / sqrt(rhoSq);
      A_v[i*4+1] = (-(Rot_X_v[1] - pos_v[1])) / sqrt(rhoSq);
      A_v[i*4+2] = (-(Rot_X_v[2] - pos_v[2])) / sqrt(rhoSq);
      A_v[i*4+3] = 1.0;

    } // end of numSat loop

    // Exit if bad rank

    // Perform BLS to solve for error states
    //   dx = inv(A' * A) * A' * dz
    tempMat1 = nml_mat_from(numSat,4,numSat*4,A_v); // A: nx4
    tempMat2 = nml_mat_transp(tempMat1);            // A': 4xn
    tempMat3 = nml_mat_dot(tempMat2, tempMat1);     // A'A: 4x4
    tempMat4 = nml_mat_lup_solve(tempMat3);         // lup(A'A): 4x4

    // Check for viable inverse (if det equals zero, rank too low)
    det = nml_mat_det(tempMat4);
    if (fabs(det)<detTol) {
      printf("Exiting estRcvrPosn, determinant=%lf\n", det);
      goto errorDetected;
    }

    tempMat5 = nml_mat_inv(tempMat4);               // inv(A'A): 4x4
    tempMat6 = nml_mat_dot(tempMat5,tempMat2);      // inv(A'A)A': 4xn
    omc      = nml_mat_from(numSat,1,numSat,omc_v); // omc: nx1
    x        = nml_mat_dot(tempMat6,omc);           // inv(A'A)A'dz: 4x1

    // Apply error states to whole states
    tempPos = nml_mat_from(4,1,4,pos_v);             // tempPos: 4x1
    pos = nml_mat_add(tempPos,x);                   // pos: 4x1

    // Reset pos_v for next iteration
    pos_v[0] = pos->data[0][0];
    pos_v[1] = pos->data[1][0];
    pos_v[2] = pos->data[2][0];
    pos_v[3] = pos->data[3][0];

    // Calculate norm of X
    normX = sqrt( (x->data[0][0]*x->data[0][0]) +
                  (x->data[1][1]*x->data[1][1]) +
                  (x->data[2][2]*x->data[2][2]) +
                  (x->data[3][3]*x->data[3][3]) );

    // Increment iter counter
    iter = iter + 1;

  } // end the estimator loop

  // Set pvtflag
  sdrstat.pvtflag = 1;

  // Save pos values to xyzdt_v
  xyzdt_v[0] = pos->data[0][0];
  xyzdt_v[1] = pos->data[1][0];
  xyzdt_v[2] = pos->data[2][0];
  xyzdt_v[3] = pos->data[3][0];

  // Calculate DOP values (DOP is sqrt of the trace)
  *gdop =  sqrt( tempMat5->data[0][0] + tempMat5->data[1][1] +
                 tempMat5->data[2][2] + tempMat5->data[3][3] );

  // Write residuals to sdrstat so GUI can display
  for (int i=0; i<numSat; i++) {
    sdrstat.vk1_v[i] = omc_v[i];
  }

  // Free matrices
  nml_mat_free(x);
  nml_mat_free(A);
  nml_mat_free(pos);
  nml_mat_free(omc);
  nml_mat_free(tempMat1);
  nml_mat_free(tempMat2);
  nml_mat_free(tempMat3);
  nml_mat_free(tempMat5);
  nml_mat_free(tempMat6);
  nml_mat_free(tempPos);
  nml_mat_lup_free(tempMat4);

  // Free vectors (arrays)
  //free(X_v);
  free(A_v);
  free(omc_v);
  free(az);
  free(el);

  // Normal return with good solution
  return 0;

errorDetected:
  // Save output values as zero
  xyzdt_v[0] = 0.0;
  xyzdt_v[1] = 0.0;
  xyzdt_v[2] = 0.0;
  xyzdt_v[3] = 0.0;
  *gdop = 0.0;

  // Free matrices
  nml_mat_free(x);
  nml_mat_free(A);
  nml_mat_free(pos);
  nml_mat_free(omc);
  nml_mat_free(tempMat1);
  nml_mat_free(tempMat2);
  nml_mat_free(tempMat3);
  nml_mat_free(tempMat5);
  nml_mat_free(tempMat6);
  nml_mat_free(tempPos);
  nml_mat_lup_free(tempMat4);

  // Free vectors (arrays)
  //free(X_v);
  free(A_v);
  free(omc_v);
  free(az);
  free(el);

  // Error return without good solution
  sdrstat.pvtflag = 0;
  return -1;

} // end function

//-----------------------------------------------------------------------------
// Function to account for week crossover
//-----------------------------------------------------------------------------
extern void check_t(double time, double *corrTime) {

 // Initial guess for latitude
 double halfweek = 302400.0;

 // Initialize corrTime
 *corrTime = time;

 // Iterative process to find the geodetic latitude
 if (time > halfweek) {
  *corrTime = time - (2 * halfweek);
 }  else if (time < -halfweek) {
  *corrTime = time + (2 * halfweek);
 }
} // end function

//-----------------------------------------------------------------------------
// Function to convert ECEF to geodetic coordinates
//-----------------------------------------------------------------------------
extern void ecef2lla(double x, double y, double z,
                    double *lambda, double *phi, double *height)
{
 // WGS-84 ellipsoid parameters
 double A = 6378137.0;  // Semi-major axis (meters)
 double F = 1.0 / 298.257223563;  // Flattening
 double E = sqrt( (2*F) - (F*F) );

 // Initial guess for latitude
 *lambda = atan2(y,x);
 double p = sqrt(x * x + y * y);

 // Initial value of phi assuming h = 0
 *height = 0.0;
 *phi = atan2(z, p*(1.0 - E*E));
 double N = A / sqrt( (1.0 - sin(*phi)) * (1.0 - sin(*phi)) );
 double delta_h = 1000000;

 // Iterative process to find the geodetic latitude
 while (delta_h > 0.01) {
   double prev_h = *height;
   *phi = atan2(z, p*(1-(E*E) * (N / (N + *height) ) ) );
   N = A / sqrt( (1 - ( (E * sin(*phi)) * (E * sin(*phi)) ) ) );
   *height = (p / cos(*phi)) - N;
   delta_h = fabs(*height - prev_h);
 }
} // end function

//-----------------------------------------------------------------------------
// Function to calculate satellite position
//-----------------------------------------------------------------------------
extern int satPos(sdreph_t *sdreph, double transmitTime, double svPos[3],
            double *svClkCorr)
{
  // Sample ephemeris parameters (in order of 1 to 21). These values are
  // from the Thompson paper "Computing GPS Velocity and Acceleration from
  // the Broadcast Navigation Message."
  //int sat                         = sdreph->eph.sat; // not used
  double toe                      = sdreph->eph.toes;
  double toc						= toe;
  //double sqrta                    = sdreph->eph.A;
  double sqrta                    = sqrt(sdreph->eph.A);
  double e                        = sdreph->eph.e;
  double M0                       = sdreph->eph.M0;
  double omega                    = sdreph->eph.omg;
  double i0                       = sdreph->eph.i0;
  double Omega0                   = sdreph->eph.OMG0;
  double deltan                   = sdreph->eph.deln;
  double idot                     = sdreph->eph.idot;
  double Omegadot                 = sdreph->eph.OMGd;
  double cuc                      = sdreph->eph.cuc;
  double cus                      = sdreph->eph.cus;
  double crc                      = sdreph->eph.crc;
  double crs                      = sdreph->eph.crs;
  double cic                      = sdreph->eph.cic;
  double cis                      = sdreph->eph.cis;
  double af0                      = sdreph->eph.f0;
  double af1                      = sdreph->eph.f1;
  double af2                      = sdreph->eph.f2;
  double tgd                      = sdreph->eph.tgd[0];

  /*
  printf("sat: %d\n", sat);
  printf("transmitTime: %e, toe: %e, toc: %e\n", transmitTime, toe, toc);
  printf("sqrtA: %e, e: %e, m0: %e, omega: %e\n", sqrta, e, M0, omega);
  printf("i0: %e, Omega0: %e, deltan: %e, idot: %e\n", i0, Omega0, deltan, idot);
  printf("Omegadot: %e, cuc: %e, cus: %e, crc: %e\n", Omegadot, cuc, cus, crc);
  printf("crs: %e, cic: %e, cis: %e, af0: %e\n", crs, cic, cis, af0);
  printf("af1: %e, af2: %e, tgd: %e\n", af1,af2,tgd);
  //*/

  // Local parameters
  double t = transmitTime; // rcvr_tow
  double A;
  double n0, n;
  double tk, tc;
  double Mk, vk, ik, phik, u_k, rk, Omega_k;
  double delta_uk, delta_rk, delta_ik;
  double x_kp, y_kp;
  double xk, yk, zk;
  double E0, Ek, dtr, dt;
  int ii;

  // Initialize svPos to zero
  svPos[0] = 0.0; svPos[1] = 0.0; svPos[2] = 0.0;

  // Velocity terms (if velocity calculated)
  //long double ekdot, vkdot, ukdot, ikdot, rkdot, omegakdot;
  //long double xpkdot, ypkdot;
  //long double xkdot, ykdot, zkdot;

  // Semi-major axis
  A = sqrta * sqrta;           //sqrta is the square root of A

  // Computed mean motion
  n0 = sqrt( MU / (A*A*A) );

  // Time from ephemeris, reference clock
  tk = t - toe;              // t is the time of the pos. & vel. request.
  if (tk > 302400) {
    tk = tk - (2 * 302400);
  }
  else if (tk < -302400) {
  	tk = tk + (604800);
  }

  // Corrected mean motion
  n = n0 + deltan;

  // Mean anomaly at t
  Mk = M0 + n*tk;

  // Kepler equation for eccentric anomaly
  E0 = Mk;

  for (ii=0; ii<3; ii++){
    Ek = E0 + (Mk - E0 + e * sin(E0)) / (1 - e * cos(E0));
    E0 = Ek;
  }

  //In the line, below, tak is the true anomaly (which is nu in the ICD-200).
  vk = 2 * atan( sqrt( (1.0 + e) / (1.0 - e) ) * tan(Ek / 2) );

  // Argument of latitude
  phik = omega + vk;

  // Second harmonic perturbations
  delta_uk = cus * sin(2.0 * phik) + cuc * cos(2.0 * phik);
  delta_rk = crs * sin(2.0 * phik) + crc * cos(2.0 * phik);
  delta_ik = cis * sin(2.0 * phik) + cic * cos(2.0 * phik);

  // Corrected argument of latitude
  u_k = phik + delta_uk;

  // Corrected radius
  rk = A * (1.0 - (e * cos(Ek))) + delta_rk;

  // Corrected inclination
  ik = i0 + (idot * tk) + delta_ik;

  // Cos and sin values for rk
  x_kp = rk * cos(u_k);
  y_kp = rk * sin(u_k);

  // Longitude of ascending node
  Omega_k = Omega0 + (Omegadot - OMEGAEDOT) * tk - (OMEGAEDOT * toe);

  // SV position in ECEF
  xk = (x_kp * cos(Omega_k)) - (y_kp * sin(Omega_k) * cos(ik));
  yk = (x_kp * sin(Omega_k)) + (y_kp * cos(Omega_k) * cos(ik));
  zk =                          y_kp * sin(ik);

  if ( isnan(xk) || isnan(yk) || isnan(zk) ) {
    goto errorDetected;
  }

  // Load svPos vector
  svPos[0] = xk; svPos[1] = yk; svPos[2] = zk;

  //-----------------------------------------------------------------------------
  // SV clock correction calculation
  //-----------------------------------------------------------------------------

  // Time from ephemeris, SV clock
  tc = t - toc;              // t is the time of the pos. & vel. request.
  if (tc > 302400) {
    tc = tc - (2 * 302400);
  }
  else if (tc < -302400) {
    tc = tc + (604800);
  }

  // Relativistic SV clock correction correction
  dtr = -2 * sqrt(MU) / (CTIME * CTIME) * e * sqrta * sin(Ek);

  // SV clock bias with group delay (tgd) removed and relativistic correction
  // (dtr) added. the tgd term is included for L1-only signals, so it would
  // be removed if dual L1/L2 signals are used.
  dt = af0 + (af1 * tc) + (af2 * (tc * tc)) - tgd + dtr;
  *svClkCorr = dt;

  // Normal return
  return 0;

// Exit with error detected
errorDetected:
  return -1;

} // end of function

//-----------------------------------------------------------------------------
// Rotation matrix function
//-----------------------------------------------------------------------------
extern void rot(double R[9], double angle, int axis) {
 R[0] = 1.0; R[4] = 1.0; R[8] = 1.0;
 R[1] = 0.0; R[3] = 0.0; R[6] = 0.0;
 R[2] = 0.0; R[5] = 0.0; R[7] = 0.0;

 double cang, sang;
 cang = cos(angle * M_PI/180.0);
 sang = sin(angle * M_PI/180.0);

 if (axis == 1) {
  R[4] = cang;
  R[8] = cang;
  R[5] = sang;
  R[7] = -sang;
 }
 if (axis == 2) {
  R[0] = cang;
  R[8] = cang;
  R[2] = -sang;
  R[6] = sang;
 }
 if (axis == 3) {
  R[0] = cang;
  R[4] = cang;
  R[3] = -sang;
  R[1] = sang;
 }
} // end function

//-----------------------------------------------------------------------------
// precheckObs
//
// Input: index (the index of row we wish to remove, starting at 0
//-----------------------------------------------------------------------------
extern void precheckObs()
{
  // Check pseudoranges to ensure within reasonable range. Acceptable ranges
  // typically 60 to 92ms, or 64 to 85ms if tighter limits desired. Also
  // check for low SNR

  int index = 0;
  int nsat_temp = 0;
  double obs_v_temp[10*MAXSAT] = {0.0};

  mlock(hobsvecmtx);
  int nsat = sdrstat.nsat;
  nsat_temp = nsat;
  for (int i=0;i<nsat; i++){
    if ( ((int)sdrstat.obs_v[i*10+0]==0) ||
         (sdrstat.obs_v[i*10+4]<LOW_PR) ||
         (sdrstat.obs_v[i*10+4]>HIGH_PR) ||
         (sdrstat.obs_v[i*10+6]<GPS_WEEK) ||
         (sdrstat.obs_v[i*10+7]<SNR_PVT_THRES) ) {
         //(sdrstat.elapsedTime>ET_TIMER && sdrstat.obs_v[i*10+9]<SV_EL_MASK) ){
      nsat_temp = nsat_temp - 1;
      //printf("Row of obs_v removed due to invalid obs\n");
    } else {
      // If obs is valid, load it into the temp obs array
      for (int j=0; j<10; j++) {
        obs_v_temp[index+j] = sdrstat.obs_v[i*10+j];
      }
      index = index + 10;
    }
  } // end for

  // Update nsat, and load temp obs array back into obs_v
  sdrstat.nsat = nsat_temp;
  memcpy(sdrstat.obs_v, obs_v_temp, sizeof(obs_v_temp));
  unmlock(hobsvecmtx);

} // end function

//-----------------------------------------------------------------------------
// checkSvEl
//
// Input: index (the index of row we wish to remove, starting at 0
//-----------------------------------------------------------------------------
extern int checkSvEl()
{
  // Check SV elevations to make sure they are above the elevation mask

  int nsat_temp = 0;
  double obs_v_temp[10*MAXSAT] = {0.0};
  int lowElFlag = 0;

  mlock(hobsvecmtx);
  int nsat = sdrstat.nsat;
  nsat_temp = 0;

  // First check to see if any obs have low elevations
  for (int i=0;i<nsat; i++){
    if ((sdrstat.obs_v[i*10+9]>=SV_EL_MASK) ||
        (fabs(sdrstat.obs_v[i*10+9])<1e-10)) {
      // elevation is good, or is 0.0 (not solved for yet), so do nothing
    } else {
      lowElFlag = 1;
      //printf("Low elevation detected for G%02d El:%.1f\n",
      //        (int)(sdrstat.obs_v[i*10+0]),sdrstat.obs_v[i*10+9]);
    }
  }
  unmlock(hobsvecmtx);

  if (lowElFlag==0) {
    goto endFunction; // Elevations are all good, so return from function
  }

// If low el obs have been detected, remove them
if (lowElFlag==1) {
  mlock(hobsvecmtx);
  // If we get here, at least one SV is below the elevation mask, Check all.
  // If the el angle appears to be zero, then leave the obs in the set to
  // give the estimator time to get non-zero value (note: this last step may
  // not be needed).
  //printf("\n\n");
  for (int i=0; i<nsat; i++) {
    if ((sdrstat.obs_v[i*10+9]>=SV_EL_MASK) || (fabs(sdrstat.obs_v[i*10+9])<1e-10)) {
      // Obs is valid, so copy it to temp variable
      for (int j=0; j<10; j++) {
        obs_v_temp[(nsat_temp*10)+j] = sdrstat.obs_v[i*10+j];
        //printf("G%02d  El:%.1f  SV_EL_MASK:%.1f\n",
        //   (int)sdrstat.obs_v[i*10+0], sdrstat.obs_v[i*10+9], SV_EL_MASK);
      }
      nsat_temp = nsat_temp + 1;
    } // end if
  } // end for

  // Update nsat, and load temp obs array back into obs_v
  sdrstat.nsat = nsat_temp;
  // Zero out obs_v
  memset(sdrstat.obs_v, 0, sizeof(sdrstat.obs_v));
  //for (int n=0;n<32*10;n++) {
  //  sdrstat.obs_v[n] = 0.0;
  //}
  // Load updated obs into obs_v
  for (int j=0; j<nsat; j++) {
    for (int k=0; k<10; k++) {
      sdrstat.obs_v[j*10+k] = obs_v_temp[j*10+k];
    }
  }
  //memcpy(sdrstat.obs_v, obs_v_temp, sizeof(obs_v_temp));
  unmlock(hobsvecmtx);
}

  return 0;

endFunction:
  return 0;

} // end function

/*
//-----------------------------------------------------------------------------
// checkSvEl
//
// Input: index (the index of row we wish to remove, starting at 0
//-----------------------------------------------------------------------------
extern void checkSvEl()
{
  // Check SV elevations to make sure they are above the elevation mask

  int index = 0;
  int nsat_temp = 0;
  double obs_v_temp[10*MAXSAT] = {0.0};

  mlock(hobsvecmtx);
  int nsat = sdrstat.nsat;
  nsat_temp = nsat;

  for (int i=0;i<nsat; i++){
    // If the SV is below the elevation mask, remove it from the obs set. Also,
    // if the el angle appears to be zero, then leave the obs in the set to
    // give the estimator time to get non-zero value (note: this las step may
    // not be needed).
    if (sdrstat.obs_v[i*10+9]<SV_EL_MASK && fabs(sdrstat.obs_v[i*10+9])>1e-12) {
      nsat_temp = nsat_temp - 1;
      //printf("Row of obs_v removed due to invalid obs\n");
      printf("checkSvEl: PRN %d obs removed with SvEl %.1f degrees\n",
           (int)sdrstat.obs_v[i*10+0], sdrstat.obs_v[i*10+9]);
    } else {
      // If obs is valid, load it into the temp obs array
      for (int j=0; j<10; j++) {
        obs_v_temp[index+j] = sdrstat.obs_v[i*10+j];
      }
      index = index + 10;
    } // end if else
  } // end for

  // Update nsat, and load temp obs array back into obs_v
  sdrstat.nsat = nsat_temp;
  memcpy(sdrstat.obs_v, obs_v_temp, sizeof(obs_v_temp));
  unmlock(hobsvecmtx);

} // end function
*/

//-----------------------------------------------------------------------------
// precheckEph
//
// Input: index (the index of row we wish to remove, starting at 0
//-----------------------------------------------------------------------------
extern void precheckEph()
{
  // Check GPS week and eph values to make sure they seem somewhat reasonable,
  // i.e. non-zero. Need to add smarts to the checks below, as in what to do
  // if new eph comes along and it is bad.
  double tol = 1e-15; // tolerance for checking if non-zero

  int index = 0;
  int nsat_temp = 0;
  double obs_v_temp[10*MAXSAT] = {0.0};

  mlock(hobsvecmtx);
  int nsat = sdrstat.nsat;
  nsat_temp = nsat;
  for (int n=0; n<nsat; n++){
    if ( ((int)sdrstat.obs_v[n*10+6]==0) ||
         (sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.toes<1.0) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.A)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.e)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.M0)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.omg)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.i0)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.OMG0)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.deln)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.idot)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.OMGd)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.cuc)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.cus)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.crc)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.crs)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.cic)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.cis)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.f0)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.f1)<tol) ||
         (fabs(sdrch[(int)sdrstat.obs_v[n*10]-1].nav.sdreph.eph.tgd[0])<tol) ){
      nsat_temp = nsat_temp - 1;
      //printf("Row of obs_v removed due to invalid eph\n");
    } else {
      for (int j=0; j<10; j++) {
        obs_v_temp[index+j] = sdrstat.obs_v[n*10+j];
      }
      index = index + 10;
    } // end else
  } // end for

  // Update nsat and obs_v
  sdrstat.nsat = nsat_temp;
  memcpy(sdrstat.obs_v, obs_v_temp, sizeof(obs_v_temp));
  unmlock(hobsvecmtx);

} // end function

//-----------------------------------------------------------------------------
// Estimate tropo correction
//-----------------------------------------------------------------------------
extern int tropo(double sinel, double hsta, double p, double tkel,
                 double hum, double hp, double htkel, double hhum,
                 double *ddr) {

//TROPO  Calculation of tropospheric correction.
//       The range correction ddr in m is to be subtracted from
//       pseudo-ranges and carrier phases
//
// ddr = tropo(sinel, hsta, p, tkel, hum, hp, htkel, hhum);
//
//   Inputs:
//       sinel   - sin of elevation angle of satellite
//       hsta    - height of station in km
//       p       - atmospheric pressure in mb at height hp
//       tkel    - surface temperature in degrees Kelvin at height htkel
//       hum     - humidity in % at height hhum
//       hp      - height of pressure measurement in km
//       htkel   - height of temperature measurement in km
//       hhum    - height of humidity measurement in km
//
//   Outputs:
//       ddr     - range correction (meters)
//
// Reference
// Goad, C.C. & Goodman, L. (1974) A Modified Tropospheric
// Refraction Correction Model. Paper presented at the
// American Geophysical Union Annual Fall Meeting, San
// Francisco, December 12-17

// A Matlab reimplementation of a C code from driver.
// Kai Borre 06-28-95
//
// CVS record:
// $Id: tropo.m,v 1.1.1.1.2.4 2006/08/22 13:46:00 dpl Exp $
//==========================================================================

double a_e    = 6378.137;     // semi-major axis of earth ellipsoid
double b0     = 7.839257e-5;
double tlapse = -6.5;
double tkhum  = tkel + tlapse * (hhum - htkel);
double atkel  = 7.5*(tkhum-273.15) / (237.3+tkhum-273.15);
double e0     = 0.0611 * hum * pow(10,atkel);
double tksea  = tkel - tlapse * htkel;
double em     = -978.77 / (2.8704e6*tlapse*1.0e-5);
double tkelh  = tksea + tlapse*hhum;
double e0sea  = e0 * pow((tksea/tkelh),(4*em));
double tkelp  = tksea + tlapse*hp;
double psea   = p * pow((tksea/tkelp),em);
double alpha[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double rn[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

if (sinel < 0) {
    sinel = 0;
}

double tropo   = 0;
int done    = 0;
double refsea  = 77.624e-6 / tksea;
double htop    = 1.1385e-5 / refsea;
refsea  = refsea * psea;
double ref     = refsea * pow(((htop-hsta)/htop),4);

while (1) {
    double rtop = pow((a_e+htop),2) - pow((a_e+hsta),2)*(1-pow(sinel,2));

    // check to see if geometry is crazy
    if (rtop < 0) {
        rtop = 0;
    }

    rtop = sqrt(rtop) - (a_e+hsta)*sinel;
    double a    = -sinel/(htop-hsta);
    double b    = -b0*(1-pow(sinel,2)) / (htop-hsta);
    //double rn   = zeros(8,1);

    for (int i=0; i<7; i++) {
        rn[i] = pow(rtop,(i+2));
    }

    alpha[0] = 2 * a;
    alpha[1] = 2 * pow(a,2) + 4 * b/3;
    alpha[2] = a * (pow(a,2) + 3 * b);
    alpha[3] = pow(a,4)/5 + 2.4 * pow(a,2) * b + 1.2 * pow(b,2);
    alpha[4] = 2 * a * b * (pow(a,2) + 3 * b) / 3;
    alpha[5] = pow(b,2) * (6 * pow(a,2) + 4 * b) * 1.428571e-1;
    alpha[6] =     0.0;
    alpha[7] =     0.0;

    if (pow(b,2) > 1.0e-35) {
        alpha[6] = a*pow(b,3)/2;
        alpha[7] = pow(b,4)/9;
    }

    double dr = rtop;
    // Calculate dr = dr + alpha * rn;
    for (int i=0; i<7; i++) {
      dr = dr + alpha[i] * rn[i];
    }

    tropo = tropo + dr * ref * 1000;

    if (done == 1) {
        *ddr = tropo;
        break;
    }

    done    = 1;
    refsea  = (371900.0e-6 / tksea - 12.92e-6) / tksea;
    htop    = 1.1385e-5 * (1255 / tksea + 0.05) / refsea;
    ref     = refsea * e0sea * pow(((htop-hsta)/htop),4);

} // end while loop

// end of function
return 0;

} // end function

//-----------------------------------------------------------------------------
// Estimate LLA given ECEF
//-----------------------------------------------------------------------------
extern int togeod(double a, double finv, double X, double Y, double Z,
                  double *dphi, double *dlambda, double *h) {
//TOGEOD   Subroutine to calculate geodetic coordinates latitude, longitude,
//         height given Cartesian coordinates X,Y,Z, and reference ellipsoid
//         values semi-major axis (a) and the inverse of flattening (finv).
//
//[dphi, dlambda, h] = togeod(a, finv, X, Y, Z);
//
//  The units of linear parameters X,Y,Z,a must all agree (m,km,mi,ft,..etc)
//  The output units of angular quantities will be in decimal degrees
//  (15.5 degrees not 15 deg 30 min). The output units of h will be the
//  same as the units of X,Y,Z,a.
//
//   Inputs:
//       a           - semi-major axis of the reference ellipsoid
//       finv        - inverse of flattening of the reference ellipsoid
//       X,Y,Z       - Cartesian coordinates
//
//   Outputs:
//       dphi        - latitude
//       dlambda     - longitude
//       h           - height above reference ellipsoid

//  Copyright (C) 1987 C. Goad, Columbus, Ohio
//  Reprinted with permission of author, 1996
//  Fortran code translated into MATLAB
//  Kai Borre 03-30-96
//
// CVS record:
// $Id: togeod.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
//==========================================================================

  *h       = 0;
  double tolsq   = 1.e-10;
  double maxit   = 50;
  double esq = 0.0;
  double sinphi = 0.0;

  // compute radians-to-degree factor
  //double rtd     = 180/PI;

  // compute square of eccentricity
  if (finv < 1.e-20) {
    esq = 0;
  } else {
    esq = (2 - 1/finv) / finv;
  }

  double oneesq  = 1 - esq;

  // first guess
  // P is distance from spin axis
  double P = sqrt(X*X+Y*Y);
  // direct calculation of longitude

  if (P > 1.e-20) {
    *dlambda = atan2(Y,X) * R2D;
  } else {
    *dlambda = 0;
  }

  if (*dlambda < 0) {
    *dlambda = *dlambda + 360;
  }

  // r is distance from origin (0,0,0)
  double r = sqrt(P*P + Z*Z);

  if (r > 1.e-20) {
    sinphi = Z/r;
  } else {
    sinphi = 0;
  }

  *dphi = asin(sinphi);

  // initial value of height  =  distance from origin minus
  // approximate distance from origin to surface of ellipsoid
  if (r < 1.e-20) {
    *h = 0;
    goto errorDetected;
  }

  *h = r - a * (1 - sinphi * sinphi / finv);

  // iterate
  for (int i=0; i<maxit; i++) {
    double sinphi  = sin(*dphi);
    double cosphi  = cos(*dphi);

    // compute radius of curvature in prime vertical direction
    double N_phi   = a / sqrt(1 - esq * sinphi * sinphi);

    // compute residuals in P and Z
    double dP      = P - (N_phi + *h) * cosphi;
    double dZ      = Z - (N_phi * oneesq + *h) * sinphi;

    // update height and latitude
    *h       = *h + (sinphi*dZ + cosphi*dP);
    *dphi    = *dphi + (cosphi*dZ - sinphi*dP)/(N_phi + *h);

    // test for convergence
    if (dP*dP + dZ*dZ < tolsq) {
        break;
    }

    // Not Converged--Warn user
    if (i == maxit) {
        printf("Problem in TOGEOD, did not converge in %2d iterations\n", i);
    }
  } // for i = 1:maxit

  *dphi = *dphi * R2D;

  // end togeod.m
  return 0;

errorDetected:
    printf("togeod--- errorDetected\n");
    return -1;

} // end function

//-----------------------------------------------------------------------------
// Estimate Az and El given ECEF vector
//-----------------------------------------------------------------------------
extern int topocent(double X[], double dx[], double *Az, double *El, double *D) {
//TOPOCENT  Transformation of vector dx into topocentric coordinate
//          system with origin at X.
//          Both parameters are 3 by 1 vectors.
//
//[Az, El, D] = topocent(X, dx);
//
//   Inputs:
//       X           - vector origin corrdinates (in ECEF system [X; Y; Z;])
//       dx          - vector ([dX; dY; dZ;]).
//
//   Outputs:
//       D           - vector length. Units like units of the input
//       Az          - azimuth from north positive clockwise, degrees
//       El          - elevation angle, degrees

//Kai Borre 11-24-96
//Copyright (c) by Kai Borre
//
// CVS record:
// $Id: topocent.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
//==========================================================================

int ret;
//double F_v[9] = {0};
double local_v[3] = {0};
double phi;
double lambda;
double h;

ret = togeod(6378137, 298.257223563, X[0], X[1], X[2], &phi, &lambda, &h);
if (ret!=0) {
 printf("togeod function: error\n");
}

double cl  = cos(lambda * D2R);
double sl  = sin(lambda * D2R);
double cb  = cos(phi * D2R);
double sb  = sin(phi * D2R);

// Definition of F
//double F   = [-sl -sb*cl cb*cl;
//        cl -sb*sl cb*sl;
//        0    cb   sb];

// Calculate local_v = F' * dx;
local_v[0] = -sl * dx[0] + cl * dx[1] + 0.0 * dx[2];
local_v[1] = -sb*cl * dx[0] - sb*sl * dx[1] + cb * dx[2];
local_v[2] = cb*cl * dx[0] + cb*sl * dx[1] + sb * dx[2];

double E   = local_v[0];
double N   = local_v[1];
double U   = local_v[2];

double hor_dis = sqrt(E*E + N*N);

if (hor_dis < 1.e-20) {
  *Az = 0.0;
  *El = 90;
} else {
  *Az = atan2(E, N)/D2R;
  *El = atan2(U, hor_dis)/D2R;
}

if (*Az < 0) {
  *Az = *Az + 360;
}

*D   = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

return 0;
}  // end of function
