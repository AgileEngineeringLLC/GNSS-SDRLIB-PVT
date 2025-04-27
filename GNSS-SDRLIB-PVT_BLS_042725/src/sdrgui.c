//----------------------------------------------------------------------------
//
//  sdrgui.c
//
//  Copyright 2025, Don Kelly, don.kelly@mac.com
//----------------------------------------------------------------------------

#include "sdr.h"

//-----------------------------------------------------------------------------
// GUI Thread
//-----------------------------------------------------------------------------
extern void *guithread(void *arg) {

  double elapsedTime = 0.0;

  // Buffers for PRN data and LLA
  char buffer1[80] = " ";
  char buffer2[80] = " ";
  char buffer3[80] = " ";
  char buffer4[80] = " ";
  char buffer5[100] = " ";

  // Buffers for obs info
  char buffer6[110];
  char buffer7[110];
  char buffer8[110];
  char buffer9[110];
  char buffer10[110];
  char buffer11[110];
  char buffer12[110];
  char buffer13[110];
  char buffer14[110];
  char buffer15[110];

  char buffer16[100]; // UTC display

  // Temp string
  char str1[10];

  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    fprintf(stderr, "SDL_Init Error: %s\n", SDL_GetError());
  }

  SDL_Window* window = SDL_CreateWindow("GPS Navigation Status",
         SCREEN_LOC_X, SCREEN_LOC_Y, SCREEN_WIDTH, SCREEN_HEIGHT, 0);
                         //SDL_WINDOWPOS_CENTERED,
                         //SDL_WINDOWPOS_CENTERED,
  if (window == NULL) {
    fprintf(stderr, "SDL_CreateWindow Error: %s\n", SDL_GetError());
    SDL_Quit();
  }

  SDL_Renderer* renderer = SDL_CreateRenderer(window, -1,
                              SDL_RENDERER_ACCELERATED);
  if (renderer == NULL) {
    fprintf(stderr, "SDL_CreateRenderer Error: %s\n", SDL_GetError());
    SDL_DestroyWindow(window);
    SDL_Quit();
  }

  if (TTF_Init() != 0) {
    fprintf(stderr, "TTF_Init Error: %s\n", TTF_GetError());
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
  }

  TTF_Font* font = TTF_OpenFont(sdrini.fontfile, FONT_SIZE);
  if (!font) {
    fprintf(stderr, "TTF_OpenFont Error: %s\n", TTF_GetError());
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
  }

  // Color options for text
  SDL_Color textColor =  {0,0,0,255}; // Black
  SDL_Color textColor1 = {0,192,0,255}; // Green
  //SDL_Color textColor2 = {255,0,0,255}; // Red
  SDL_Color textColor3 = {160,32,255,255}; // Purple
  SDL_Color textColor4 = {0,32,255,255}; // Blue
  //SDL_Color textColor5 = {255,224,32,255}; // Yellow
  //SDL_Color textColor6 = {255,255,255,255}; // White
  //SDL_Color textColor7 = {128,128,128,255}; // Gray

  // Declare textSurfaces
  SDL_Surface* textSurface1 =
      TTF_RenderText_Solid(font, buffer1, textColor);
  SDL_Surface* textSurface2 =
      TTF_RenderText_Solid(font, buffer2, textColor);
  SDL_Surface* textSurface3 =
      TTF_RenderText_Solid(font, buffer3, textColor);
  SDL_Surface* textSurface4 =
      TTF_RenderText_Solid(font, buffer4, textColor);
  SDL_Surface* textSurface5 =
      TTF_RenderText_Solid(font, buffer5, textColor);
  SDL_Surface* textSurface6 =
      TTF_RenderText_Solid(font, buffer6, textColor);
  SDL_Surface* textSurface7 =
      TTF_RenderText_Solid(font, buffer7, textColor);
  SDL_Surface* textSurface8 =
      TTF_RenderText_Solid(font, buffer8, textColor);
  SDL_Surface* textSurface9 =
      TTF_RenderText_Solid(font, buffer9, textColor);
  SDL_Surface* textSurface10 =
      TTF_RenderText_Solid(font, buffer10, textColor);
  SDL_Surface* textSurface11 =
      TTF_RenderText_Solid(font, buffer11, textColor);
  SDL_Surface* textSurface12 =
      TTF_RenderText_Solid(font, buffer12, textColor);
  SDL_Surface* textSurface13 =
      TTF_RenderText_Solid(font, buffer13, textColor);
  SDL_Surface* textSurface14 =
      TTF_RenderText_Solid(font, buffer14, textColor);
  SDL_Surface* textSurface15 =
      TTF_RenderText_Solid(font, buffer15, textColor);
  SDL_Surface* textSurface16 =
      TTF_RenderText_Solid(font, buffer16, textColor);

  // Verify textSurfaces created
  if (!textSurface1) {
    fprintf(stderr, "TTF_RenderText_Solid Error: %s\n", TTF_GetError());
    TTF_CloseFont(font);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
  }

  // Create textTextures from textSurfaces
  SDL_Texture* textTexture1 =
        SDL_CreateTextureFromSurface(renderer, textSurface1);
  SDL_Texture* textTexture2 =
        SDL_CreateTextureFromSurface(renderer, textSurface2);
  SDL_Texture* textTexture3 =
        SDL_CreateTextureFromSurface(renderer, textSurface3);
  SDL_Texture* textTexture4 =
        SDL_CreateTextureFromSurface(renderer, textSurface4);
  SDL_Texture* textTexture5 =
        SDL_CreateTextureFromSurface(renderer, textSurface5);
  SDL_Texture* textTexture6 =
        SDL_CreateTextureFromSurface(renderer, textSurface6);
  SDL_Texture* textTexture7 =
        SDL_CreateTextureFromSurface(renderer, textSurface7);
  SDL_Texture* textTexture8 =
        SDL_CreateTextureFromSurface(renderer, textSurface8);
  SDL_Texture* textTexture9 =
        SDL_CreateTextureFromSurface(renderer, textSurface9);
  SDL_Texture* textTexture10 =
        SDL_CreateTextureFromSurface(renderer, textSurface10);
  SDL_Texture* textTexture11 =
        SDL_CreateTextureFromSurface(renderer, textSurface11);
  SDL_Texture* textTexture12 =
        SDL_CreateTextureFromSurface(renderer, textSurface12);
  SDL_Texture* textTexture13 =
        SDL_CreateTextureFromSurface(renderer, textSurface13);
  SDL_Texture* textTexture14 =
        SDL_CreateTextureFromSurface(renderer, textSurface14);
  SDL_Texture* textTexture15 =
        SDL_CreateTextureFromSurface(renderer, textSurface15);
  SDL_Texture* textTexture16 =
        SDL_CreateTextureFromSurface(renderer, textSurface16);

  // Verify textTextures
  if (!textTexture1) {
    fprintf(stderr, "SDL_CreateTextureFromSurface error: %s\n",
              SDL_GetError());
    SDL_FreeSurface(textSurface1);
    TTF_CloseFont(font);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
  }

  // Define text boxes
  SDL_Rect textRect1;
  textRect1.x = 10;
  textRect1.y = 10;
  SDL_QueryTexture(textTexture1, NULL, NULL, &textRect1.w, &textRect1.h);

  SDL_Rect textRect2;
  textRect2.x = 10;
  textRect2.y = 100;
  SDL_QueryTexture(textTexture2, NULL, NULL, &textRect2.w, &textRect2.h);

  SDL_Rect textRect3;
  textRect3.x = 10;
  textRect3.y = 135;
  SDL_QueryTexture(textTexture3, NULL, NULL, &textRect3.w, &textRect3.h);

  SDL_Rect textRect4;
  textRect4.x = 10;
  textRect4.y = 170;
  SDL_QueryTexture(textTexture4, NULL, NULL, &textRect4.w, &textRect4.h);

  SDL_Rect textRect5;
  textRect5.x = 10;
  textRect5.y = 230;
  SDL_QueryTexture(textTexture5, NULL, NULL, &textRect5.w, &textRect5.h);

  SDL_Rect textRect6;
  textRect6.x = 10;
  textRect6.y = 290;
  SDL_Rect textRect7;
  textRect7.x = 10;
  textRect7.y = 330;
  SDL_Rect textRect8;
  textRect8.x = 10;
  textRect8.y = 370;
  SDL_Rect textRect9;
  textRect9.x = 10;
  textRect9.y = 410;
  SDL_Rect textRect10;
  textRect10.x = 10;
  textRect10.y = 450;
  SDL_Rect textRect11;
  textRect11.x = 10;
  textRect11.y = 490;
  SDL_Rect textRect12;
  textRect12.x = 10;
  textRect12.y = 530;
  SDL_Rect textRect13;
  textRect13.x = 10;
  textRect13.y = 570;
  SDL_Rect textRect14;
  textRect14.x = 10;
  textRect14.y = 610;
  SDL_Rect textRect15;
  textRect15.x = 10;
  textRect15.y = 650;

  SDL_Rect textRect16;
  textRect16.x = 10;
  textRect16.y = 50;

  SDL_Event event;
  //int quit = 0;

  // While loop for keeping window active
  while (!sdrstat.stopflag) {
//  while (!quit) {
    // Inner while loop for monitoring to quit
    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT) {
        sdrstat.stopflag=1;
        //quit = 1;
      }
    } // end of while (quit flag)

    // Try this to quit SDL2 window if keyboard 'q' received
    if (sdrstat.stopflag) {
      SDL_Quit();
    }

    // Render the window
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);

    // Update Runtime
    sprintf(buffer1, "Runtime:           %.1f", elapsedTime);
    SDL_Surface* textSurface1 =
        TTF_RenderText_Solid(font, buffer1, textColor);
    SDL_Texture* textTexture1 =
          SDL_CreateTextureFromSurface(renderer, textSurface1);
    SDL_QueryTexture(textTexture1, NULL, NULL, &textRect1.w, &textRect1.h);
    SDL_RenderCopy(renderer, textTexture1, NULL, &textRect1);
    SDL_DestroyTexture(textTexture1);
    SDL_FreeSurface(textSurface1);

    // Pull data from sdrstat and sdrch
    int prn[32] = {0};
    int flagacq[32] = {0};
    int flagsync[32] = {0};
    int flagdec[32] = {0};
    int nsat = 0;
    double lat = 0.0;
    double lon = 0.0;
    double hgt = 0.0;
    double gdop = 0.0;
    double clkBias = 0.0;
    double obs_v[MAXSAT*10] = {0.0};
    double vk1_v[MAXSAT] = {0.0};

    mlock(hobsvecmtx);
    for (int i=0; i<32; i++) {
      prn[i] = sdrch[i].prn;
      flagacq[i] = sdrch[i].flagacq;
      flagsync[i] = sdrch[i].nav.flagsync;
      flagdec[i] = sdrch[i].nav.flagdec;
    }
    nsat = sdrstat.nsat;
    lat = sdrstat.lat;
    lon = sdrstat.lon;
    hgt = sdrstat.hgt;
    gdop = sdrstat.gdop;
    clkBias = sdrstat.xyzdt[3];
    for (int m=0; m<nsat; m++) {
      for (int p=0; p<10; p++) {
        obs_v[m*10+p] = sdrstat.obs_v[m*10+p];
      }
      vk1_v[m] = sdrstat.vk1_v[m];
    }
    unmlock(hobsvecmtx);

    // UTC display
    int gps_week;
    double gps_tow;
    gps_week = (int)obs_v[6];

    // Correct rcvr TOW with rcvr clock bias for precise UTC
    mlock(hobsvecmtx);
    //gps_tow = sdrstat.obs_v[5] + sdrstat.clkBias/CTIME;
    gps_tow = obs_v[5] + clkBias/CTIME;
    unmlock(hobsvecmtx);
    time_t utc_time_seconds = gps_to_utc(gps_week, gps_tow);
    struct tm utc_tm;
    gmtime_r(&utc_time_seconds, &utc_tm);
    //strftime(buffer16, 100, "UTC Time:         %Y-%m-%d %H:%M:%S UTC", &utc_tm);
    sprintf(buffer16,"UTC Time:         %04d-%02d-%02d %02d:%02d:%02d.%03d",
       utc_tm.tm_year + 1900, utc_tm.tm_mon + 1, utc_tm.tm_mday,
       utc_tm.tm_hour, utc_tm.tm_min, utc_tm.tm_sec, (int)(gps_tow * 1000) % 1000);

    SDL_Surface* textSurface16 =
        TTF_RenderText_Solid(font, buffer16, textColor);
    SDL_Texture* textTexture16 =
          SDL_CreateTextureFromSurface(renderer, textSurface16);
    SDL_QueryTexture(textTexture16, NULL, NULL, &textRect16.w, &textRect16.h);
    SDL_RenderCopy(renderer, textTexture16, NULL, &textRect16);
    SDL_DestroyTexture(textTexture16);
    SDL_FreeSurface(textSurface16);

    // Update acquired PRNs
    sprintf(buffer2, "Acquired:         ");
    for (int i=0; i<32; i++) {
       if (flagacq[i]==1) {
         sprintf(str1, "%2d  ", prn[i]);
         strcat(buffer2, str1);
       }
    }
    SDL_Surface* textSurface2 =
        TTF_RenderText_Solid(font, buffer2, textColor4);
    SDL_Texture* textTexture2 =
          SDL_CreateTextureFromSurface(renderer, textSurface2);
    SDL_QueryTexture(textTexture2, NULL, NULL, &textRect2.w, &textRect2.h);
    SDL_RenderCopy(renderer, textTexture2, NULL, &textRect2);
    SDL_DestroyTexture(textTexture2);
    SDL_FreeSurface(textSurface2);

    // Update tracked PRNs
    sprintf(buffer3, "Tracked :          ");
    for (int i=0; i<32; i++) {
       if (flagsync[i]==1) {
         sprintf(str1, "%2d  ", prn[i]);
         strcat(buffer3, str1);
       }
    }
    SDL_Surface* textSurface3 =
        TTF_RenderText_Solid(font, buffer3, textColor4);
    SDL_Texture* textTexture3 =
          SDL_CreateTextureFromSurface(renderer, textSurface3);
    SDL_QueryTexture(textTexture3, NULL, NULL, &textRect3.w, &textRect3.h);
    SDL_RenderCopy(renderer, textTexture3, NULL, &textRect3);
    SDL_DestroyTexture(textTexture3);
    SDL_FreeSurface(textSurface3);

    // Update nav decoded PRNs
    // Update acquired PRNs
    sprintf(buffer4, "NavDecode:    ");
    for (int i=0; i<32; i++) {
       if (flagdec[i]==1) {
         sprintf(str1, "%2d  ", prn[i]);
         strcat(buffer4, str1);
       }
    }
    SDL_Surface* textSurface4 =
        TTF_RenderText_Solid(font, buffer4, textColor4);
    SDL_Texture* textTexture4 =
          SDL_CreateTextureFromSurface(renderer, textSurface4);
    SDL_QueryTexture(textTexture4, NULL, NULL, &textRect4.w, &textRect4.h);
    SDL_RenderCopy(renderer, textTexture4, NULL, &textRect4);
    SDL_DestroyTexture(textTexture4);
    SDL_FreeSurface(textSurface4);

    // Update LLA values
    mlock(hobsvecmtx);
    sprintf(buffer5,
            "Lat: %.7f   Lon: %.7f   Alt: %.1f   GDOP: %.1f   CB: %.5e   SVs: %02d",
            lat, lon, hgt, gdop, clkBias/CTIME, nsat);
    unmlock(hobsvecmtx);
    SDL_Surface* textSurface5 =
       TTF_RenderText_Solid(font, buffer5, textColor3);
    SDL_Texture* textTexture5 =
       SDL_CreateTextureFromSurface(renderer, textSurface5);
    SDL_QueryTexture(textTexture5, NULL, NULL, &textRect5.w, &textRect5.h);
    SDL_RenderCopy(renderer, textTexture5, NULL, &textRect5);
    SDL_DestroyTexture(textTexture5);
    SDL_FreeSurface(textSurface5);

    // List info for valid measurements
    for (int i=0; i<nsat; i++) {
      if (i==0) {
        sprintf(buffer6,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface6 =
          TTF_RenderText_Solid(font, buffer6, textColor1);
        SDL_Texture* textTexture6 =
          SDL_CreateTextureFromSurface(renderer, textSurface6);
        SDL_QueryTexture(textTexture6, NULL, NULL, &textRect6.w, &textRect6.h);
        SDL_RenderCopy(renderer, textTexture6, NULL, &textRect6);
        SDL_DestroyTexture(textTexture6);
        SDL_FreeSurface(textSurface6);
      }
      if (i==1) {
        sprintf(buffer7,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface7 =
          TTF_RenderText_Solid(font, buffer7, textColor1);
        SDL_Texture* textTexture7 =
          SDL_CreateTextureFromSurface(renderer, textSurface7);
        SDL_QueryTexture(textTexture7, NULL, NULL, &textRect7.w, &textRect7.h);
        SDL_RenderCopy(renderer, textTexture7, NULL, &textRect7);
        SDL_DestroyTexture(textTexture7);
        SDL_FreeSurface(textSurface7);
      }
      if (i==2) {
        sprintf(buffer8,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface8 =
          TTF_RenderText_Solid(font, buffer8, textColor1);
        SDL_Texture* textTexture8 =
          SDL_CreateTextureFromSurface(renderer, textSurface8);
        SDL_QueryTexture(textTexture8, NULL, NULL, &textRect8.w, &textRect8.h);
        SDL_RenderCopy(renderer, textTexture8, NULL, &textRect8);
        SDL_DestroyTexture(textTexture8);
        SDL_FreeSurface(textSurface8);
      }
      if (i==3) {
        sprintf(buffer9,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface9 =
          TTF_RenderText_Solid(font, buffer9, textColor1);
        SDL_Texture* textTexture9 =
          SDL_CreateTextureFromSurface(renderer, textSurface9);
        SDL_QueryTexture(textTexture9, NULL, NULL, &textRect9.w, &textRect9.h);
        SDL_RenderCopy(renderer, textTexture9, NULL, &textRect9);
        SDL_DestroyTexture(textTexture9);
        SDL_FreeSurface(textSurface9);
      }
      if (i==4) {
        sprintf(buffer10,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface10 =
          TTF_RenderText_Solid(font, buffer10, textColor1);
        SDL_Texture* textTexture10 =
          SDL_CreateTextureFromSurface(renderer, textSurface10);
        SDL_QueryTexture(textTexture10, NULL, NULL, &textRect10.w, &textRect10.h);
        SDL_RenderCopy(renderer, textTexture10, NULL, &textRect10);
        SDL_DestroyTexture(textTexture10);
        SDL_FreeSurface(textSurface10);
      }
      if (i==5) {
        sprintf(buffer11,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface11 =
          TTF_RenderText_Solid(font, buffer11, textColor1);
        SDL_Texture* textTexture11 =
          SDL_CreateTextureFromSurface(renderer, textSurface11);
        SDL_QueryTexture(textTexture11, NULL, NULL, &textRect11.w, &textRect11.h);
        SDL_RenderCopy(renderer, textTexture11, NULL, &textRect11);
        SDL_DestroyTexture(textTexture11);
        SDL_FreeSurface(textSurface11);
      }
      if (i==6) {
        sprintf(buffer12,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface12 =
          TTF_RenderText_Solid(font, buffer12, textColor1);
        SDL_Texture* textTexture12 =
          SDL_CreateTextureFromSurface(renderer, textSurface12);
        SDL_QueryTexture(textTexture12, NULL, NULL, &textRect12.w, &textRect12.h);
        SDL_RenderCopy(renderer, textTexture12, NULL, &textRect12);
        SDL_DestroyTexture(textTexture12);
        SDL_FreeSurface(textSurface12);
      }
      if (i==7) {
        sprintf(buffer13,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface13 =
          TTF_RenderText_Solid(font, buffer13, textColor1);
        SDL_Texture* textTexture13 =
          SDL_CreateTextureFromSurface(renderer, textSurface13);
        SDL_QueryTexture(textTexture13, NULL, NULL, &textRect13.w, &textRect13.h);
        SDL_RenderCopy(renderer, textTexture13, NULL, &textRect13);
        SDL_DestroyTexture(textTexture13);
        SDL_FreeSurface(textSurface13);
      }
      if (i==8) {
        sprintf(buffer14,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface14 =
          TTF_RenderText_Solid(font, buffer14, textColor1);
        SDL_Texture* textTexture14 =
          SDL_CreateTextureFromSurface(renderer, textSurface14);
        SDL_QueryTexture(textTexture14, NULL, NULL, &textRect14.w, &textRect14.h);
        SDL_RenderCopy(renderer, textTexture14, NULL, &textRect14);
        SDL_DestroyTexture(textTexture14);
        SDL_FreeSurface(textSurface14);
      }
      if (i==9) {
        sprintf(buffer15,"G%02.0f  TOW=%.1f  Week=%.0f  SNR=%.1f  PR=%.1f  Az=%05.1f  El=%04.1f  vk1=%7.1f",
          obs_v[i*10+0],
          obs_v[i*10+5],
          obs_v[i*10+6],
          obs_v[i*10+7],
          obs_v[i*10+4],
          obs_v[i*10+8],
          obs_v[i*10+9],
          vk1_v[i]);
        SDL_Surface* textSurface15 =
          TTF_RenderText_Solid(font, buffer15, textColor1);
        SDL_Texture* textTexture15 =
          SDL_CreateTextureFromSurface(renderer, textSurface15);
        SDL_QueryTexture(textTexture15, NULL, NULL, &textRect15.w, &textRect15.h);
        SDL_RenderCopy(renderer, textTexture15, NULL, &textRect15);
        SDL_DestroyTexture(textTexture15);
        SDL_FreeSurface(textSurface15);
      }
    } // end for

    // Render all
    SDL_RenderPresent(renderer);

    // Add delay
    SDL_Delay(GUI_RATE);

    // Update elapsed time for next iteration
    elapsedTime = elapsedTime + 0.2;

  } // end of main while loop

  // Close out and free memory (May not need since done in while loop)
  SDL_DestroyTexture(textTexture1);
  SDL_DestroyTexture(textTexture2);
  SDL_DestroyTexture(textTexture3);
  SDL_DestroyTexture(textTexture4);
  SDL_DestroyTexture(textTexture5);
  SDL_DestroyTexture(textTexture6);
  SDL_DestroyTexture(textTexture7);
  SDL_DestroyTexture(textTexture8);
  SDL_DestroyTexture(textTexture9);
  SDL_DestroyTexture(textTexture10);
  SDL_DestroyTexture(textTexture11);
  SDL_DestroyTexture(textTexture12);
  SDL_DestroyTexture(textTexture13);
  SDL_DestroyTexture(textTexture14);
  SDL_DestroyTexture(textTexture15);
  SDL_DestroyTexture(textTexture16);

  SDL_FreeSurface(textSurface1);
  SDL_FreeSurface(textSurface2);
  SDL_FreeSurface(textSurface3);
  SDL_FreeSurface(textSurface4);
  SDL_FreeSurface(textSurface5);
  SDL_FreeSurface(textSurface6);
  SDL_FreeSurface(textSurface7);
  SDL_FreeSurface(textSurface8);
  SDL_FreeSurface(textSurface9);
  SDL_FreeSurface(textSurface10);
  SDL_FreeSurface(textSurface11);
  SDL_FreeSurface(textSurface12);
  SDL_FreeSurface(textSurface13);
  SDL_FreeSurface(textSurface14);
  SDL_FreeSurface(textSurface15);
  SDL_FreeSurface(textSurface16);

  // Final wrap-up
  TTF_CloseFont(font);
  TTF_Quit();
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();

  // Announce thread closed
  printf("GUI thread closed.\n");

  // Main return
  return 0;
}
