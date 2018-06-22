#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "bcl.h"
#include "conversion.h"
#include "BestMatch.h"
#include "patchMatch.h"
#include "pyramid.h"


int L = 1;
int neighborSize = 5;
int sum = 0;
int colorTransfert = 0;
int bestmatch = 1;
int patchmatch = 0;
int nbMaxIter = 3;
float w = 1.;
float levelweight = 1.;
float K = 5.;

void
summarize(pnm A, pnm Aprim, pnm B, pnm Bprim)
{
  int Aw = pnm_get_width(A);
  int Ah = pnm_get_height(A);
  int Bw = pnm_get_width(B);
  int Bh = pnm_get_height(B);

  int width = max(Aw, Bw) * 2 + 10;
  int height = Ah + Bh + 10;

  pnm out = pnm_new(width, height, PnmRawPpm);

  for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
      for(int k=0; k<3; k++)
        pnm_set_component(out, i, j, k, 255/2);
    }
  }

  unsigned short val;
  //A
  for(int i=0; i<Ah; i++){
    for(int j=0; j<Aw; j++){
      for(int k=0; k<3; k++){
        val = pnm_get_component(A, i, j, k);
        pnm_set_component(out, i, j, k, val);
      }
    }
  }
  //A'
  for(int i=0; i<Ah; i++){
    for(int j=Aw+10; j<2*Aw+10; j++){
      for(int k=0; k<3; k++){
        val = pnm_get_component(Aprim, i, j-Aw-10, k);
        pnm_set_component(out, i, j, k, val);
      }
    }
  }
  //B
  for(int i=Ah+10; i<Ah+Bh+10; i++){
    for(int j=0; j<Bw; j++){
      for(int k=0; k<3; k++){
        val = pnm_get_component(B, i-Ah-10, j, k);
        pnm_set_component(out, i, j, k, val);
      }
    }
  }
  //B'
  for(int i=Ah+10; i<Ah+Bh+10; i++){
    for(int j=Bw+10; j<2*Bw+10; j++){
      for(int k=0; k<3; k++){
        val = pnm_get_component(Bprim, i-Ah-10, j-Bw-10, k);
        pnm_set_component(out, i, j, k, val);
      }
    }
  }

  pnm_save(out, PnmRawPpm, "summarize.ppm");
  pnm_free(out);
}

void
process(char* A_name, char* Aprim_name, char* B_name)
{
  time_t debut = time(NULL);

  printf("Chargement Image\n");
  pnm ims = pnm_load(A_name);
  pnm ims_filtre = pnm_load(Aprim_name);
  pnm imt = pnm_load(B_name);

  if(ims == NULL || ims_filtre == NULL || imt == NULL){
    fprintf(stderr,"Erreur chargement images\n");
    exit(EXIT_FAILURE);
  }

  int srcW = pnm_get_width(ims);
  int srcH = pnm_get_height(ims);
  int dstW = pnm_get_width(imt);
  int dstH = pnm_get_height(imt);

  pnm imt_filtre = pnm_new(dstW, dstH, PnmRawPpm);

  unsigned short* sRGB = pnm_get_image(ims);
  unsigned short* sfRGB = pnm_get_image(ims_filtre);
  unsigned short* tRGB = pnm_get_image(imt);
  unsigned short* data_out = malloc(sizeof(unsigned short) * 3 * dstW * dstH);

  printf("Conversion RGB -> YIQ\n");

  float* tmp_s = convertUnsignedShort2Float(sRGB, 3 * srcW, srcH);
  float* tmp_sf = convertUnsignedShort2Float(sfRGB, 3 * srcW, srcH);
  float* tmp_t = convertUnsignedShort2Float(tRGB, 3 * dstW, dstH);

  float* s_yiq = convertRGB2YIQ(srcW, srcH, tmp_s);
  float* sf_yiq = convertRGB2YIQ(srcW, srcH, tmp_sf);
  float* t_yiq = convertRGB2YIQ(dstW, dstH, tmp_t);

  float* data_source = channel(srcW, srcH, s_yiq, 0);
  float* data_source_filter = channel(srcW, srcH, sf_yiq, 0);
  float* data_target = channel(dstW, dstH, t_yiq, 0);

  printf("Création des pyramides A A\' B B\' \n");

  Pyramid* A = init_pyramid(srcW, srcH, 0, L);
  Pyramid* Aprim = init_pyramid(srcW, srcH, 0, L);
  Pyramid* B = init_pyramid(dstW, dstH, 0, L);
  Pyramid* Bprim = init_pyramid(dstW, dstH, 1, L);

  printf("Début Algo\n");

  int pixel, p, q;
  float val;
  printf("\n");
  float** kernel = initKernel(neighborSize);

  for(int l=0; l<L; l++){
    printf("Level %d\n", l);

    srcW = A[l].cols;
    srcH = A[l].rows;
    dstW = B[l].cols;
    dstH = B[l].rows;

    free(A[l].lum);
    free(Aprim[l].lum);
    free(B[l].lum);

    if( l == L-1){
      A[l].lum = data_source;
      Aprim[l].lum = data_source_filter;
      B[l].lum = data_target;
    }
    else{
      A[l].lum = pyramid_level(A[L-1].cols, A[L-1].rows, data_source, L-1-l);
      Aprim[l].lum = pyramid_level(A[L-1].cols, A[L-1].rows, data_source_filter, L-1-l);
      B[l].lum = pyramid_level(B[L-1].cols, B[L-1].rows, data_target, L-1-l);
    }

    for(int i=0; i<dstH; i++){
      for(int j=0; j<dstW; j++){
        pixel = rand() % (srcH * srcW);
        val = rand() % 256;
        Bprim[l].lum[i*dstW+j] = val;
        Bprim[l].s[i*dstW+j] = pixel;
      }
    }

    if( bestmatch ){
      for(int i=0; i<dstH; i++){

        printf("Ligne %d / %d\r", i, dstH-1);
        fflush(stdout);

        for(int j=0; j<dstW; j++){
          q = i*dstW+j;

          p = BestMatch(A, Aprim, i, j, B, Bprim, l, neighborSize, w, levelweight, L, K, kernel);

          /*
          if(w == 1)
            Bprim[l].lum[q] = B[l].lum[q] + (Aprim[l].lum[p] - A[l].lum[p]);
          else{
          */

            Bprim[l].lum[q] = Aprim[l].lum[p];
            Bprim[l].s[q] = p;

          /*
          }
          */
        }
      }
    }

    if( patchmatch ){
      int istart,iend, jstart, jend, step, p_patchmatch, p_coh, x, y;
      float dist_coh, dist_patch;
      for(int k=0; k<nbMaxIter; k++){
        printf("Itération %d:\n",k);

        if(k%2==0){
          istart = 0; iend = dstH; jstart = 0; jend = dstW; step = 1;
        }else{
          istart = dstH-1; iend = -1; jstart = dstW-1; jend = -1; step = -1;
        }
        for(int i=istart; i!=iend; i+=step){

          if(k%2==0)
            printf("Ligne %d / %d\r", i, dstH-1);
          else
            printf("Ligne %d / %d\r", dstH-1-i, dstH-1);
          fflush(stdout);

          for(int j=jstart; j!=jend; j+=step){
            q = i*dstW+j;

            p_patchmatch = patchMatch(A, Aprim, i, j, B, Bprim, l, neighborSize, w, levelweight, kernel, Bprim[l].s, k);
            x = p_patchmatch/A[l].cols;
            y = p_patchmatch - x * A[l].cols;
            dist_patch = dist(x, y, A, Aprim, i, j, B, Bprim, l, neighborSize, w, levelweight, kernel);

            p_coh = Ashikhmin(A, Aprim, i, j, B, Bprim, l, neighborSize, w, levelweight, kernel);
            x = p_coh/A[l].cols;
            y = p_coh - x * A[l].cols;
            dist_coh = dist(x, y, A, Aprim, i, j, B, Bprim, l, neighborSize, w, levelweight, kernel);

            if(dist_coh < dist_patch * (1 + powf(2, l-L)*K) )
              p = p_coh;
            else
              p = p_patchmatch;

            Bprim[l].lum[q] = Aprim[l].lum[p];
            Bprim[l].s[q] = p;
          }
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  freeKernel(kernel, neighborSize);

  printf("Fin algo\n");

  printf("Conversion Yiq -> RGB \n");

  dstH = Bprim[L-1].rows;
  dstW = Bprim[L-1].cols;

  float* it;
  float* qt;

  if(w != 1 || colorTransfert){

    printf("Transfert de couleur\n");

    float* isf = channel(A[L-1].cols, A[L-1].rows, sf_yiq, 1);
    float* qsf = channel(A[L-1].cols, A[L-1].rows, sf_yiq, 2);

    it = malloc(sizeof(float) * dstW * dstH);
    qt = malloc(sizeof(float) * dstW * dstH);

    for(int i=0; i<dstH; i++){
      for(int j=0; j<dstW; j++){
        pixel = Bprim[L-1].s[i*dstW+j];
        it[i*dstW+j] = isf[pixel];
        qt[i*dstW+j] = qsf[pixel];
      }
    }
  }
  else{
    it = channel(dstW, dstH, t_yiq, 1);
    qt = channel(dstW, dstH, t_yiq, 2);
  }

  float* yiq = putChannel(dstW, dstH, Bprim[L-1].lum, it, qt);

  float* rgb = convertYIQ2RGB(dstW, dstH, yiq);

  convertFloat2UnsignedShort(rgb, data_out, 3*dstW, dstH);

  printf("Sauvegarde image final\n");

  for(int i=0; i<dstH; i++){
    for(int j=0; j<dstW; j++){
      for(int k=0; k<3; k++){
	       pnm_set_component(imt_filtre, i, j, k, data_out[i*3*dstW+(j*3)+k]);
      }
    }
  }

  pnm_save(imt_filtre, PnmRawPpm, "output.ppm");

  if( sum )
    summarize(ims, ims_filtre, imt, imt_filtre);

  printf("Liberation de la mémoire\n");

  free(rgb);
  free(yiq);
  free(qt);
  free(it);
  free_pyramid(A, 0, L);
  free_pyramid(Aprim, 0, L);
  free_pyramid(B, 0, L);
  free_pyramid(Bprim, 1, L);
  free(t_yiq);
  free(sf_yiq);
  free(s_yiq);
  free(tmp_t);
  free(tmp_sf);
  free(tmp_s);
  free(tRGB);
  free(sfRGB);
  free(sRGB);
  free(ims);
  free(ims_filtre);
  free(imt);
  pnm_free(imt_filtre);


  time_t fin = time(NULL);
  printf("Temps écoulé: %f \n", difftime(fin,debut));

  printf("\nFIN\n");
}

void
Usage(char* command)
{
  printf("Usage : \n");
  printf("\t%s <A> <A'> <B> [options]\n",command);
  printf("\tOptions:\n");
  printf("\t\t-h/--help : Affiche l'aide\n");
  printf("\t\t-w: Poid distance métrique sur ||A-B|| compris [0,1]\n");
  printf("\t\t\tw = 0 => Synthèse, w = 1 => filtre, w compris ]0,1[ => transfert\n");
  printf("\t\t-lw/--levelWeight: Poid entre chaque niveau de la pyramide gaussienne\n");
  printf("\t\t-ns/--neighborSize: Taille du patch au niveau l\n");
  printf("\t\t-l/--levelPyramid: Taille de la pyramide gaussienne\n");
  printf("\t\t-k/--keppa : poid favorisant la coherence\n");
  printf("\t\t-sum/--sumarize: créer une image mosaïque de A A\' B B\'\n");
  printf("\t\t-ct/--colortransfert: force le transfert de couleur de A\' vers B\'\n");
  printf("\t\t-pm/--patchMatch: Utilisation de patchMatch et Ashikhmin comme méthode de recherche\n");
  printf("\t\t-bm/--bestMatch: Utilisation de Wei Levoy et de Ashikhmin comme méthode de recherche (méthode utilisé par défaut)\n");
  printf("\t\t-nb/--maxIteration: Le nombre d'itération pour patchMatch\n");
  printf("\t\t\tAttention maxIteration doit être utilisé que si option -pm avant\n");

  exit(EXIT_FAILURE);
}

int
main(int argc, char** argv)
{
  if(argc < 4) Usage(argv[0]);

  if(argc > 4){
    for(int k=4; k<argc; k++){

      if(!strcmp(argv[k],"-w")){
        k++;
        w = atof(argv[k]);
      }

      else if(!strcmp(argv[k], "-lw") || !strcmp(argv[k], "--levelWeight")){
        k++;
        levelweight = atof(argv[k]);
      }

      else if(!strcmp(argv[k], "-ns") || !strcmp(argv[k], "--neighborSize")){
        k++;
        neighborSize = atoi(argv[k]);
      }

      else if(!strcmp(argv[k], "-l") || !strcmp(argv[k], "--levelPyramid")){
        k++;
        L = atoi(argv[k]);
      }

      else if(!strcmp(argv[k], "-k") || !strcmp(argv[k], "--keppa")){
        k++;
        K = atof(argv[k]);
      }

      else if(!strcmp(argv[k], "-sum") || !strcmp(argv[k], "--summarize")){
        sum = 1;
      }

      else if(!strcmp(argv[k], "-ct") || !strcmp(argv[k], "--colortransfert")){
        colorTransfert = 1;
      }

      else if(!strcmp(argv[k], "-pm") || !strcmp(argv[k], "--patchMatch")){
        patchmatch = 1;
        bestmatch = 0;
      }

      else if(!strcmp(argv[k], "-bm") || !strcmp(argv[k], "--bestmatch")){
        bestmatch = 1;
        patchmatch = 0;
      }

      else if(!strcmp(argv[k], "-nb") || !strcmp(argv[k], "--maxIteration")){
        if( patchmatch ){
          k++;
          nbMaxIter = atoi(argv[k]);
        }else{
          fprintf(stderr,"Attention, vous devez choisir la méthode patchMatch pour utiliser cette option\n");
          exit(EXIT_FAILURE);
        }
      }

      else{ Usage(argv[0]); }
    }
  }

  process(argv[1], argv[2], argv[3]);

  return EXIT_SUCCESS;
}
