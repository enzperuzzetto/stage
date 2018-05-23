#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "bcl.h"
#include "BestMatch.h"
#include "conversion.h"


void
process(char* ims_name, char* ims_filter_name, char* imt_name)
{
  pnm imsRGB = pnm_load(ims_name);
  pnm ims_filterRGB = pnm_load(ims_filter_name);
  int cols_s = pnm_get_width(imsRGB);
  int rows_s = pnm_get_height(imsRGB);

  pnm imtRGB = pnm_load(imt_name);
  int cols_t = pnm_get_width(imtRGB);
  int rows_t = pnm_get_height(imtRGB);

  pnm imt_filter = pnm_new(cols_t, rows_t, PnmRawPpm);

  
  //*********** covert RGB to YIQ ***************//
  
  printf("Convertion en YIQ\n");
  unsigned short* sRGB = pnm_get_image(imsRGB);
  unsigned short* s_fRGB = pnm_get_image(ims_filterRGB);
  unsigned short* tRGB = pnm_get_image(imtRGB);

  float* tmp_s = convertUnsignedShort2Float(sRGB, 3 * cols_s, rows_s);
  float* tmp_sf = convertUnsignedShort2Float(s_fRGB, 3 * cols_s, rows_s);
  float* tmp_t = convertUnsignedShort2Float(tRGB, 3 * cols_t, rows_t);
  
  float* s_yiq = convertRGB2YIQ(cols_s, rows_s, tmp_s);
  float* sf_yiq = convertRGB2YIQ(cols_s, rows_s, tmp_sf);
  float* t_yiq = convertRGB2YIQ(cols_t, rows_t, tmp_t);

  unsigned short* data_out = malloc(sizeof(unsigned short) * 3 * cols_t * rows_t);
 
  float* data_source = channel(cols_s, rows_s, s_yiq, 0);
  float* data_source_filter = channel(cols_s, rows_s, sf_yiq, 0);
  float* data_target = channel(cols_t, rows_t, t_yiq, 0);
  float* data_target_filter = malloc(sizeof(float) * cols_t * rows_t);

  printf("fait\n");
   
  //*********** Compute L level of pyramid *****************************//
  printf("Démarrage algo \n");
    
  Pyramid* source = init_pyramid(cols_s, rows_s, 0);
  Pyramid* source_filter = init_pyramid(cols_s, rows_s, 0);
  Pyramid* target = init_pyramid(cols_t, rows_t, 0);
  Pyramid* target_filter = init_pyramid(cols_t, rows_t, 1);

    
  for(int l=0; l<L; l++){
    printf("Début level %d\n",l);
     
    if( l == L-1) {
      free(source[l].data);
      free(source_filter[l].data);
      free(target[l].data);
      
      source[l].data = data_source;
      source_filter[l].data = data_source_filter;
      target[l].data = data_target;

    }else{
      free(source[l].data);
      free(source_filter[l].data);
      free(target[l].data);
      free(target_filter[l].data);
      
      source[l].data = pyramid_level(cols_s, rows_s, data_source, L-1-l);
      source_filter[l].data = pyramid_level(cols_s, rows_s, data_source_filter, L-1-l);
           
      target[l].data = pyramid_level(cols_t, rows_t, data_target, L-1-l);
      target_filter[l].data = pyramid_level(cols_t, rows_t, data_target_filter, L-1-l);
    }
    //************* Compute Features Means & standard deviation **************************
    
    mean_img(source[l].cols, source[l].rows, source[l].data, source[l].mean, NEIGHBOOR_SIZE_FINER, 0);
    mean_img(source_filter[l].cols, source_filter[l].rows, source_filter[l].data, source_filter[l].mean, NEIGHBOOR_SIZE_FINER, 1);
    mean_img(target[l].cols, target[l].rows, target[l].data, target[l].mean, NEIGHBOOR_SIZE_FINER, 0);      
    
    sd_img(source[l].cols, source[l].rows, source[l].data, source[l].mean, source[l].sd, NEIGHBOOR_SIZE_FINER, 0);
    sd_img(source_filter[l].cols, source_filter[l].rows, source_filter[l].data, source_filter[l].mean, source_filter[l].sd, NEIGHBOOR_SIZE_FINER, 1);
    sd_img(target[l].cols, target[l].rows, target[l].data, target[l].mean, target[l].sd, NEIGHBOOR_SIZE_FINER, 0);

      
    //*********** Compute feature l-1 with neighboor size coarser ***********************// 
    if( l > 0){
      mean_img(source[l-1].cols, source[l-1].rows, source[l-1].data, source[l-1].mean, NEIGHBOOR_SIZE_COARSER, 0);
      mean_img(source_filter[l-1].cols, source_filter[l-1].rows, source_filter[l-1].data, source_filter[l-1].mean, NEIGHBOOR_SIZE_COARSER, 0);
      mean_img(target[l-1].cols, target[l-1].rows, target[l-1].data, target[l-1].mean, NEIGHBOOR_SIZE_COARSER, 0);
      mean_img(target_filter[l-1].cols, target_filter[l-1].rows, target_filter[l-1].data, target_filter[l-1].mean, NEIGHBOOR_SIZE_COARSER, 0);
	
      sd_img(source[l-1].cols, source[l-1].rows, source[l-1].data, source[l-1].mean, source[l-1].sd, NEIGHBOOR_SIZE_COARSER, 0);
      sd_img(source_filter[l-1].cols, source_filter[l-1].rows, source_filter[l-1].data, source_filter[l-1].mean, source_filter[l-1].sd, NEIGHBOOR_SIZE_COARSER, 0);
      sd_img(target[l-1].cols, target[l-1].rows, target[l-1].data, target[l-1].mean, target[l-1].sd, NEIGHBOOR_SIZE_COARSER, 0);
      sd_img(target_filter[l-1].cols, target_filter[l-1].rows, target_filter[l-1].data, target_filter[l-1].mean, target_filter[l-1].sd, NEIGHBOOR_SIZE_COARSER, 0);
    }
  
    //************* Best Match   **********************************************************

    int p, q;

    for(int i=0; i<target_filter[l].rows; i++){
      for(int j=0; j<target_filter[l].cols; j++){

	q = i*target_filter[l].cols+j;
	        
	target_filter[l].mean[q] = mean(target_filter[l].cols, target_filter[l].rows, target_filter[l].data, q, NEIGHBOOR_SIZE_FINER, 1);
	target_filter[l].sd[q] = standard_deviation(target_filter[l].cols, target_filter[l].rows, target_filter[l].data, target_filter[l].mean[q], q, NEIGHBOOR_SIZE_FINER, 1);
	
	p = bestMatch(source, source_filter, target, target_filter, l, q);

	

	target_filter[l].data[q] = source_filter[l].data[p];
	target_filter[l].s[q] = p;
      }
    }
      
    printf("Fait \n");
  }

  //*************************************************************************************

  printf("Convertion en RGB\n");
       
  //************* Copy I & Q B to B' ************************************//
#if 1
  float* it = channel(cols_t, rows_t, t_yiq, 1);
  float* qt = channel(cols_t, rows_t, t_yiq, 2);
#endif
  
#if 0 //transfer de couleur
  float* isf = channel(cols_s, rows_s, sf_yiq, 1);
  float* qsf = channel(cols_s, rows_s, sf_yiq, 2);

  float* it = malloc(sizeof(float) * cols_t * rows_t);
  float* qt = malloc(sizeof(float) * cols_t * rows_t);

  int pixel = 0;
  for(int i=0; i<rows_t; i++){
    for(int j=0; j<cols_t; j++){
      pixel = target_filter[L-1].s[i*cols_t+j];
      it[i*cols_t+j] = isf[pixel];
      qt[i*cols_t+j] = qsf[pixel];
    }
  }
#endif
  
  float* yiq = putChannel(cols_t, rows_t, target_filter[L-1].data, it, qt);

  float* rgb = convertYIQ2RGB(cols_t, rows_t, yiq);
  
  convertFloat2UnsignedShort(rgb, data_out, 3*cols_t, rows_t);


  free_pyramid(source, 0);
  free_pyramid(source_filter, 0);
  free_pyramid(target, 0);
  free_pyramid(target_filter, 1);

  for(int i=0; i<rows_t; i++){
    for(int j=0; j<cols_t; j++){
      for(int k=0; k<3; k++){
	pnm_set_component(imt_filter, i, j, k, data_out[i*3*cols_t+(j*3)+k]);
      }
    }
  }
  

  printf("fait\n");
 
  printf("Sauvegarde\n");
  pnm_save(imt_filter, PnmRawPpm, "filter.ppm");
  printf("Fait\n");


  free(imsRGB);
  free(ims_filterRGB);
  free(imtRGB);
  pnm_free(imt_filter);
  free(sRGB);
  free(s_fRGB);
  free(tRGB);
  free(tmp_s);
  free(tmp_sf);
  free(tmp_t);
  free(s_yiq);
  free(sf_yiq);
  free(t_yiq);
  free(data_out);
  free(data_target_filter);
  free(it);
  free(qt);
  free(yiq);
  free(rgb);
}

void
usage (char *s){
  fprintf(stderr, "Usage: %s <ims> <ims_filter> <imt> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 3
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
    
  process(argv[1], argv[2], argv[3]);
  
  return EXIT_SUCCESS;
}
