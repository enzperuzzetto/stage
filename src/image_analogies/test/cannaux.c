#include <stdio.h>
#include <stdlib.h>

#include <bcl.h>

void
process(char* ims_name)
{
  pnm ims = pnm_load(ims_name);
  int cols_s = pnm_get_width(ims);
  int rows_s = pnm_get_height(ims);

  pnm r = pnm_new(cols_s, rows_s, PnmRawPpm);
  pnm g = pnm_new(cols_s, rows_s, PnmRawPpm);
  pnm b = pnm_new(cols_s, rows_s, PnmRawPpm);

  unsigned short* data = pnm_get_channel(ims, NULL, 0);
  pnm_set_channel(r, data, 0);

  pnm_save(r, PnmRawPpm, "rouge.ppm");  

  data = pnm_get_channel(ims, NULL, 1);
  pnm_set_channel(g, data, 1);

  pnm_save(g, PnmRawPpm, "vert.ppm");

  data = pnm_get_channel(ims, NULL, 2);
  pnm_set_channel(b, data, 2);

  pnm_save(b, PnmRawPpm, "bleu.ppm");

  free(data);
  free(ims);
  pnm_free(r);
  pnm_free(g);
  pnm_free(b);
}


void
usage (char *s){
  fprintf(stderr, "Usage: %s <ims> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 1
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
  process(argv[1]);
  return EXIT_SUCCESS;
}
