#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
#include <time.h>

typedef struct {
  bool *site;
  int *cluster;
  int *cluster_sizes;
  int stride;
  int sx;
  int sy;
  int c_n;
  int c_c;
} sitegrid;

typedef struct {
  bool *bondh;
  bool *bondv;
  int *cluster;
  int *cluster_sizes;
  int stride;
  int sx;
  int sy;
  int c_n;
  int c_c;
} bondgrid;

typedef struct {
  int largest_cluster;
  bool percolates;
} gridstats;

sitegrid seed_sitegrid(int size, float p)
{
  srand(time(NULL));
  sitegrid sg;
  sg.site = malloc(size * size * sizeof(bool));
  sg.cluster = calloc(size * size, sizeof(int));
  sg.cluster_sizes = calloc(size * size, sizeof(int));
  sg.sx = size;
  sg.sy = size;
  sg.stride = size;
  sg.c_n = 1;
  sg.c_c = 0;
  for (int i=0; i<size*size; i++) {
    sg.site[i] = ((float)rand() / (float)RAND_MAX) < p;
  }
  return sg;
}

void reset_sitegrid(sitegrid sg)
{
  for (int x=0; x<sg.sx; x++) {
    for (int y=0; y<sg.sy; y++) {
      int ind = (y * sg.stride) + x;
      sg.cluster[ind] = 0;
    }
  }
  free(sg.cluster_sizes);
  sg.cluster_sizes = calloc(sg.sx * sg.sy, sizeof(int));
}

void free_sitegrid(sitegrid sg)
{
  free(sg.site);
  free(sg.cluster);
  free(sg.cluster_sizes);
}

void print_sitegrid(sitegrid sg, bool clusters)
{
  for (int y=0; y<sg.sy; y++) {
    for (int x=0; x<sg.sx; x++) {
      if (clusters) {
        printf("%02d ", sg.cluster[x + y*sg.stride]);
      } else {
        printf("%s", sg.site[x + y*sg.stride] ? "O" : " ");
      }
    }
    printf("\n");
  }
}

int sitegrid_dfs(sitegrid sg)
{
  int cluster = sg.c_c + 1;
  int *stack = malloc(2 * sg.sx * sg.sy * sizeof(int));
  for (int x=0; x<sg.sx; x++) {
    for (int y=0; y<sg.sy; y++) {
      int ind = (y * sg.stride) + x;
      if (!sg.site[ind] || sg.cluster[ind] != 0) continue;
      int sp = 0;
      sg.cluster[ind] = cluster;
      stack[sp++] = x;
      stack[sp++] = y;
      int clust_sz = 1;
      while (sp > 0) {
        int cy = stack[--sp];
        int cx = stack[--sp];
        if (cx > 0) {
          int left = (cy * sg.stride) + cx - 1;
          if (sg.site[left] && sg.cluster[left] == 0) {
            sg.cluster[left] = cluster;
            clust_sz++;
            stack[sp++] = cx - 1;
            stack[sp++] = cy;
          }
        }
        if (cx < sg.sx - 1) {
          int right = (cy * sg.stride) + cx + 1;
          if (sg.site[right] && sg.cluster[right] == 0) {
            sg.cluster[right] = cluster;
            clust_sz++;
            stack[sp++] = cx + 1;
            stack[sp++] = cy;
          }
        }
        if (cy > 0) {
          int top = ((cy - 1) * sg.stride) + cx;
          if (sg.site[top] && sg.cluster[top] == 0) {
            sg.cluster[top] = cluster;
            clust_sz++;
            stack[sp++] = cx;
            stack[sp++] = cy - 1;
          }
        }
        if (cy < sg.sy - 1) {
          int bottom = ((cy + 1) * sg.stride) + cx;
          if (sg.site[bottom] && sg.cluster[bottom] == 0) {
            sg.cluster[bottom] = cluster;
            clust_sz++;
            stack[sp++] = cx;
            stack[sp++] = cy + 1;
          }
        }
      }
      sg.cluster_sizes[cluster - 1] = clust_sz;
      cluster += sg.c_n;
    }
  }
  free(stack);
  return cluster - sg.c_n;
}

sitegrid sitegrid_divide(sitegrid sg, int x, int y, int sx, int sy, int c_c, int c_n)
{
  sitegrid new;
  int offset = (y * sg.stride) + x;
  new.site = sg.site + offset;
  new.cluster = sg.cluster + offset;
  new.cluster_sizes = sg.cluster_sizes;
  new.stride = sg.stride;
  new.sx = sx;
  new.sy = sy;
  new.c_n = sg.c_n * c_n;
  new.c_c = sg.c_c + (sg.c_n * c_c);
  return new;
}

void sitegrid_merge_clusters(sitegrid sg, int *alias, int max, int pos, bool vert)
{
  int adj = vert ? 1 : sg.stride;
  int jump = vert ? sg.stride : 1;
  int len = vert ? sg.sy : sg.sx;
  for (int i=0; i<len; i++) {
    int ind = vert ? (i * sg.stride) + pos : (pos * sg.stride) + i;
    int ind2 = ind + adj;
    if (!sg.site[ind] || !sg.site[ind2]) continue;
    int c1 = alias[sg.cluster[ind] - 1], c2 = alias[sg.cluster[ind2] - 1];
    if (c1 == c2) continue;
    if (c1 > c2) {
      int t = c2;
      c2 = c1;
      c1 = t;
    }
    for (int i=c2 - 1; i<max; i++) {
      if (alias[i] == c2) {
        alias[i] = c1;
      }
    }
    sg.cluster_sizes[c1 - 1] += sg.cluster_sizes[c2 - 1];
    sg.cluster_sizes[c2 - 1] = 0;
  }
}

bool sitegrid_percolates(sitegrid sg, int *alias)
{
  for (int x=0; x<sg.sx; x++) {
    int ind = x;
    int ind2 = ((sg.sy - 1) * sg.stride) + x;
    if (!sg.site[ind] || !sg.site[ind2]) continue;
    if (alias[sg.cluster[ind] - 1] == alias[sg.cluster[ind2] - 1]) {
      return true;
    }
  }
  for (int y=0; y<sg.sy; y++) {
    int ind = y * sg.stride;
    int ind2 = (y * sg.stride) + (sg.sx - 1);
    if (!sg.site[ind] || !sg.site[ind2]) continue;
    if (alias[sg.cluster[ind] - 1] == alias[sg.cluster[ind2] - 1]) {
      return true;
    }
  }
  return false;
}

gridstats sitegrid_test(sitegrid sg, char speedup)
{ 
  int max = 0;
  int subsize = sg.sx >> speedup;
  for (int i = 0; i < 1<<(2*speedup); i++) {
    int x = i & ((1 << speedup) - 1);
    int y = (i & (((1 << speedup) - 1) << speedup)) >> speedup;
    sitegrid sub = sitegrid_divide(sg, subsize * x, subsize * y, subsize, subsize, i, 1 << (2*speedup));
    int m = sitegrid_dfs(sub);
    if (m > max) max = m;
  }
  int *alias = malloc(max * sizeof(int));
  for (int i=0; i<max; i++) alias[i] = i + 1;
  for (int i = 1; i < 1<<speedup; i++) {
    sitegrid_merge_clusters(sg, alias, max, (i * subsize) - 1, false);
    sitegrid_merge_clusters(sg, alias, max, (i * subsize) - 1, true);
  }
  gridstats out;
  out.largest_cluster = 0;
  for (int i=0; i<max; i++) {
    if (sg.cluster_sizes[i] > out.largest_cluster) {
      out.largest_cluster = sg.cluster_sizes[i];
    }
  }
  out.percolates = sitegrid_percolates(sg, alias);
  return out;
}

int main(int argc, char *argv[])
{
  int size = atoi(argv[1]);
  float p = atof(argv[2]);
  printf("Testing %i x %i grid with p = %f\n\n", size, size, p);
  sitegrid sg = seed_sitegrid(size, p);
  gridstats stats = sitegrid_test(sg, 0);
  printf(
    "Largest cluster: %i\nPercolates: %s\n\n",
    stats.largest_cluster,
    stats.percolates ? "Yes" : "No"
  );
}