/*
 * test1.c
 *
 *  Created on: Mar 5, 2013
 *      Author: vova
 */

#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include "ccarray.h"

#define UNUSED(x)               ((void)(x))
#define MAX_HEADER_LENGTH       2048
#define MAX_INPUT_LINE_LENGTH   4096
#define PI                      M_PI

pthread_mutex_t lock; //Исключающая блокировка

int count_map =0;

/** known compression types */
typedef
enum compression_t {
    compression_unknown = -1,
    compression_none,
    compression_bzip,
    compression_gzip
} compression_t;

/** object data */
typedef
struct obj_t {
  double x, y, mux, muy;       /*< in radians */
  char * line;
} obj_t;

/** RA and DC measure units */
enum {
  radians,
  degrees,
};

struct ARG_data
{
	int dx,dy;
   int minx,miny;
   int maxx,maxy;
   int mincountmed;
   int nx,ny;
   int cw;
   int fit;
};
typedef struct ARG_data ARG_DATA;


typedef
struct matrix_stars{
		double	x;
		double	y;
		double	mux;
		double	muy;
		double	med_mux;
		double	med_muy;
		int		count;
		int		med_count;
	}matrix_stars;

static matrix_stars map[1000][1000];

/** uses fopen() if input file seems to be uncompresssed, and popen() if input file seems compressed */
static FILE * open_file( const char * fname, compression_t * compression )
{
  char cmd[PATH_MAX] = {0};
  FILE * input = NULL;

  if ( *compression == compression_unknown )
  {
    const char * suffix;

    if ( (suffix = strstr(fname, ".bz2")) && *(suffix + 4) == 0 ) {
      *compression = compression_bzip;
    }
    else if ( (suffix = strstr(fname, ".bz")) && *(suffix + 3) == 0 ) {
      *compression = compression_bzip;
    }
    else if ( (suffix = strstr(fname, ".gz")) && *(suffix + 3) == 0 ) {
      *compression = compression_gzip;
    }
    else {
      *compression = compression_none;
    }
  }

  switch ( *compression )
  {
  case compression_bzip:
    snprintf(cmd, sizeof(cmd) - 1, "bzip2 -dc '%s'", fname);
    if ( !(input = popen(cmd, "r"))) {
      fprintf(stderr, "popen('%s') fails: %s\n", cmd, strerror(errno));
    }
    break;

  case compression_gzip:
    snprintf(cmd, sizeof(cmd) - 1, "gzip -dc '%s'", fname);
    if ( !(input = popen(cmd, "r"))) {
      fprintf(stderr, "popen('%s') fails: %s\n", cmd, strerror(errno));
    }
    break;

  case compression_none:
    if ( !(input = fopen(fname, "r")) ) {
      fprintf(stderr, "fopen('%s') fails: %s\n", fname, strerror(errno));
    }
    break;

  default:
    fprintf(stderr, "BUG IN CODE: invalid compression tag=%d\n", *compression);
    break;
  }

  return input;
}

/* load lines from input file */
static int load_objects(FILE * input, int x, int y, int mux, int muy, char header[MAX_HEADER_LENGTH],
    ccarray_t * objects)
{
  char line[MAX_INPUT_LINE_LENGTH];
  int ic;
  char * pc;
  size_t size;
  size_t capacity;

  obj_t * obj;


  header[MAX_HEADER_LENGTH-1] = 0;

  if ( !fgets(header, MAX_HEADER_LENGTH, input) ) {
    fprintf(stderr,"fgets(header) fails: %d (%s)\n", errno, strerror(errno));
    return -1;
  }

  if ( header[MAX_HEADER_LENGTH - 1] != 0 ) {
    fprintf(stderr,"too long header line in this file\n");
    return -1;
  }

  /* remove trailing new line */
  header[strlen(header) - 1] = 0;


  size = ccarray_size(objects);
  capacity = ccarray_capacity(objects);

  while ( size < capacity && !feof(input) )
  {
    obj = ccarray_peek(objects, size);

    line[MAX_INPUT_LINE_LENGTH-1] = 0;
    if ( !fgets(line, MAX_INPUT_LINE_LENGTH, input) ) {
      break;
    }

    if ( line[MAX_INPUT_LINE_LENGTH - 1] != 0 ) {
      fprintf(stderr,"too long input line in this file\n");
      return -1;
    }

    /* remove trailing new line */
    line[strlen(line) - 1] = 0;



    for ( pc = line, ic = 1; ic < x; ++ic ) {
      if ( !(pc = strchr(pc + 1, '\t')) ) {
        break;
      }
    }
    if ( ic != x || sscanf(pc, " %lf", &obj->x) != 1 ) {
      continue;
    }


    for ( pc = line, ic = 1; ic < y; ++ic ) {
      if ( !(pc = strchr(pc + 1, '\t')) ) {
        break;
      }
    }
    if ( ic != y || sscanf(pc, " %lf", &obj->y) != 1 ) {
      continue;
    }

    for ( pc = line, ic = 1; ic < mux; ++ic ) {
	  if ( !(pc = strchr(pc + 1, '\t')) ) {
		break;
	  }
	}
	if ( ic != mux || sscanf(pc, " %lf", &obj->mux) != 1 ) {
	  continue;
	}

	for ( pc = line, ic = 1; ic < muy; ++ic ) {
	  if ( !(pc = strchr(pc + 1, '\t')) ) {
		break;
	  }
	}
	if ( ic != muy || sscanf(pc, " %lf", &obj->muy) != 1 ) {
	  continue;
	}

    obj->line = strdup(line);
    ccarray_set_size(objects, ++size );
  }

  return 0;
}

/** uses fclose() if input file was uncompressed, and pclose() for compressed case */
static void close_file( FILE * input, compression_t compression )
{
  if ( input && input != stdin ) {
    if ( compression > compression_none ) {
      pclose(input);
    }
    else {
      fclose(input);
    }
  }
}


static int partition (double * m, int a, int b)
	{
	  int i = a;
	  int j;
	  for (j = a; j <= b; j++ )    // просматриваем с a по b
	   {
		 if (m[j] <= m[b])            // если элемент m[j] не превосходит m[b],
		  {
			//swap(m[i],m[j]);         // меняем местами m[j] и m[a], m[a+1], m[a+2] и так далее...
			double temp=m[i];						  // то есть переносим элементы меньшие m[b] в начало,
			m[i]=m[j];						  // а затем и сам m[b] «сверху»
			m[j]=temp;
			i++;                      // таким образом последний обмен: m[b] и m[i], после чего i++
		  }
	   }
	  return i-1;                     // в индексе i хранится <новая позиция элемента m[b]> + 1
	}

static void quicksort (double * m, int a, int b) // a - начало подмножества, b - конец
	{                                        // для первого вызова: a = 0, b = <элементов в массиве> - 1
	 if (a >= b) return;
	 int c = partition (m, a, b);
	 //printf("a=%d\tb=%d\tc=%d\n",a,b,c);

	 quicksort (m, a, c-1);
	 quicksort (m, c+1, b);

	}

static int get_index(double z, double max, double min, int n)
{
	if(z>max ){
		return n-1;
	}else if (z<min){
		return 0;
	}else {return (int)((z-min)/(max-min)*(n-1));}
}
///////////////////////////////////////////////////////////////////////////////////////////
static int med_array(int i, int j, int minx,int miny,int maxx,int maxy, int count_stars,int fit)
{
	int x = minx;
    int y = miny;
    int n = 0;
    int all_n=0;
    int t = miny;
    int max_array = (maxx-minx+1)*(maxy-miny+1);
    //fprintf(stderr,"count_stars = %d\n",count_stars);
    double arr_dRA [max_array];
    double arr_dDEC [max_array];

    //fprintf(stderr,"count array= %d\n",count_stars);
    if(fit<2)
    {
    while( x <= maxx)
      {y=t;
      while( y <= maxy)
        {
			  if(map[x][y].count!=0)
			  {
				 arr_dRA[n]  =   map[x][y].mux;
				 arr_dDEC[n] =   map[x][y].muy;
				 all_n       +=  map[x][y].count;
				 ++n;
			  }
		  ++y;
        }
      ++x;
      }
    }else
    {
    	while( x <= maxx)
		  {y=t;
		  while( y <= maxy)
			{
    	 	 if(map[x][y].med_count!=0)
			 {
				  arr_dRA[n]  =   map[x][y].med_mux;
				  arr_dDEC[n] =   map[x][y].med_muy;
				  all_n       +=  map[x][y].med_count;
				  ++n;
			 }
    	 	 ++y;
			}
		  ++x;
		  }
    }

      if(n >= count_stars)
        {
    	  quicksort(arr_dRA,0,n-1);
    	  quicksort(arr_dDEC,0,n-1);

    	// устанавливаем блокировку
		//pthread_mutex_lock(&lock);
			// изменяем глобальную переменную

    	if(n%2==0)
          {
			int nn = n/2;
			map[i][j].med_mux   = (arr_dRA[nn]+arr_dRA[nn-1])/2;
			map[i][j].med_muy   = (arr_dDEC[nn]+arr_dDEC[nn-1])/2;
			map[i][j].med_count = n;
          }
			else
			{
			  int nn = (int)((n-1)/2);
			  map[i][j].med_mux   = arr_dRA[nn];
			  map[i][j].med_muy   = arr_dDEC[nn];
			  map[i][j].med_count = n;
			}
  			// снимаем блокировку
   			//pthread_mutex_unlock(&lock);

        }

      return n;
}
static int med_array_w(int i, int j, int minx,int miny,int maxx,int maxy, int count_stars, int cw, int fit)
{
	int x = minx;
    int y = miny;
    int n = 0;

    int t = miny;
    //max_count objects for array
	int max_array = (int)(((maxx-minx)+1)*((maxy-miny)+1)*(cw+1)/2)+cw;
	//max distance objects with weight one
	int max_nw = (int)(sqrt((maxx-minx)*(maxx-minx)+(maxy-miny)*(maxy-miny))/2);

    //fprintf(stderr,"count_stars = %d\n",count_stars);
    double arr_dRA [max_array];
    double arr_dDEC [max_array];

    //fprintf(stderr,"count array= %d\n",count_stars);
    if(fit<2)
    {
    while( x <= maxx)
      {y=t;
      while( y <= maxy)
        {
           if(map[x][y].count!=0)
             {
        	    int xx = (int)abs(i-x);
				int yy = (int)abs(j-y);
				//distance carent objects to  i,j
				int nb = (int)((1-(sqrt(xx*xx+yy*yy)/max_nw))*cw)+1;
				//fprintf(stderr,"nb = %d\n",nb);
				while(nb > 0)
				{
					arr_dRA[n]  =   map[x][y].mux;
					arr_dDEC[n] =   map[x][y].muy;
					--nb;
					++n;
				}
             }
            ++y;
        }
      ++x;
      }
    }
    else
    {
        while( x <= maxx)
          {y=t;
          while( y <= maxy)
            {
               if(map[x][y].med_count!=0)
                 {
            	    int xx = (int)abs(i-x);
    				int yy = (int)abs(j-y);
    				//distance carent objects to  i,j
    				int nb = (int)((1-(sqrt(xx*xx+yy*yy)/max_nw))*cw)+1;
    				//fprintf(stderr,"nb = %d\n",nb);
    				while(nb > 0)
    				{
    					arr_dRA[n]  =   map[x][y].med_mux;
    					arr_dDEC[n] =   map[x][y].med_muy;
    					--nb;
    					++n;
    				}
                 }
                ++y;
            }
          ++x;
          }
    }

      if(n >= count_stars)
        {
    	  quicksort(arr_dRA,0,n-1);
    	  quicksort(arr_dDEC,0,n-1);

    	// устанавливаем блокировку
		//pthread_mutex_lock(&lock);
			// изменяем глобальную переменную

    	if(n%2==0)
          {
			int nn = n/2;
			map[i][j].med_mux   = (arr_dRA[nn]+arr_dRA[nn-1])/2;
			map[i][j].med_muy   = (arr_dDEC[nn]+arr_dDEC[nn-1])/2;
			map[i][j].med_count = n;
          }
			else
			{
			  int nn = (int)((n-1)/2);
			  map[i][j].med_mux   = arr_dRA[nn];
			  map[i][j].med_muy   = arr_dDEC[nn];
			  map[i][j].med_count = n;
			}
  			// снимаем блокировку
   			//pthread_mutex_unlock(&lock);
        }

return n;
}
void * create_map_p(void *arg)
{
	ARG_DATA * a = (ARG_DATA *) arg;
	int minx = a->minx;
	int miny = a->miny;
	int maxx = a->maxx;
	int maxy = a->maxy;
	int count_stars = a->mincountmed;
	int dx = a->dx;
	int dy = a->dy;

	int i;
	int j;

	if(a->fit==0)
	{
	for( i=minx; i<maxx;++i)
	    for( j=miny; j<maxy;++j)
	    {
	    	//fprintf(stderr,"count_stars = %d\n",count_stars);

	    	int mindx = i-dx;
			int mindy = j-dy;
			int maxdx = i+dx;
			int maxdy = j+dy;

			if(mindx<0){mindx=0;}
			if(mindy<0){mindy=0;}
			if(maxdx>=a->nx){maxdx=(a->nx-1);}
			if(maxdy>=a->ny){maxdy=(a->ny-1);}


	    	//fprintf(stderr,"count_stars = %d\n",count_map);
	    	//fprintf(stderr,"pthread_create = %d\n",i+j);
			int nstars =0;
			if(a->cw<2)
			{
				nstars =med_array(i,j,mindx,mindy,maxdx,maxdy,count_stars,a->fit);
			}else
			{
				nstars =med_array_w(i,j,mindx,mindy,maxdx,maxdy,count_stars,a->cw,a->fit);
			}
	    	if(nstars>=count_stars)
	    	{
	    		// устанавливаем блокировку
	    		pthread_mutex_lock(&lock);
	    		// изменяем глобальную переменную
	    		++count_map;
	    		// снимаем блокировку
	    		pthread_mutex_unlock(&lock);
	    		//fprintf(stderr,"good creat median filter for %d %d slot!\n",i,j);

	    	}
	    	//else{fprintf(stderr,"Can not creat median filter for %d stars!\n",nstars);}
	    }
	}else
		{
		for( i=minx; i<maxx;++i)
			    for( j=miny; j<maxy;++j)
			    {
			    	if(map[i][j].med_count==0)
			    	{
						//fprintf(stderr,"count_stars = %d\n",count_stars);

						int mindx = i-dx;
						int mindy = j-dy;
						int maxdx = i+dx;
						int maxdy = j+dy;

						if(mindx<0){mindx=0;}
						if(mindy<0){mindy=0;}
						if(maxdx>=a->nx){maxdx=(a->nx-1);}
						if(maxdy>=a->ny){maxdy=(a->ny-1);}


						//fprintf(stderr,"count_stars = %d\n",count_map);
						//fprintf(stderr,"pthread_create = %d\n",i+j);
						int nstars =0;
						if(a->cw<2)
							{
							nstars =med_array(i,j,mindx,mindy,maxdx,maxdy,count_stars,a->fit);
							}else
							{
							nstars =med_array_w(i,j,mindx,mindy,maxdx,maxdy,count_stars,a->cw,a->fit);
							}
						if(nstars>=count_stars)
						{
							// устанавливаем блокировку
							pthread_mutex_lock(&lock);
							// изменяем глобальную переменную
							if(++count_map%10000==0)
								{
								fprintf(stderr,"creat median filter slot! %d \r",count_map);

								}
							// снимаем блокировку
							pthread_mutex_unlock(&lock);
							//fprintf(stderr,"good creat median filter for %d %d slot!\n",i,j);

						}
									    	}
			    	//else{fprintf(stderr,"Can not creat median filter for %d stars!\n",nstars);}
			    }

		}

	//fprintf(stderr,"good creat median filter for %d slot!\n",count_map);
	//fprintf(stderr,"good creat median filter!\n");

 return NULL;
}

////////////////////////////////////////////////////////////////////////////////////
static int create_map(int p,int nx,int ny,int dx,int dy, int count_stars, int cw, int fit)
{
	//fprintf(stderr,"count_stars = %d\n",count_stars);
	//инициализация исключающей блокировки
	pthread_mutex_init(&lock, NULL);

	//массив идентификаторов потоков
	pthread_t thr[p];
	//fprintf(stderr,"count_stars = %d\n",count_stars);
	//int mn=0;

	int nxp = (int)(nx/p);
	ARG_DATA arg[p];
	int i;
	for( i=0; i<p;++i)
	    {
	    	//fprintf(stderr,"count_stars = %d\n",count_stars);

	    	arg[i].minx = i*nxp;
	    	arg[i].miny = 0;
	    	arg[i].maxx = (i+1)*nxp;
	    	arg[i].maxy = ny;
	    	arg[i].mincountmed = count_stars;
	    	arg[i].dx	= dx;
	    	arg[i].dy	= dy;
	    	arg[i].nx	= nx;
	    	arg[i].ny	= ny;
	    	arg[i].cw	= cw;
	    	arg[i].fit	= fit;



	    	int pth=0;
	    	//fprintf(stderr,"count_stars = %d\n",count_map);
	    	//fprintf(stderr,"pthread_create = %d\n",i+j);
	    	pth = pthread_create(&thr[i], NULL, create_map_p,(void *)&arg[i]);
	    	//fprintf(stderr,"pthread_create = %d\n",i);
	    	if(pth==1)
	    	{fprintf(stderr,"ERROR Create thread %d\n",i);}
	    }

	for(i = 0; i<p; ++i)
	   {
		   pthread_join(thr[i], NULL);
	   }

	//fprintf(stderr,"good creat %d slot for %d thread!!\n",arg[i].count_map,i);


 return count_map;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void show_usage( FILE * output, int argc, char * argv[] )
{
  fprintf(output, "2D Median filter UTILITY\n");
  fprintf(output, "USAGE:\n");
  fprintf(output, "  %s OPTIONS FILE1 FILE2\n", basename(argv[0]));
  fprintf(output, "OPTIONS:\n");
  fprintf(output, "  x1=	<double>    one-based colum number of x column in first file\n");
  fprintf(output, "  y1=	<double>    one-based colum number of y column in first file\n");
  fprintf(output, "  dmux1=	<double>    one-based colum number of dmux column in first file\n");
  fprintf(output, "  dmuy1=	<double>    one-based colum number of dmuy column in first file\n");
  fprintf(output, "  x2=	<double>    one-based colum number of x column in second file\n");
  fprintf(output, "  y2=	<double>    one-based colum number of y column in second file\n");
  fprintf(output, "  ra2=	<double>    one-based colum number of mux column in second file\n");
  fprintf(output, "  dec2=	<double>    one-based colum number of muy column in second file\n");
  fprintf(output, "  nx=	<integer>   count bins on x (max=1000)\n");
  fprintf(output, "  ny=	<integer>   count bins on y (max=1000)\n");
  fprintf(output, "  dx=	<integer>   count bins on shift x\n");
  fprintf(output, "  dy=	<integer>   count bins on shift y\n");
  fprintf(output, "  cs=	<integer>   minimum counts stars for median filter\n");
  fprintf(output, "  cw=	<integer>   counts of weighing for median filter\n");
  fprintf(output, "  cp=	<integer>   used counts of THREADS\n");
  fprintf(output, "  -int               Interpolation used\n");
  fprintf(output, "  -ext               extended used\n");
  fprintf(output, "  -i                 Invert match\n");
  fprintf(output, "  -d                 Print corrected difference\n");
  fprintf(output, "  -v                 Print some diagnostic messages to stderr (verbose mode)\n");
  fprintf(output, "  \n");

  UNUSED(argc);
}

int main(int argc, char *argv[])
{
	/* TODO: Fix this bug */

	int beverbose = 0;
	int usedfiles = 0;
	int interpolation = 0;
	int extended = 0;
	int invert_match = 0;
	int printcor_d =0;

	 compression_t compression[2] =
	    { compression_unknown, compression_unknown };


	const char * fname[2] =
	    { NULL, NULL };

	  FILE * fp[2] =
	    { NULL, NULL };

	ccarray_t * list[2] =
	{ NULL, NULL };

	size_t capacity[2] =
	{ 15000000, 15000000 };

	char head[2][MAX_HEADER_LENGTH];
	size_t size1, size2, pos1;
	const obj_t * obj1, * obj2;

	int x[2] ={ -1, -1 };
	int y[2] ={ -1, -1 };

	int dra[2] ={ -1, -1 };
	int ddec[2] ={ -1, -1 };

	int nx = 0; // count bins on x (max_nx=1000)
	int ny = 0; // count bins on y (max_ny=1000)

	int dx = 0; // count bins on shift x
	int dy = 0; // count bins on shift y

	int cs = 0; // min count stars for median filter
	int cw = 0; // count of weighting for median filter
	int cp = 0; // used count of THREAD

	double max_x = 0,min_x=1000000,max_y=0,min_y=1000000;

	/* parse command line */

	int i;
	for ( i = 1; i < argc; ++i )
	  {
		 if ( strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-help") == 0 ) {
		      show_usage(stdout, argc, argv);
		      return 0;
		    }

		/* read used columns from file1*/
		if ( strncmp(argv[i], "x1=", 3) == 0 )
	        {
	          if ( sscanf(argv[i] + 3, "%d", &x[0]) != 1 || x[0] < 1 )
	          {
	            fprintf(stderr, "Invalid value of %s\n", argv[i]);
	            return 1;
	          }
	        }
	    else if ( strncmp(argv[i], "y1=", 3) == 0 )
        {
          if ( sscanf(argv[i] + 3, "%d", &y[0]) != 1 || y[0] < 1 )
          {
            fprintf(stderr, "Invalid value of %s\n", argv[i]);
            return 1;
          }
        }
	    else if ( strncmp(argv[i], "dra1=", 5) == 0 )
        {
          if ( sscanf(argv[i] + 5, "%d", &dra[0]) != 1 || dra[0] < 1 )
          {
            fprintf(stderr, "Invalid value of %s\n", argv[i]);
            return 1;
          }
        }
	    else if ( strncmp(argv[i], "ddec1=", 6) == 0 )
        {
          if ( sscanf(argv[i] + 6, "%d", &ddec[0]) != 1 || ddec[0] < 1 )
          {
            fprintf(stderr, "Invalid value of %s\n", argv[i]);
            return 1;
          }
        }/* read used columns from file1*/
	    else if ( strncmp(argv[i], "x2=", 3) == 0 )
        {
          if ( sscanf(argv[i] + 3, "%d", &x[1]) != 1 || x[1] < 1 )
          {
            fprintf(stderr, "Invalid value of %s\n", argv[i]);
            return 1;
          }
        }
	    else if ( strncmp(argv[i], "y2=", 3) == 0 )
        {
          if ( sscanf(argv[i] + 3, "%d", &y[1]) != 1 || y[1] < 1 )
          {
            fprintf(stderr, "Invalid value of %s\n", argv[i]);
            return 1;
          }
        }
	    else if ( strncmp(argv[i], "ra2=", 4) == 0 )
        {
          if ( sscanf(argv[i] + 4, "%d", &dra[1]) != 1 || dra[1] < 1 )
          {
            fprintf(stderr, "Invalid value of %s\n", argv[i]);
            return 1;
          }
        }
	    else if ( strncmp(argv[i], "dec2=", 5) == 0 )
        {
          if ( sscanf(argv[i] + 5, "%d", &ddec[1]) != 1 || ddec[1] < 1 )
          {
            fprintf(stderr, "Invalid value of %s\n", argv[i]);
            return 1;
          }
        } /* read parameters median filters */
	    else if ( strncmp(argv[i], "--nx=", 5) == 0 )
        {
          if ( sscanf(argv[i] + 5, "%d", &nx) != 1 || nx < 1 )
          {
            fprintf(stderr, "Invalid value of %s\n", argv[i]);
            return 1;
          }
        }
	    else if ( strncmp(argv[i], "--ny=", 5) == 0 )
			{
			  if ( sscanf(argv[i] + 5, "%d", &ny) != 1 || ny < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
	    else if ( strncmp(argv[i], "--dx=", 5) == 0 )
			{
			  if ( sscanf(argv[i] + 5, "%d", &dx) != 1 || dx < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
	    else if ( strncmp(argv[i], "--dy=", 5) == 0 )
			{
			  if ( sscanf(argv[i] + 5, "%d", &dy) != 1 || dy < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
	    else if ( strncmp(argv[i], "--cs=", 5) == 0 )
			{
			  if ( sscanf(argv[i] + 5, "%d", &cs) != 1 || cs < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
	    else if ( strncmp(argv[i], "--cw=", 5) == 0 )
			{
			  if ( sscanf(argv[i] + 5, "%d", &cw) != 1 || cw < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
	    else if ( strncmp(argv[i], "--cp=", 5) == 0 )
			{
			  if ( sscanf(argv[i] + 5, "%d", &cp) != 1 || cp < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
	    else if ( strcmp(argv[i], "-v") == 0 ) {
	          beverbose = 1;
	        }
	    else if ( strcmp(argv[i], "-i") == 0 ) {
	          invert_match = 1;
	        }
	    else if ( strcmp(argv[i], "-d") == 0 ) {
	    	printcor_d = 1;
	        }
	    else if ( strcmp(argv[i], "-int") == 0 ) {
	    	  interpolation = 1;
	    	}
	    else if ( strcmp(argv[i], "-ext") == 0 ) {
			  extended = 1;
			}
	    else if ( !fname[0] ) {
	          fname[0] = argv[i];
	          usedfiles = 1;
	        }
	    else if ( !fname[1] ) {
	          fname[1] = argv[i];
	          usedfiles = 1;
	        }
	    else
	        {
	          fprintf(stderr, "Invalid argument %s. Try %s --help\n", argv[i], argv[0]);
	          return 1;
	        }
	  }


	/* check command line inputs */
	  if ( !fname[0] || !fname[1] ) {
	    fprintf(stderr,"Two input file names expected\n");
	    show_usage(stderr, argc, argv);
	    return -1;
	  }

	  if ( x[0] < 1 ) {
	      fprintf(stderr,"x1 argument is mandatory\n");
	      show_usage(stderr, argc, argv);
	      return -1;
	    }
	  if ( y[0] < 1 ) {
	  	      fprintf(stderr,"y1 argument is mandatory\n");
	  	      show_usage(stderr, argc, argv);
	  	      return -1;
	  	    }
	  if ( dra[0] < 1 ) {
		  fprintf(stderr,"dmux1 argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( ddec[0] < 1 ) {
		  fprintf(stderr,"dmuy1 argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}

	  if ( x[1] < 1 ) {
		  fprintf(stderr,"x2 argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( y[1] < 1 ) {
		  fprintf(stderr,"y2 argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( dra[1] < 1 ) {
		  fprintf(stderr,"dmux2 argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( ddec[1] < 1 ) {
		  fprintf(stderr,"dmuy2 argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}

	  if ( nx < 1 || nx > 1000) {
	    fprintf(stderr,"nx argument is mandatory\n");
	    show_usage(stderr, argc, argv);
	    return -1;
	  }

	  if ( ny < 1 || ny > 1000) {
	    fprintf(stderr,"ny argument is mandatory\n");
	    show_usage(stderr, argc, argv);
	    return -1;
	  }

	  if ( dx < 1 ) {
	    fprintf(stderr,"dx argument is mandatory\n");
	    show_usage(stderr, argc, argv);
	    return -1;
	  }

	  if ( dy < 1 ) {
	    fprintf(stderr,"dy argument is mandatory\n");
	    show_usage(stderr, argc, argv);
	    return -1;
	  }

	  if ( cs < 1 ) {
	    fprintf(stderr,"cs argument is mandatory\n");
	    show_usage(stderr, argc, argv);
	    return -1;
	  }

	  if ( cw < 1 ) {
	  	    fprintf(stderr,"cw argument is mandatory\n");
	  	    show_usage(stderr, argc, argv);
	  	    return -1;
	  	  }

	  if ( cp < 1 ) {
	  	    fprintf(stderr,"cp argument is mandatory\n");
	  	    show_usage(stderr, argc, argv);
	  	    return -1;
	  	  }

	if (usedfiles)
	{
		/* check if input files are readable */
		  for ( i = 0; i < 2; ++i )
		  {
		    if ( access(fname[i], R_OK) != 0 ) {
		      fprintf(stderr, "Can't read %s: %s\n", fname[i], strerror(errno));
		      return -1;
		    }
		  }

		 /* allocate memory storage */
		  for ( i = 0; i < 2; ++i )
		  {
		    if ( !(list[i] = ccarray_create(capacity[i], sizeof(obj_t))) ) {
		      fprintf(stderr, "ccarray_create(capacity=%zu) fails: %s\n", capacity[i], strerror(errno));
		      return -1;
		    }
		  }

		/* load input files */
		  for ( i = 0; i < 2; ++i )
		    {
		      if ( beverbose ) {
		        fprintf(stderr,"loading %s....\n", fname[i]);
		      }

		      if ( !(fp[i] = open_file(fname[i], &compression[i])) ) {
		        fprintf(stderr, "Can't read '%s': %s\n", fname[i], strerror(errno));
		        return -1;
		      }
		      if ( load_objects(fp[i], x[i], y[i], dra[i], ddec[i], head[i], list[i]) != 0 ) {
		        fprintf(stderr, "Can't load %s\n", fname[i]);
		        return -1;
		      }
		      close_file( fp[i], compression[i] );

		      if ( beverbose ) {
		        fprintf(stderr,"%s: %zu objects\n", fname[i], ccarray_size(list[i]));
		      }

		      if ( ccarray_size(list[i]) == capacity[i] && !feof(fp[i]) ) {
		        fprintf(stderr, "load_objects() fails for %s: Too many objects. Try increase capacity1 (currently %zu)\n",
		            fname[i], capacity[i]);
		        return -1;
		      }
		    }

		  size1 = ccarray_size(list[0]);
		  size2 = ccarray_size(list[1]);

		  int nstars = size1;

		  if (nstars<100){
			  fprintf(stderr,"Imposibiliti to creating median map for %d stars\n",nstars);
			  return 0;
		  }
		  /* search pairs */
		    if ( beverbose ) {
		    	fprintf(stderr,"Used file1 '%s' file2 '%s' and argumernts: x1=%d\ty1=%d\tdra1=%d\tddec1=%d\tx2=%d\ty2=%d\tdra2=%d\tddec2=%d\tnx=%d\tny=%d\tdx=%d\tdy=%d\tcs=%d\tcw=%d\tcp=%d\n",fname[0],fname[1],x[0],y[0],dra[0],ddec[0],x[1],y[1],dra[1],ddec[1],nx,ny,dx,dy,cs,cw,cp);
		    	fprintf(stderr,"creating median map for %d stars\n",nstars);
		    }

		    /* found max and min*/
		    for ( pos1 = 0; pos1 < size1; ++pos1 )
		      {
		    	obj1 = ccarray_peek(list[0], pos1);

		    	if(obj1->x>max_x)
					{max_x = obj1->x;}
				if(obj1->x<min_x)
					{min_x = obj1->x;}
				if(obj1->y>max_y)
					{max_y = obj1->y;}
				if(obj1->y<min_y)
					{min_y = obj1->y;}
		      }

		    if ( beverbose ) {
		    fprintf(stderr,"Min_x=%lf\tMax_x=%lf\tMin_y=%lf\tMax_y=%lf\n",min_x,max_x,min_y,max_y);
		    }
		    int xx,yy;
		    	for(xx=0;xx<nx;xx++)
		    		for(yy=0;yy<ny;yy++)
		    			{
		    			map[xx][yy].x 		 = 0;
		    			map[xx][yy].y 	     = 0;
		    			map[xx][yy].mux 	 = 0;
		    			map[xx][yy].muy 	 = 0;
		    			map[xx][yy].med_mux  = 0;
		    			map[xx][yy].med_muy  = 0;
		    			map[xx][yy].med_count= 0;
		    			map[xx][yy].count 	 = 0;

		    			}
		    	//write stars data

		    	int ncor =0;
		    	int ddx=dx,ddy=dy;

		    	if(nstars<5000)
		    	{
		    		ddx=dx;
		    		ddy=dy;
		    		nx=nx/5;
		    		ny=ny/5;

		    	}else
		    	{
		    		if(nstars<30000)
		    		{
		    			ddx=2*dx;
		    			ddy=2*dy;
		    		}
		    	}
		    	for ( pos1 = 0; pos1 < size1; ++pos1 )
					{
						obj1 = ccarray_peek(list[0], pos1);

						xx = get_index(obj1->x,max_x, min_x, nx);
						yy = get_index(obj1->y,max_y, min_y, ny);

						map[xx][yy].count ++;
						map[xx][yy].x  		+= (obj1->x-map[xx][yy].x)/map[xx][yy].count;
						map[xx][yy].y  		+= (obj1->y-map[xx][yy].y)/map[xx][yy].count;
						map[xx][yy].mux 	+= (obj1->mux-map[xx][yy].mux)/map[xx][yy].count;
						map[xx][yy].muy 	+= (obj1->muy-map[xx][yy].muy)/map[xx][yy].count;

					}

				if ( beverbose )
				{
					fprintf(stderr,"Used dx = %d dy = %d for correction without weight %d\n",ddx, ddy, cw);
					fprintf(stderr,"Used %d THREADS for maping correction\n",cp);
				}

				ncor = create_map(cp,nx,ny,ddx,ddy,cs,cw,0);
				if ( beverbose ) {
						fprintf(stderr,"Create %d correction dx =%d dy =%d\n",ncor,ddx,ddy);
					}
				if (interpolation==1)
					{
						int d =1;
						//extended fitting windows
						while (d<dx && d<dy && ncor<nx*ny)
						{
							ncor = create_map(cp,nx,ny,ddx+d,ddy+d,cs,cw,1);
							if ( beverbose ) {
								fprintf(stderr,"Create extended %d correction dx =%d dy =%d\n",ncor,ddx+d,ddy+d);
							}
							++d;
						}
						//interpolation
						if(ncor>0.15*nx*ny)
						{
							d=1;
							while ( ncor<nx*ny && d<ddx)
								{
									ncor = create_map(cp,nx,ny,dx/2+d,dy/2+d,cs,cw,2);
									if ( beverbose ) {
										fprintf(stderr,"Create interpolation %d correction for dx =%d dy =%d\n",ncor,dx/2+d,dy/2+d);
									}
									++d;
								}
						}else{	if ( beverbose ) {
							fprintf(stderr,"ERROR Can not Create interpolation for %d\n",ncor);
						}
					}
					}else{
						if(extended==1)
						{
							int d=1;
							while ( ncor<nx*ny)
							{
								ncor = create_map(cp,nx,ny,ddx+d,ddy+d,cs,cw,1);
								if ( beverbose )
								{
									fprintf(stderr,"Create extended %d correction dx =%d dy =%d\n",ncor,ddx+d,ddy+d);
								}
								++d;
							}
						}
					}


		    	/*corection data*/
				if(!invert_match)
				{
					if(printcor_d)
					{
						/* printf new head */
						printf("cra\tcdec\tcor_ra\tcor_dec\tmed_count\tcor_dra\tcor_ddec\t%s\n", head[1]);

						for ( pos1 = 0; pos1 < size2; ++pos1 )
						{
							obj2 = ccarray_peek(list[1], pos1);

							xx = get_index(obj2->x,max_x, min_x, nx);
							yy = get_index(obj2->y,max_y, min_y, ny);

						  printf("%+14.10f\t%+14.10f\t%+9.6f\t%+9.6f\t%+9.6f\t%+9.6f\t%d\t%s\n", (obj2->mux+map[xx][yy].med_mux/3600/180*PI/cos(obj2->muy)),(obj2->muy+map[xx][yy].med_muy/3600/180*PI), map[xx][yy].med_mux,map[xx][yy].med_muy,map[xx][yy].mux-map[xx][yy].med_mux,map[xx][yy].muy-map[xx][yy].med_muy,map[xx][yy].med_count,obj2->line);
						}
					}else
					{
						/* printf new head */
						printf("cra\tcdec\tcor_ra\tcor_dec\tmed_count\t%s\n", head[1]);

						for ( pos1 = 0; pos1 < size2; ++pos1 )
						{
							obj2 = ccarray_peek(list[1], pos1);

							xx = get_index(obj2->x,max_x, min_x, nx);
							yy = get_index(obj2->y,max_y, min_y, ny);

						  printf("%+14.10f\t%+14.10f\t%+9.6f\t%+9.6f\t%d\t%s\n", (obj2->mux+map[xx][yy].med_mux/3600/180*PI/cos(obj2->muy)),(obj2->muy+map[xx][yy].med_muy/3600/180*PI), map[xx][yy].med_mux,map[xx][yy].med_muy,map[xx][yy].med_count,obj2->line);
						}

					}
				}else{
					/* printf new head */
				    printf("i\tj\tmed_mux\tmed_muy\tmed_count\n");

					int k,m;
					for(k=0;k<nx;++k)
						for(m=0;m<ny;++m)
						{
							if(map[k][m].med_count!=0)
							{
								printf("%d\t%d\t%+9.3f\t%+9.3f\t%d\n", k,m, map[k][m].med_mux,map[k][m].med_muy,map[k][m].med_count);
							}

						}
				}

	}
	return 0;
}
