#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lodepng.cpp"
//#include <iostream>
void freeImage2D(unsigned char **mat);
unsigned char *convert2Dto1D(unsigned char **imgR, unsigned char **imgG,unsigned char **imgB, unsigned int h, unsigned int w);
unsigned char *convert2Dto1D(unsigned char **img2D, unsigned int h, unsigned int w);
unsigned char **getChannel2D(unsigned char *img1D, unsigned int h, unsigned int w,unsigned int index);
unsigned char *readPNG(const char* filename, unsigned int &width, unsigned int &height,unsigned int &bitdepth, unsigned int &bitsXpixel, unsigned int &channels,unsigned int &isGrey, unsigned int &haveAlpha);
int savePNG(const char *filename, unsigned char *img, int width, int height);

unsigned char **createMatrix(int h,int w) {
	unsigned char **mat;
	int i;
	mat    = (unsigned char **) malloc( h*sizeof(unsigned char *));
    if(mat==NULL) return(NULL);
    mat[0] = (unsigned char *) malloc( h*w*sizeof(unsigned char));
    if(mat[0]==NULL) return(NULL);
    for(i=1; i<h; ++i)
        mat[i] = mat[i-1] + w;
    
	return(mat);
}


double bicubicpol(double x,double y,double p[4][4]){

	double a00, a01, a02, a03;
 	double a10, a11, a12, a13;
 	double a20, a21, a22, a23;
 	double a30, a31, a32, a33;
 	double x2 = x * x;
	double x3 = x2 * x;
	double y2 = y * y;
	double y3 = y2 * y;

	a00 = p[1][1];
	a01 = -.5*p[1][0] + .5*p[1][2];
	a02 = p[1][0] - 2.5*p[1][1] + 2*p[1][2] - .5*p[1][3];
	a03 = -.5*p[1][0] + 1.5*p[1][1] - 1.5*p[1][2] + .5*p[1][3];
	a10 = -.5*p[0][1] + .5*p[2][1];
	a11 = .25*p[0][0] - .25*p[0][2] - .25*p[2][0] + .25*p[2][2];
	a12 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + .5*p[2][0] - 1.25*p[2][1] + p[2][2] - .25*p[2][3];
	a13 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .25*p[2][0] + .75*p[2][1] - .75*p[2][2] + .25*p[2][3];
	a20 = p[0][1] - 2.5*p[1][1] + 2*p[2][1] - .5*p[3][1];
	a21 = -.5*p[0][0] + .5*p[0][2] + 1.25*p[1][0] - 1.25*p[1][2] - p[2][0] + p[2][2] + .25*p[3][0] - .25*p[3][2];
	a22 = p[0][0] - 2.5*p[0][1] + 2*p[0][2] - .5*p[0][3] - 2.5*p[1][0] + 6.25*p[1][1] - 5*p[1][2] + 1.25*p[1][3] + 2*p[2][0] - 5*p[2][1] + 4*p[2][2] - p[2][3] - .5*p[3][0] + 1.25*p[3][1] - p[3][2] + .25*p[3][3];
	a23 = -.5*p[0][0] + 1.5*p[0][1] - 1.5*p[0][2] + .5*p[0][3] + 1.25*p[1][0] - 3.75*p[1][1] + 3.75*p[1][2] - 1.25*p[1][3] - p[2][0] + 3*p[2][1] - 3*p[2][2] + p[2][3] + .25*p[3][0] - .75*p[3][1] + .75*p[3][2] - .25*p[3][3];
	a30 = -.5*p[0][1] + 1.5*p[1][1] - 1.5*p[2][1] + .5*p[3][1];
	a31 = .25*p[0][0] - .25*p[0][2] - .75*p[1][0] + .75*p[1][2] + .75*p[2][0] - .75*p[2][2] - .25*p[3][0] + .25*p[3][2];
	a32 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + 1.5*p[1][0] - 3.75*p[1][1] + 3*p[1][2] - .75*p[1][3] - 1.5*p[2][0] + 3.75*p[2][1] - 3*p[2][2] + .75*p[2][3] + .5*p[3][0] - 1.25*p[3][1] + p[3][2] - .25*p[3][3];
	a33 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .75*p[1][0] + 2.25*p[1][1] - 2.25*p[1][2] + .75*p[1][3] + .75*p[2][0] - 2.25*p[2][1] + 2.25*p[2][2] - .75*p[2][3] - .25*p[3][0] + .75*p[3][1] - .75*p[3][2] + .25*p[3][3];

		
	return (a00 + a01 * y + a02 * y2 + a03 * y3) +
		    (a10 + a11 * y + a12 * y2 + a13 * y3) * x +
		    (a20 + a21 * y + a22 * y2 + a23 * y3) * x2 +
		    (a30 + a31 * y + a32 * y2 + a33 * y3) * x3;
		       
}


unsigned char **imageInterpolate(int f,int h, int w, unsigned char **imgC) {
	double arr[4][4];
    int i,j;
    unsigned char **img2g = createMatrix(f*h,f*w);
  
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			arr[i][j]=0;
	
     for(int i=0; i<f*(h-1)-f+1; i++) {

        for(int j=0; j<f*(w-2)-f+1; j++) {
            
            for( int l = 0; l < 4; l++ ) {
            	
            	for( int k = 0; k < 4; k++ ) {
            			
            		if((int)i/f + l < h && (int)j/f + k < w) {
            			
            			arr[l][k]=(double)imgC[(int)i/f][(int) (j/f)];
            			
					}

        		}
        	}
        	
        	img2g[i][j] = (unsigned char)bicubicpol((double)(i%f)/f,(double)(j%f)/f ,arr);
        }
        
	}
	return img2g;
}

int main(int argc, char *argv[]) {
	
    if(argc < 2)
    {
        printf("De el nombre de la imagen.\n");
        return 0;
    }
    const char* filename = argv[1];
    char nameout[50];
    int f=atoi(argv[2]);
   

    unsigned char *image1D;
    unsigned int    w, h, bitdepth, bitsXpixel, channels, isGrey, haveAlpha;

    // Lectura de la imagen
    image1D = readPNG(filename, w, h, bitdepth, bitsXpixel, channels, isGrey, haveAlpha);

    printf("Dimensions:      %d x %d\n", h, w);
    printf("Depth Bits:      %d\n", bitdepth);
    printf("Bits per pixel:  %d\n", bitsXpixel);
    printf("Channels:        %d\n", channels);
    printf("Is Gray Scale?:  %d\n", isGrey);
    printf("Alpha Channel?:  %d\n", haveAlpha);
	sprintf(nameout, "%dx%s_%dx%d.png ",f, "_bc",h,w);
    // Arreglos 2D de cada componente de la imagen. Si la imagen esta
    // en escala de grises, los tres arreglos son iguales.
    // Este formato es mas util para trabajar con la imagen como si
    // fuera una matriz pues los valores de los pixeles se pueden obtener
    // como img[i][j] para i=0,...,h-1 y j=0,...,w-1
    unsigned char **imgR = getChannel2D(image1D, h, w, 0);
    //printf("1");
    unsigned char **imgG = getChannel2D(image1D, h, w, 1);
    //printf("2");
    unsigned char **imgB = getChannel2D(image1D, h, w, 2);
	//printf("3");    
    unsigned char **img2gR;
    unsigned char **img2gG;
    unsigned char **img2gB;
    //printf("4");

	img2gR = imageInterpolate(f,h,w,imgR);
	img2gG = imageInterpolate(f,h,w,imgG);
    img2gB = imageInterpolate(f,h,w,imgB);
   
  	
    // Se convierte cada canal a un arreglo 1D por separado para poder grabarlos como PNG
    //f*(h)-f+1, f*(w)-f+1
    unsigned char *imgR1 = convert2Dto1D(img2gR, f*h, f*w);
    unsigned char *imgG1 = convert2Dto1D(img2gG, f*h, f*w);
    unsigned char *imgB1 = convert2Dto1D(img2gB, f*h, f*w);

    // Se combinan los tres canales en un solo arreglo 1D para formar una imagen a color
   unsigned char *imgRGB1 = convert2Dto1D(img2gR, img2gG, img2gB, f*h, f*w);


    // Se graban todos los arreglos 1D como imagen PNG.
    //savePNG("canal_R.png", imgR1, w, h);
    //savePNG(nameout, imgG1, f*(w)-f+1, f*(h)-f+1);
    //savePNG("canal_B.png", imgB1, w, h);
    //f*(w)-f+1, f*(h)-f+1
	//savePNG("p.png", imgR1, f*w, f*h);
    savePNG(nameout, imgRGB1, f*w, f*h);

    free(image1D);
    freeImage2D(imgR);
    freeImage2D(imgG);
    freeImage2D(imgB);
    freeImage2D(img2gR);
    freeImage2D(img2gG);
    freeImage2D(img2gB);
    free(imgR1);
    free(imgG1);
    free(imgB1);
    free(imgRGB1);
    //free(img2gR);
    //free(img2gG);
   // free(img2gB);
    

    return(0);
}


// Escritura de un arreglo 1D en una imagen PNG
int savePNG(const char *filename, unsigned char *img, int width, int height)
{
    unsigned char *png;
    size_t pngsize;
    int error = lodepng_encode24(&png, &pngsize, img, width, height);
    if(!error)
    {
        lodepng_save_file(png, pngsize, filename);
    }

    if(error)
        printf("\tError %u al grabar la imagen %s: %s\n", error, filename,
                lodepng_error_text(error));

    printf("Finished: %s\n", filename);
    free(png);
    return(error);
}


unsigned char *readPNG(const char* filename, unsigned int &width, unsigned int &height,
             unsigned int &bitdepth, unsigned int &bitsXpixel, unsigned int &channels,
             unsigned int &isGrey, unsigned int &haveAlpha) {
  unsigned error;
  unsigned char* image;
  unsigned char* png = 0;
  size_t pngsize;
  LodePNGState state;

  lodepng_state_init(&state);

  error = lodepng_load_file(&png, &pngsize, filename);
  if(!error) error = lodepng_decode(&image, &width, &height, &state, png, pngsize);
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
  free(png);

  LodePNGColorMode& color = state.info_png.color;

  bitdepth   = color.bitdepth;
  bitsXpixel = lodepng_get_bpp(&color);
  channels   = lodepng_get_channels(&color);
  isGrey     = lodepng_is_greyscale_type(&color);
  haveAlpha  = lodepng_can_have_alpha(&color);

  lodepng_state_cleanup(&state);
  return(image);
}

// A partir de una imagen img1D codificada como mediante un arreglo 1D, esta funcion
// devuelve un arreglo 2D que tiene la componente de color indicada por index:
// index = 0  para recuperar el canal rojo
// index = 1  para recuperar el canal verde
// index = 2  para recuperar el canal azul
unsigned char **getChannel2D(unsigned char *img1D, unsigned int h, unsigned int w,
                             unsigned int index)
{
    unsigned int    i, j, k, l;
    unsigned int    channels = 3;
    unsigned char **mat;

    // Reservamos memoria
    mat    = (unsigned char **) malloc( h*sizeof(unsigned char *));
    if(mat==NULL) return(NULL);
    mat[0] = (unsigned char *) malloc( h*w*sizeof(unsigned char));
    if(mat[0]==NULL) return(NULL);
    for(i=1; i<h; ++i)
        mat[i] = mat[i-1] + w;

    // Lectura de los datos
    l = (channels + 1);
    for(i=0; i<h; i++) {
        k = i*w*(channels+1);
        for(j=0; j<w; j++) {
            mat[i][j] = img1D[j*l + index + k];
        }
    }

    return(mat);
}

// Convierte un arraglo 2D que tiene la informacion de una imagen a un arreglo 1D
// como se requiere para poderlo grabarlo como imagen PNG en escala de grises
unsigned char *convert2Dto1D(unsigned char **img2D, unsigned int h, unsigned int w)
{
    unsigned char *img1D = (unsigned char *) malloc(sizeof(unsigned char)*h*w*3);
    unsigned int   i, j, k, l;
    unsigned char  val;

    for(i=0; i<h; i++) {
        k = 3*w*i;
        for(j=0; j<w; j++) {
            val = (unsigned char) img2D[i][j];
            l   = 3*j + k;
            img1D[l]   = val;
            img1D[l+1] = val;
            img1D[l+2] = val;
        }
    }
    return(img1D);
}

// Combina la informacion de tres arreglos 2D para formar un arreglo 1D para una imagen a color
unsigned char *convert2Dto1D(unsigned char **imgR, unsigned char **imgG,
                             unsigned char **imgB, unsigned int h, unsigned int w)
{
    unsigned char *img1D = (unsigned char *) malloc(sizeof(unsigned char)*h*w*3);
    unsigned int   i, j, k, l;

    for(i=0; i<h; i++) {
        k = 3*w*i;
        for(j=0; j<w; j++) {
            l   = 3*j + k;
            img1D[l]   = imgR[i][j];
            img1D[l+1] = imgG[i][j];
            img1D[l+2] = imgB[i][j];
        }
    }
    return(img1D);
}


// Libera la memoria del arreglo bidimensional
void freeImage2D(unsigned char **mat) {
    free(mat[0]);
    free(mat);
}



