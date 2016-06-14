/*                                                                
**  Copyright (C) 2004-2008,2016  Smithsonian Astrophysical Observatory 
*/

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */


#include <dslib.h>
#include <dsnan.h>
#include <histlib.h>
#include <math.h>
#include <float.h>
#include <cxcregion.h>


typedef enum { LAPLACE, ROBERTS, PREWITT, SOBEL, ROBISON, KIRSCH, GAUSSIAN } Shapes;

typedef enum { X_GRAD, Y_GRAD, MAGNITUDE, ANGLE } Direction;


/* Using the dmtools/dmimgio routines removes lots of duplicate code that was
 * originally here.  Also allow us to keep track of NULL/NaN value pixels
 * more easily
 */
#include "dmimgio.h"


float *evaluate_kernel(Shapes shape, float width, Direction xy, long *kx, long *ky);
Shapes get_shape( char *kernel ); 
Direction get_dir( char *output_val );
double *slide_convovle(void *data, dmDataType dt, long *lAxes, short *mask, Shapes shape, float width, Direction dir );
double *combine_gradients( Direction dir, double *outdata_x, double *outdata_y, long *lAxes);
int dmimggrad(void);


Shapes get_shape( char *kernel )
{
    Shapes shape;

    if (strcmp(kernel, "laplace") == 0) {
        shape = LAPLACE;
    } else if (strcmp(kernel, "roberts") == 0) {
        shape = ROBERTS;
    } else if (strcmp(kernel, "prewitt") == 0) {
        shape = PREWITT;
    } else if (strcmp(kernel, "sobel") == 0) {
        shape = SOBEL;
    } else if (strcmp(kernel, "robison") == 0) {
        shape = ROBISON;
    } else if (strcmp(kernel, "kirsch") == 0) {
        shape = KIRSCH;
    } else if (strcmp(kernel, "gaussian") ==0) {
        shape = GAUSSIAN;
    } else {
        shape = -1;
    }
    
    return (shape);
}




float *evaluate_kernel(Shapes shape, float width, Direction xy, long *kx, long *ky)
{
    float *retval;
    *kx = 3;
    *ky = 3;
    retval = (float *) calloc(9, sizeof(float));


    switch (shape) {
    case LAPLACE:
        {
            float op[9] = { 0, 1, 0, 1, -4, 1, 0, 1, 0 };
            memcpy(retval, op, 9 * sizeof(float));
            break;
        }
    case ROBERTS:
        *kx = 2;
        *ky = 2;
        if (xy == X_GRAD) {
            float op[4] = { 1, 0, 0, -1 };
            memcpy(retval, op, 4 * sizeof(float));
        } else {
            float op[4] = { 0, 1, -1, 0 };
            memcpy(retval, op, 4 * sizeof(float));
        }
        break;
    case PREWITT:
        if (xy == X_GRAD) {
            float op[9] = { 1, 1, 1, 0, 0, 0, -1, -1, -1 };
            memcpy(retval, op, 9 * sizeof(float));
        } else {
            float op[9] = { -1, 0, 1, -1, 0, 1, -1, 0, 1 };
            memcpy(retval, op, 9 * sizeof(float));
        }
        break;
    case SOBEL:
        if (xy == X_GRAD) {
            float op[9] = { 1, 2, 1, 0, 0, 0, -1, -2, -1 };
            memcpy(retval, op, 9 * sizeof(float));
        } else {
            float op[9] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
            memcpy(retval, op, 9 * sizeof(float));
        }
        break;
    case ROBISON:
        if (xy == X_GRAD) {
            float op[9] = { 1, 1, 1, 1, -2, 1, -1, -1, -1 };
            memcpy(retval, op, 9 * sizeof(float));
        } else {
            float op[9] = { -1, 1, 1, -1, -2, 1, -1, 1, 1 };
            memcpy(retval, op, 9 * sizeof(float));
        }
        break;
    case KIRSCH:
        if (xy == X_GRAD) {
            float op[9] = { 3, 3, 3, 3, 0, 3, -5, -5, -5 };
            memcpy(retval, op, 9 * sizeof(float));
        } else {
            float op[9] = { -5, 3, 3, -5, 0, 3, -5, 3, 3 };
            memcpy(retval, op, 9 * sizeof(float));
        }
        break;
    case GAUSSIAN:
        {
            float nsigma = 4.0;
            long half_width = (nsigma*width)+0.5;
            *kx = (2*half_width)+1;
            *ky = *kx;
            free(retval);
            retval = (float *) calloc( ((*kx)*(*ky)), sizeof(float));
            long ix, iy, dx, dy;
            for (ix=0;ix< *kx;ix++) {
                for (iy=0;iy< *ky;iy++) {
                    dx = ix-half_width;
                    dy = iy-half_width;
                    long dd;
                    if ( xy == X_GRAD ) {
                        dd = dx;
                    } else {
                        dd = dy;
                    }
                    
                    retval[ix+iy*(*kx)] = (dd)*exp( -1.0*(dx*dx+dy*dy)/(width*width));
                    //printf("%g\n", retval[ix+iy*(*kx)] );
                } // end for iy
            } // end for ix


            break;
        }



    default:
        retval = NULL;
    }


    return (retval);

}

double *slide_convovle(void *data, dmDataType dt, long *lAxes, short *mask, Shapes shape, float width, Direction dir )
{

    long dkx, dky;    
    long xx, yy;
    long ix, iy;
    double *outdata_x;
    float *kernel;
    long kx, ky;
    outdata_x = (double *) calloc(lAxes[0] * lAxes[1], sizeof(double));
    if (NULL == (kernel = evaluate_kernel(shape, width, dir, &kx, &ky))) {
        err_msg("ERROR: Unsupported gradient function ");
        return (NULL);
    }

    dkx = kx / 2 - 1;
    dky = ky / 2 - 1;
    for (xx = lAxes[0]; xx--;) {
        for (yy = lAxes[1]; yy--;) {
            long pix;
            pix = 0;

            for (ix = kx; ix--;) {
                for (iy = ky; iy--;) {
                    if (kernel[pix]) {
                        double dater;
                        long ax, ay;
                        ax = xx + ix - dkx;
                        ay = yy + iy - dky;

                        if ((ax < 0) || (ay < 0) ||
                            (ax >= lAxes[0]) || (ay >= lAxes[1])) {
                            continue;
                        }
                        dater = get_image_value(data, dt, ax, ay, lAxes, mask);
                        if (ds_dNAN(dater)) {
                            continue;
                        }

                        outdata_x[xx + yy * lAxes[0]] += (kernel[pix] * dater);
                    }       /* end kernel[pix] */
                    pix++;
                }           /* end for iy */
            }               /* end for ix */
        }                   /* end for yy */
    }                       /* end for xx */

    return (outdata_x);

}

Direction get_dir( char *output_val )
{
    Direction dir;
    
    switch (*output_val) {  // switch on first letter in direction string
        case 'x':
            dir = X_GRAD;
            break;
        case 'y':
            dir = Y_GRAD;
            break;
        case 'm':
            dir = MAGNITUDE;
            break;
        case 'a':
            dir = ANGLE;
            break;
        default:
            err_msg("ERROR: Unsupported output value '%s'\n", output_val);
            return (-1);
    }
    return (dir);
}


double *combine_gradients( Direction dir, double *outdata_x, double *outdata_y, long *lAxes)
{
    long ii;

    switch (dir) {
        case X_GRAD:
            break;
        case Y_GRAD:
            outdata_x = outdata_y;
            break;
        case MAGNITUDE:{
                for (ii = (lAxes[0] * lAxes[1]); ii--;) {
                    outdata_x[ii] =
                        sqrt(pow(outdata_x[ii], 2.0) +
                             pow(outdata_y[ii], 2.0));
                }
                break;
            }
        case ANGLE:{
                double dat;
                for (ii = (lAxes[0] * lAxes[1]); ii--;) {
                    dat =
                        180.0 * atan2(outdata_y[ii], outdata_x[ii]) / 3.141592;
                    outdata_x[ii] = dat;
                }
                break;
            }
    }                           /* END switch */

    return (outdata_x);
    
}





int dmimggrad(void)
{

    char infile[DS_SZ_PATHNAME];
    char outfile[DS_SZ_PATHNAME];
    char operator[30];
    char output_val[10];
    short clobber;
    short verbose;

    void *data;
    long *lAxes;
    regRegion *dss;
    long null;
    short has_null;
    short *mask;
    dmDataType dt;
    dmBlock *inBlock;
    dmDescriptor *xdesc, *ydesc;

    dmBlock *outBlock;
    dmDescriptor *outDesc;

    double *outdata_x;
    double *outdata_y;

    Direction dir;
    Shapes shape;
    double kwidth;

    /* Get the parameters */
    clgetstr("infile", infile, DS_SZ_FNAME);
    clgetstr("outfile", outfile, DS_SZ_FNAME);
    clgetstr("gradient", operator, 30);
    clgetstr("value", output_val, 10);
    kwidth = clgetd("width");
    clobber = clgetb("clobber");
    verbose = clgeti("verbose");


    /* Boiler plate */
    if (NULL == (inBlock = dmImageOpen(infile))) {
        err_msg("ERROR: Cannot open image '%s'\n", infile);
        return (-1);
    }
    if (dmUNKNOWNTYPE == (dt = get_image_data(inBlock, &data, &lAxes,
                                              &dss, &null, &has_null))) {
        err_msg
            ("ERROR: Cannot get image data or unknown image data-type for "
             "file '%s'\n", infile);
        return (-1);
    }
    get_image_wcs(inBlock, &xdesc, &ydesc);
    mask = get_image_mask(inBlock, data, dt, lAxes, dss, null, has_null,
                          xdesc, ydesc);


    if (ds_clobber(outfile, clobber, NULL) != 0) {
        return (-1);
    }

    dir = get_dir( output_val );
    shape = get_shape( operator );
    if ( (GAUSSIAN == shape) && (INDEFD == kwidth)) {
        err_msg("ERRROR: Please specify the width for this gradient");
        return (-1);
    }
    

    /* Compute the X-gradient (if needed) */
    if (Y_GRAD != dir) {
        outdata_x = slide_convovle( data, dt, lAxes, mask, shape, kwidth, X_GRAD );
    } else {                    /* end X-dir */
        outdata_x = NULL;
    }

    /* Okay, now repeat for Y direction with some special cases */
    if (X_GRAD != dir) {
        if ((strcmp(operator, "laplace") == 0) && (outdata_x)) {
            outdata_y = outdata_x;
        } else {
            outdata_y = slide_convovle( data, dt, lAxes, mask, shape, kwidth, Y_GRAD );
        }
    } else {
        outdata_y = NULL;
    }
        
    double *outdata = combine_gradients( dir, outdata_x, outdata_y, lAxes );


    if (NULL == (outBlock = dmImageCreate(outfile, dmDOUBLE, lAxes, 2))) {
        err_msg("ERROR: Cannot create output image '%s'\n", outfile);
        return (-1);
    }
    outDesc = dmImageGetDataDescriptor(outBlock);
    dmBlockCopy(inBlock, outBlock, "HEADER");
    ds_copy_full_header(inBlock, outBlock, "dmimggrad", 0);
    put_param_hist_info(outBlock, "dmimggrad", NULL, 0);
    dmBlockCopyWCS(inBlock, outBlock);
    dmSetArray_d(outDesc, outdata, (lAxes[0] * lAxes[1]));
    dmImageClose(outBlock);
    dmImageClose(inBlock);

    if ( outdata_x ) {
        free( outdata_x );
    }
    if ( ( NULL != outdata_y ) && ( outdata_x != outdata_y )) {
        free( outdata_y);
    }
    
    return (0);

}
