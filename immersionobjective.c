#include <windows.h>
#include <math.h>
#include <string.h>
#include "usersurf.h"

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD);
int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn);

BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved) {
    return TRUE;
}

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD) {
    int i;
    double p2, alpha, power, a, b, c, rad, casp, t, zc;
    double f_length, n_imm, cosphi, r0x, r0y, r0z, scaling, x_vert, y_vert;

    switch(FD->type) {
        case 0:
      	    switch(FD->numb) {
                case 0:
                    // Surface name
		            strcpy_s(UD->string, 12, "ImmersionOL");
                    break;
                case 1:
                    // Rotational symmetry
                    strcpy_s(UD->string, 2, "1");
                    break;
                case 2:
                    // Gradient-index
                    UD->string[0] = '\0';
                    break;
            }
            break;
        case 1:
            switch(FD->numb) {
		        case 1:
            	    strcpy_s(UD->string, 28, "Objective lens focal length");
            	    break;
				case 2:
            	    strcpy_s(UD->string, 34, "Immersion medium refractive index");
                    break;
                default:
                    UD->string[0] = '\0';
            	    break;
            }
      	    break;
        case 2:
            switch(FD->numb) {
                default:
            	    UD->string[0] = '\0';
            	    break;
            }
      	    break;
        case 3:
            // Surface sag (representation)
            UD->sag1 = 0.0;
            UD->sag2 = 0.0;

			if (FD->cv == 0) return(0);
         
            p2 = UD->x * UD->x + UD->y * UD->y;
            alpha = 1 - (1+FD->k)*FD->cv*FD->cv*p2;

        		 if (fabs(alpha) < 1e-13)
		        alpha = 0;

            if (alpha < 0) return(-1);
            
            UD->sag1 = (FD->cv*p2)/(1 + sqrt(alpha));
            
            if (alpha != 1.0) UD->sag2 = (FD->cv*p2)/(1 - sqrt(alpha));
      	    break;
        case 4:
      	    // Paraxial ray trace (modified to keep consistent optical power, rest is broken)
            UD->ln =  0.0;
            UD->mn =  0.0;
            UD->nn = -1.0;
            
            // Custom part starts here
            f_length = FD->param[1];
            power = 1/f_length;
            // Custom part ends here

            if ((UD->n) != 0.0)	{
                (UD->l) = (UD->l)/(UD->n);
                (UD->m) = (UD->m)/(UD->n);

                (UD->l) = (FD->n1*(UD->l) - (UD->x)*power)/(FD->n2);
                (UD->m) = (FD->n1*(UD->m) - (UD->y)*power)/(FD->n2);

                (UD->n) = sqrt(1/(1 + (UD->l)*(UD->l) + (UD->m)*(UD->m) ) );

                (UD->l) = (UD->l)*(UD->n);
                (UD->m) = (UD->m)*(UD->n);
            }
            break;
        case 5:
      	    // Real ray trace (modified according to S-U Hwang, and Y-G Lee in 2008)
            if (FD->cv == 0.0)	{
	            UD->ln =  0.0;
   	            UD->mn =  0.0;
      	        UD->nn = -1.0;
			    
                if (Refract(FD->n1, FD->n2, &UD->l, &UD->m, &UD->n, UD->ln, UD->mn, UD->nn)) return(-FD->surf);
                
                return(0);
            }

            a = (UD->n) * (UD->n) * FD->k + 1;
            b = ((UD->n) / FD->cv) - (UD->x) * (UD->l) - (UD->y) * (UD->m);
            c = (UD->x) * (UD->x) + (UD->y) * (UD->y);
            rad = b * b - a * c;
                
            if (rad < 0) return(FD->surf);

            if (FD->cv > 0) t = c / (b + sqrt(rad));
            else           t = c / (b - sqrt(rad));
            
            // Custom part starts here
            x_vert = UD->x;
            y_vert = UD->y;
            // Custom part ends here

            (UD->x) = (UD->l) * t + (UD->x);
            (UD->y) = (UD->m) * t + (UD->y);
            (UD->z) = (UD->n) * t + (UD->z);
            UD->path = t;

            // Custom part starts here
            f_length = FD->param[1];
            n_imm = FD->param[2];
            cosphi = UD->n;

            scaling = f_length / cosphi - t;

            r0x = scaling * UD->l - x_vert;
            r0y = scaling * UD->m - y_vert;
            r0z = scaling * UD->n + f_length * (n_imm - 1);

            UD->l = r0x / sqrt(r0x * r0x + r0y * r0y + r0z * r0z);
            UD->m = r0y / sqrt(r0x * r0x + r0y * r0y + r0z * r0z);
            UD->n = r0z / sqrt(r0x * r0x + r0y * r0y + r0z * r0z);
            // Custom part ends here
            break;
        case 6:
            UD->index = FD->n2;
            UD->dndx = 0.0;
            UD->dndy = 0.0;
            UD->dndz = 0.0;
      	    break;
        case 7:
            FD->param[1] = 1.8;
            FD->param[2] = 1.515;
            break;
        case 8:
            break;
        }

        return 0;
}

int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn) {
    double nr, cosi, cosi2, rad, cosr, gamma;
    
    if (thisn != nextn) {
	    nr = thisn / nextn;
	    cosi = fabs((*l) * ln + (*m) * mn + (*n) * nn);
	    cosi2 = cosi * cosi;

	    if (cosi2 > 1) cosi2 = 1;

	    rad = 1 - ((1 - cosi2) * (nr * nr));

	    if (rad < 0) return(-1);

	    cosr = sqrt(rad);
	    gamma = nr * cosi - cosr;
	    (*l) = (nr * (*l)) + (gamma * ln);
	    (*m) = (nr * (*m)) + (gamma * mn);
	    (*n) = (nr * (*n)) + (gamma * nn);
	}

        return 0;
}

