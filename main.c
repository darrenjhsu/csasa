
#include <stdio.h>
#include <stdlib.h>
#include "mol_param.c"
#include <math.h>
#define PI 3.14159265359 

float dist(float* coord, int ii, int jj) {
    return sqrt((coord[3*ii] - coord[3*jj])*(coord[3*ii] - coord[3*jj])+
            (coord[3*ii+1] - coord[3*jj+1])*(coord[3*ii+1] - coord[3*jj+1])+
            (coord[3*ii+2] - coord[3*jj+2])*(coord[3*ii+2] - coord[3*jj+2]));
}



float calc_V(float* vdW, int* Ele, float* coord, float sol_s, int num_atom) {
    int num_raster = 256;
    int* closelist;
    int* closeidxat;
    float* V; 
    closelist = (int *)malloc(200 * num_atom *sizeof(int));
    closeidxat = (int *)malloc(num_atom * sizeof(int));
    V = (float *)malloc(num_atom * sizeof(float));
    for (int ii = 0; ii < num_atom; ii++) {
        closeidxat[ii] = 0;
        V[ii] = 0.0;
        for (int jj = 0; jj < 200; jj++) {
            closelist[ii * 200 + jj] = 0;
        }
    }
    /*for (int ii = 0; ii < num_atom; ii++) {
        printf("Fucking V[%d] = %.3f\n",ii, V[ii]);
    }
    */
    for (int ii = 0; ii < num_atom; ii++) {
        for (int jj = ii + 1; jj < num_atom; jj++) {
            if (dist(coord, ii, jj) < 5.4) {
                closelist[ii*200+closeidxat[ii]] = jj;
                closelist[jj*200+closeidxat[jj]] = ii;
                closeidxat[ii] ++;
                //printf("%d ",closeidxat[ii]);
                closeidxat[jj] ++;
                //printf("atom %d and atom %d are close\n", ii, jj);
            }
        }
    }
    /* 
    for (int ii = 0; ii < num_atom; ii++) {
        printf("Atom %d is close to %d neighbors\n", ii, closeidxat[ii]);
        printf("Atom %d close neighbors: ", ii);
        for (int jj = 0; jj < closeidxat[ii]; jj++) {
            printf("%d ", closelist[ii*200 + jj]);
        }
        printf("\n");
    }
    */
    for (int ii = 0; ii < num_atom; ii++) {
        int goodpoints = 0;
        int atom1t = Ele[ii];
        if (atom1t > 5) atom1t = 0;
        float vdW_s = vdW[atom1t];
        //printf("vdW = %.3f\n", vdW_s);
        float L = sqrt(num_raster * PI);
        float r = vdW_s + sol_s;
        // For 1 to probe points
        for (int jj = 0; jj < num_raster; jj++) {
            int pt = 1; 
            float h = 1.0 - (2.0 * (float)jj + 1.0) / (float)num_raster;
            float p = acos(h);
            float t = L * p; 
            float xu = sin(p) * cos(t);
            float yu = sin(p) * sin(t);
            float zu = cos(p);
            // vdW points
            float x = vdW_s * xu + coord[3*ii];
            float y = vdW_s * yu + coord[3*ii+1];
            float z = vdW_s * zu + coord[3*ii+2];
            // Solvent center
            float x2 = r * xu + coord[3*ii];
            float y2 = r * yu + coord[3*ii+1];
            float z2 = r * zu + coord[3*ii+2];
            /*if (ii == 0) {
                
                    printf("r = %.3f, xu = %.3f, x2 = %.3f, coord = %.3f\n",r, xu, x2, coord[3*ii]);
                
            }*/
            for (int kk = 0; kk < closeidxat[ii]; kk++) {
                //printf("%d, %d, %d\n", ii,jj,kk);
                int atom2i = closelist[ii * 200 + kk];
                int atom2t = Ele[atom2i];
                if (atom2t > 5) atom2t = 0;
                float dx = (x - coord[3*atom2i]);
                float dy = (y - coord[3*atom2i+1]);
                float dz = (z - coord[3*atom2i+2]);
                float dx2 = (x2 - coord[3*atom2i]);
                float dy2 = (y2 - coord[3*atom2i+1]);
                float dz2 = (z2 - coord[3*atom2i+2]);
                float dr2 = dx * dx + dy * dy + dz * dz; 
                float dr22 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
                // vdW points must not cross into other atom
                if (dr2 < vdW[atom2t] * vdW[atom2t]) {
                    pt = 0; //pts[jj] = 0;
                    break;
                }
                // solvent center has to be far enough
                if (dr22 < (vdW[atom2t]+sol_s) * (vdW[atom2t]+sol_s)) {
                    pt = 0; //pts[jj] = 0;
                    break;
                } 
            }
            goodpoints += pt;
            //printf("%d\n",goodpoints);
            // Generate a probe point (coord of that)
            // For nearby atoms
                // Calculate distance, if too close then break
                // If survived then add to V list (int array of num_atom)
        }
        //printf("Atom %d, goodpoints = %d, fraction = %.3f\n", ii, goodpoints, (float)goodpoints/(float)num_raster);
        V[ii]=(float)goodpoints/(float)num_raster;

//        printf("%.3f \n",(float)goodpoints / (float)num_raster);
        //printf("Vii = %.3f\n",V[ii]);
        //printf("Atom %d, V = %.3f\n", ii, V[ii]);
    }

    // Calculate surface area
    float vol = 0.0;
    for (int ii = 0; ii < num_atom; ii++) {
        int atomi = Ele[ii];
        if (atomi > 5) atomi = 0;
        vol += V[ii] * 4.0 * PI * (vdW[atomi]+sol_s) * (vdW[atomi]+sol_s);
    }

    printf("SASA is: %.3f\n", vol);
    return vol;
}

int main(int argc, char* argv[]) {
                  //   H     C     N     O     S    Fe water
    float vdW[7] = {1.07, 1.58, 0.84, 1.30, 1.68, 1.24, 1.67};
    int frames_total = atoi(argv[2]);
    int size_coord = 3 * num_atom * sizeof(float);
    float* coord;
    coord = (float *)malloc(size_coord);

    float sol_s = atof(argv[3]);
    //printf("sol_s = %.3f\n", sol_s);
    char* buf[100], buf1[100], buf2[100], buf3[100];
    float f1, f2, f3;
   

    printf("Num atoms: %d\n", num_atom);
    /* 
    for (int ii = 0; ii < num_atom; ii++) {
        printf("atomtype: %d\n", Ele[ii]);
    }
    */

    FILE *fp = fopen(argv[1],"r");
    if (fp == NULL) {
        printf("Opening file failed.\n");
        return 1;
    } else {
        printf("Opened file.\n");
    }

    FILE *fwp = fopen(argv[4],"w");
    if (fp == NULL) {
        printf("Opening file failed.\n");
        return 1;
    } else {
        printf("Opened file.\n");
    }
    // Read file by num_atom
    for (int ii = 0; ii < frames_total; ii++) {
        fscanf(fp,"%*s",buf);
        fscanf(fp,"%*s %d",buf);
        printf("Read the first two lines, ii = %d\n", ii);
        for (int jj = 0; jj < num_atom; jj++) {
            fscanf(fp,"%s %f %f %f", buf, &f1, &f2, &f3);
            //printf("Readed line %d\n", jj);
            coord[3*jj] = f1;
            coord[3*jj+1] = f2;
            coord[3*jj+2] = f3;
            //printf("Coord[jj] = %.3f, Coord[jj+1] = %.3f, Coord[jj+2] = %.3f\n",coord[3*jj], coord[3*jj+1], coord[3*jj+2]);
        }
        //printf("First coordinate is: %.3f\n",coord[0]);
        //printf("Final coordinate is: %.3f\n",coord[3*num_atom-1]);
        
        // Calculate V
        float SASA = calc_V(vdW, Ele, coord, sol_s, num_atom);
        fprintf(fwp, "Frame %d, SASA is: %.3f\n", ii, SASA);
    }
    fclose(fp);
    fclose(fwp);
    return 0;
}
