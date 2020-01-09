#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include "wtime.h"
#include "definitions.h"
#include "energy_struct.h"
#include "cuda_runtime.h"
#include "solver.h"

#define THREADS 32
#define TAMBLOCK 4


using namespace std;

/**
* Kernel del calculo de la solvation. Se debe anadir los parametros 
*/
__global__ void escalculation (int atoms_r, int atoms_l, int nlig, float *rec_x_d, float *rec_y_d, float *rec_z_d, float *lig_x_d, float *lig_y_d, float *lig_z_d, float *ql_d,float *qr_d, float *energy_d, int nconformations){

  double dist, total_elec = 0, miatomo[3], elecTerm;
  int totalAtomLig = nconformations * nlig;
  
  //#pragma omp parallel for reduction(+:total_elec), private (elecTerm, dist, miatomo)
	for (int k=0; k < totalAtomLig; k+=nlig){
	  for(int i=0;i<atoms_l;i++){					
			miatomo[0] = *(lig_x_d + k + i);
			miatomo[1] = *(lig_y_d + k + i);
			miatomo[2] = *(lig_z_d + k + i);

			for(int j=0;j<atoms_r;j++){				
				elecTerm = 0;
        dist=calculaDistancia (rec_x_d[j], rec_y_d[j], rec_z_d[j], miatomo[0], miatomo[1], miatomo[2]);
//				printf ("La distancia es %lf\n", dist);
        elecTerm = (ql_d[i]* qr_d[j]) / dist;
				total_elec += elecTerm;
//        printf ("La carga es %lf\n", total_elec);
			}
		}
		
		energy_d[k/nlig] = total_elec;
		total_elec = 0;
  }
	printf("Termino electrostatico %f\n", energy_d[0]);

	
  return;
}


/**
* Funcion para manejar el lanzamiento de CUDA 
*/
void forces_GPU_AU (int atoms_r, int atoms_l, int nlig, float *rec_x, float *rec_y, float *rec_z, float *lig_x, float *lig_y, float *lig_z, float *ql ,float *qr, float *energy, int nconformations){
	
	cudaError_t cudaStatus; //variable para recoger estados de cuda


	//seleccionamos device
	cudaSetDevice(0); //0 - Tesla K40 vs 1 - Tesla K230

	//creamos memoria para los vectores para GPU _d (device)
	float *rec_x_d, *rec_y_d, *rec_z_d, *qr_d, *lig_x_d, *lig_y_d, *lig_z_d, *ql_d, *energy_d;

	//reservamos memoria para GPU
  cudaMalloc(&rec_x_d, (int)sizeof(rec_x));
  cudaMalloc(&rec_y_d, (int)sizeof(rec_y));
  cudaMalloc(&rec_z_d, (int)sizeof(rec_z));
  cudaMalloc(&qr_d, (int)sizeof(qr));
  cudaMalloc(&lig_x_d, (int)sizeof(lig_x));
  cudaMalloc(&lig_y_d, (int)sizeof(lig_y));
  cudaMalloc(&lig_z_d, (int)sizeof(lig_z));
  cudaMalloc(&ql_d, (int)sizeof(ql));
  cudaMalloc(&energy_d, (int)sizeof(energy));
  printf("\n memoria para GPU reservada");	
 
  //pasamos datos de host to device
	cudaMemcpy(rec_x_d, rec_x, sizeof(rec_x), cudaMemcpyHostToDevice);
  cudaMemcpy(rec_y_d, rec_y, sizeof(rec_y), cudaMemcpyHostToDevice);
  cudaMemcpy(rec_z_d, rec_z, sizeof(rec_z), cudaMemcpyHostToDevice);
  cudaMemcpy(qr_d, qr, sizeof(qr), cudaMemcpyHostToDevice);
  cudaMemcpy(lig_x_d, lig_x, sizeof(lig_x), cudaMemcpyHostToDevice);
  cudaMemcpy(lig_y_d, lig_y, sizeof(lig_y), cudaMemcpyHostToDevice);
  cudaMemcpy(lig_z_d, lig_z, sizeof(lig_z), cudaMemcpyHostToDevice);
  cudaMemcpy(ql_d, ql, sizeof(ql), cudaMemcpyHostToDevice);
  cudaMemcpy(energy_d, energy, sizeof(energy), cudaMemcpyHostToDevice);
  printf("Datos pasados de host a device");
 
	//Definir numero de hilos y bloques
  dim3 block = THREADS/TAMBLOCK;
  dim3 thread = TAMBLOCK;
  printf("bloques: %d\n", THREADS/TAMBLOCK);
	printf("hilos por bloque: %d\n", TAMBLOCK);
  

	//llamamos a kernel
	escalculation <<< block,thread>>> (atoms_r, atoms_l, nlig, rec_x_d, rec_y_d, rec_z_d, lig_x_d, lig_y_d, lig_z_d, ql_d, qr_d, energy_d, nconformations);
	printf("kernel finalizado");
	//control de errores kernel
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if(cudaStatus != cudaSuccess) fprintf(stderr, "Error en el kernel %d\n", cudaStatus); 

	//Traemos info al host
  cudaMemcpy(rec_x, rec_x_d, sizeof(rec_x), cudaMemcpyDeviceToHost);
  cudaMemcpy(rec_y, rec_y_d, sizeof(rec_y), cudaMemcpyDeviceToHost);
  cudaMemcpy(rec_z, rec_z_d, sizeof(rec_z), cudaMemcpyDeviceToHost);
  cudaMemcpy(qr, qr_d, sizeof(qr), cudaMemcpyDeviceToHost);
  cudaMemcpy(lig_x, lig_x_d, sizeof(lig_x), cudaMemcpyDeviceToHost);
  cudaMemcpy(lig_y, lig_y_d, sizeof(lig_y), cudaMemcpyDeviceToHost);
  cudaMemcpy(lig_z, lig_z_d, sizeof(lig_z), cudaMemcpyDeviceToHost);
  cudaMemcpy(ql, ql_d, sizeof(ql), cudaMemcpyDeviceToHost);
  cudaMemcpy(energy, energy_d, sizeof(energy), cudaMemcpyHostToDevice);
	
  // para comprobar que la ultima conformacion tiene el mismo resultado que la primera
	printf("Termino electrostatico de conformacion %d es: %f\n", nconformations-1, energy[nconformations-1]); 

	//resultado varia repecto a SECUENCIAL y CUDA en 0.000002 por falta de precision con float
	//posible solucion utilizar double, probablemente bajara el rendimiento -> mas tiempo para calculo
	printf("Termino electrostatico %f\n", energy[0]);

	//Liberamos memoria reservada para GPU
}

/**
* Distancia euclidea compartida por funcion CUDA y CPU secuencial
*/
__device__ __host__ extern float calculaDistancia (float rx, float ry, float rz, float lx, float ly, float lz) {

  float difx = rx - lx;
  float dify = ry - ly;
  float difz = rz - lz;
  float mod2x=difx*difx;
  float mod2y=dify*dify;
  float mod2z=difz*difz;
  difx=mod2x+mod2y+mod2z;
  return sqrtf(difx);
}




/**
 * Funcion que implementa el termino electrostático en CPU
 */
void forces_CPU_AU (int atoms_r, int atoms_l, int nlig, float *rec_x, float *rec_y, float *rec_z, float *lig_x, float *lig_y, float *lig_z, float *ql ,float *qr, float *energy, int nconformations){

	double dist, total_elec = 0, miatomo[3], elecTerm;
  int totalAtomLig = nconformations * nlig;

	for (int k=0; k < totalAtomLig; k+=nlig){
	  for(int i=0;i<atoms_l;i++){					
			miatomo[0] = *(lig_x + k + i);
			miatomo[1] = *(lig_y + k + i);
			miatomo[2] = *(lig_z + k + i);

			for(int j=0;j<atoms_r;j++){				
				elecTerm = 0;
        dist=calculaDistancia (rec_x[j], rec_y[j], rec_z[j], miatomo[0], miatomo[1], miatomo[2]);
//				printf ("La distancia es %lf\n", dist);
        elecTerm = (ql[i]* qr[j]) / dist;
				total_elec += elecTerm;
//        printf ("La carga es %lf\n", total_elec);
			}
		}
		
		energy[k/nlig] = total_elec;
		total_elec = 0;
  }
	printf("Termino electrostatico %f\n", energy[0]);
}


extern void solver_AU(int mode, int atoms_r, int atoms_l,  int nlig, float *rec_x, float *rec_y, float *rec_z, float *lig_x, float *lig_y, float *lig_z, float *ql, float *qr, float *energy_desolv, int nconformaciones) {

	double elapsed_i, elapsed_o;
	
	switch (mode) {
		case 0://Sequential execution
			printf("\* CALCULO ELECTROSTATICO EN CPU *\n");
			printf("**************************************\n");			
			printf("Conformations: %d\t Mode: %d, CPU\n",nconformaciones,mode);			
			elapsed_i = wtime();
			forces_CPU_AU (atoms_r,atoms_l,nlig,rec_x,rec_y,rec_z,lig_x,lig_y,lig_z,ql,qr,energy_desolv,nconformaciones);
			elapsed_o = wtime() - elapsed_i;
			printf ("CPU Processing time: %f (seg)\n", elapsed_o);
			break;
		case 1: //OpenMP execution
			printf("\* CALCULO ELECTROSTATICO EN OPENMP *\n");
			printf("**************************************\n");			
			printf("**************************************\n");			
			printf("Conformations: %d\t Mode: %d, CMP\n",nconformaciones,mode);			
			elapsed_i = wtime();
			forces_OMP_AU (atoms_r,atoms_l,nlig,rec_x,rec_y,rec_z,lig_x,lig_y,lig_z,ql,qr,energy_desolv,nconformaciones);
			elapsed_o = wtime() - elapsed_i;
			printf ("OpenMP Processing time: %f (seg)\n", elapsed_o);
			break;
		case 2: //CUDA exeuction
			printf("\* CALCULO ELECTROSTATICO EN CUDA *\n");
      printf("**************************************\n");
      printf("Conformaciones: %d\t Mode: %d, GPU\n",nconformaciones,mode);
			elapsed_i = wtime();
			forces_GPU_AU (atoms_r,atoms_l,nlig,rec_x,rec_y,rec_z,lig_x,lig_y,lig_z,ql,qr,energy_desolv,nconformaciones);
			elapsed_o = wtime() - elapsed_i;
			printf ("GPU Processing time: %f (seg)\n", elapsed_o);			
			break; 	
	  	default:
 	    	printf("Wrong mode type: %d.  Use -h for help.\n", mode);
			exit (-1);	
	} 		
}
