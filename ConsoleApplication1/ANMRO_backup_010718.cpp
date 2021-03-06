#include "stdafx.h"
using namespace std;


int veloX(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);
int veloY(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);
int veloZ(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);

int fx(Ipp64f * field, Ipp64f *vy, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out);
int fy(Ipp64f * field, Ipp64f *vx, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out);
int fz(Ipp64f * field, Ipp64f *vx, Ipp64f *vy, int length, Ipp64f *temp, Ipp64f *out);

int _tmain(int argc, _TCHAR* argv[])
{
	Ipp64f field[3] = { 4.60929296230880, 2.149348607004536, 6.06101097852378 };
	Ipp64f tau = 0.1;
	Ipp64f final = 10;
    long steps = 1000;
	Ipp64f h = final / steps;
	
	
	DataExtractor extractor("starts.dat");
	Ipp64f * starts = extractor.getDataArray();
	int nPoints = floor((extractor.getNumberOfLines())/3);		
	
	Ipp64f *output = new Ipp64f[nPoints*steps*3];

	Ipp64f *times = new Ipp64f[steps]; //time steps
	for (int i = 0; i < steps; i++) {
		times[i] = - i;
	}
	ippsMulC_64f_I(h, times, steps);
	
	std::clock_t startT;
	Ipp64f duration;

			
		Ipp64f *vzStorage = new Ipp64f[steps*nPoints];
		Ipp64f *vz0Storage = new Ipp64f[nPoints];
		Ipp64f *DOS = new Ipp64f[nPoints];
		Ipp64f *ones = new Ipp64f[nPoints];//for inverting
		ippsSet_64f(1, ones, nPoints);

		Ipp64f *vx = new Ipp64f[nPoints];
		Ipp64f *vy = new Ipp64f[nPoints];
		Ipp64f *vz = new Ipp64f[nPoints];

		Ipp64f *argx = new Ipp64f[nPoints];
		Ipp64f *argy = new Ipp64f[nPoints];
		Ipp64f *argz = new Ipp64f[nPoints];

		Ipp64f *tempx = new Ipp64f[2*nPoints];
		Ipp64f *tempy = new Ipp64f[2*nPoints];
		Ipp64f *tempz = new Ipp64f[2*nPoints];

		Ipp64f *k1x = new Ipp64f[nPoints];
		Ipp64f *k1y = new Ipp64f[nPoints];
		Ipp64f *k1z = new Ipp64f[nPoints];		
		Ipp64f *k2x = new Ipp64f[nPoints];
		Ipp64f *k2y = new Ipp64f[nPoints];
		Ipp64f *k2z = new Ipp64f[nPoints];
		Ipp64f *k3x = new Ipp64f[nPoints];
		Ipp64f *k3y = new Ipp64f[nPoints];
		Ipp64f *k3z = new Ipp64f[nPoints];
		Ipp64f *k4x = new Ipp64f[nPoints];
		Ipp64f *k4y = new Ipp64f[nPoints];
		Ipp64f *k4z = new Ipp64f[nPoints];

		Ipp64f total = 0;

			startT = std::clock();

			for (int j = 0; j < nPoints; j++) {
				output[0 * nPoints + j] = starts[j * 3 + 0];
				output[1 * nPoints + j] = starts[j * 3 + 1];
				output[2 * nPoints + j] = starts[j * 3 + 2];
			}

			for (int i = 1; i < steps; i++) {

				ippsCopy_64f(&output[nPoints * (3*(i - 1) + 0)],argx,nPoints);//copy arguments for k1;
				ippsCopy_64f(&output[nPoints * (3*(i - 1) + 1)],argy,nPoints);
				ippsCopy_64f(&output[nPoints * (3*(i - 1) + 2)],argz,nPoints);
				veloX(argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(argx, argy, argz, nPoints, tempy, vy);
				veloZ(argx, argy, argz, nPoints, tempz, vz);	
				ippsCopy_64f(vz, &vzStorage[nPoints * (i-1)], nPoints);//store vz for conductivity later
				
				fx(field, vy, vz, nPoints, tempx, k1x); //calculate evolution in k and store in k1
				fy(field, vx, vz, nPoints, tempy, k1y);
				fz(field, vx, vy, nPoints, tempz, k1z);
								
				ippsMulC_64f(k1x,h/2,tempx,nPoints); //prep evolved k step for k2
				ippsMulC_64f(k1y,h/2,tempy,nPoints);
				ippsMulC_64f(k1z,h/2,tempz,nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3*(i - 1) + 0)], argx,nPoints); //add step to previous k point, load into arguments for k2;
				ippsAdd_64f(tempy, &output[nPoints * (3*(i - 1) + 1)], argy,nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3*(i - 1) + 2)], argz,nPoints);
				veloX(argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(argx, argy, argz, nPoints, tempy, vy);
				veloZ(argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k2x); //calculate evolution in k and store in k2
				fy(field, vx, vz, nPoints, tempy, k2y);
				fz(field, vx, vy, nPoints, tempz, k2z);
				
				ippsMulC_64f(k2x,h/2,tempx,nPoints); //prep evolved k step for k3
				ippsMulC_64f(k2y,h/2,tempy,nPoints);
				ippsMulC_64f(k2z,h/2,tempz,nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3*(i - 1) + 0)], argx,nPoints); //add step to previous k point, load into arguments for k3;
				ippsAdd_64f(tempy, &output[nPoints * (3*(i - 1) + 1)], argy,nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3*(i - 1) + 2)], argz,nPoints);
				veloX(argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(argx, argy, argz, nPoints, tempy, vy);
				veloZ(argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k3x); //calculate evolution in k and store in k3
				fy(field, vx, vz, nPoints, tempy, k3y);
				fz(field, vx, vy, nPoints, tempz, k3z);

				ippsMulC_64f(k3x,h,tempx,nPoints); //prep evolved k step for k4
				ippsMulC_64f(k3y,h,tempy,nPoints);
				ippsMulC_64f(k3z,h,tempz,nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3*(i - 1) + 0)], argx,nPoints); //add step to previous k point, load into arguments for k4;
				ippsAdd_64f(tempy, &output[nPoints * (3*(i - 1) + 1)], argy,nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3*(i - 1) + 2)], argz,nPoints);
				veloX(argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(argx, argy, argz, nPoints, tempy, vy);
				veloZ(argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k4x); //calculate evolution in k and store in k4
				fy(field, vx, vz, nPoints, tempy, k4y);
				fz(field, vx, vy, nPoints, tempz, k4z);

				ippsMulC_64f_I(2,k2x,nPoints); //scale k2
				ippsMulC_64f_I(2,k2y,nPoints);
				ippsMulC_64f_I(2,k2z,nPoints);
				ippsMulC_64f_I(2,k3x,nPoints); //scale k3
				ippsMulC_64f_I(2,k3y,nPoints);
				ippsMulC_64f_I(2,k3z,nPoints);

				ippsAdd_64f(k1x,k2x,tempx,nPoints); //add k1 + k2 to temp
				ippsAdd_64f(k1y,k2y,tempy,nPoints);
				ippsAdd_64f(k1z,k2z,tempz,nPoints);

				ippsAdd_64f_I(k3x,tempx,nPoints); //add in k3
				ippsAdd_64f_I(k3y,tempy,nPoints);
				ippsAdd_64f_I(k3z,tempz,nPoints);

				ippsAdd_64f_I(k4x,tempx,nPoints); //add in k4
				ippsAdd_64f_I(k4y,tempy,nPoints);
				ippsAdd_64f_I(k4z,tempz,nPoints);

				ippsMulC_64f_I(h/6,tempx,nPoints); //scale the entire sum
				ippsMulC_64f_I(h/6,tempy,nPoints); //scale the entire sum
				ippsMulC_64f_I(h/6,tempz,nPoints); //scale the entire sum

				ippsAdd_64f(&output[nPoints * (3*(i-1) + 0)],tempx,&output[nPoints * (3*i + 0)],nPoints); //add sum to previous output and store
				ippsAdd_64f(&output[nPoints * (3*(i-1) + 1)],tempy,&output[nPoints * (3*i + 1)],nPoints);
				ippsAdd_64f(&output[nPoints * (3*(i-1) + 2)],tempz,&output[nPoints * (3*i + 2)],nPoints);
			}
			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 0)], argx, nPoints);//get velocity for last point
			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 1)], argy, nPoints);
			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 2)], argz, nPoints);
			veloZ(argx, argy, argz, nPoints, tempz, vz);
			ippsCopy_64f(vz, &vzStorage[nPoints * (steps - 1)], nPoints);

			/*ofstream fout2;
			fout2.open("outputVz.dat");
			fout2.precision(15);

			for (int i = 0; i < steps*nPoints; i++) {

				fout2 << vzStorage[i] << endl;
			}
			fout2.close();*/


			ippsCopy_64f(&output[nPoints * (0)], argx, nPoints);//initial velocities for DOS calc;
			ippsCopy_64f(&output[nPoints * (1)], argy, nPoints);
			ippsCopy_64f(&output[nPoints * (2)], argz, nPoints);
			veloX(argx, argy, argz, nPoints, tempx, vx); //velocities for DOS are stored in vx, vy, and vz buffers.
			veloY(argx, argy, argz, nPoints, tempy, vy);
			veloZ(argx, argy, argz, nPoints, tempz, vz);

			ippsSqr_64f_I(vx, nPoints);//in-place square of velocities
			ippsSqr_64f_I(vy, nPoints);
			ippsSqr_64f_I(vz, nPoints);

			ippsAdd_64f(vx, vy, tempx, nPoints);//add all square velocities
			ippsAdd_64f_I(vz, tempx, nPoints);
			ippsSqrt_64f_I(tempx, nPoints);//square root
			ippsDiv_64f(tempx, ones, DOS, nPoints);

			ippsDivC_64f_I(tau, times, steps);//exponential stuff, negative tau is taken care of in time
			ippsExp_64f_I(times, steps);
			ippsMulC_64f_I((1E-12)*h, times, steps);
						

			ippsCopy_64f(&vzStorage[0],vz0Storage, nPoints);//save initial velocity before exp

			for (int i = 0; i < steps; i++) {
				ippsMulC_64f_I(times[i], &vzStorage[i*nPoints], nPoints); //multiply velocities by exp time factor
			}
			
			for (int i = 0; i < (steps-1); i++) {
				ippsAdd_64f_I(&vzStorage[i*nPoints],&vzStorage[(i+1)*nPoints],nPoints); //add all and accumulate in last vector
			}			
			
			ippsMul_64f_I(DOS, &vzStorage[(steps-1)*nPoints], nPoints);
			ippsMul_64f_I(vz0Storage, &vzStorage[(steps - 1)*nPoints], nPoints);//multiply by initial velocities

			ippsSum_64f(&vzStorage[(steps -1)*nPoints], nPoints, &total);//sum all elements of velocity vector

			cout << total << endl;
			
			
		duration = (std::clock() - startT) / (Ipp64f)CLOCKS_PER_SEC;
		cout << "time: " << duration << endl;

		
		/*ofstream fout;
		fout.open("output.dat");
		fout.precision(15);	

		for (int i = 0; i < nPoints*3*steps; i++) {

			fout << output[i] << endl;
		}

		fout.close();
		*/
	

		/*ofstream fout2;
		fout2.open("vztimes.dat");
		fout2.precision(15);

		for (int i = 0; i < steps; i++) {

			fout2 << times[i] << endl;
		}

		fout2.close();*/

		char a;
		cin >> a;
	
	
	delete[] output;
	return 0;
}

int veloX(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	/*return 5710 * sin(3.74767*kx);*/
	ippsMulC_64f(kx,3.74767,temp,length);
	vdSin(length,temp,out);
	ippsMulC_64f_I(5710,out,length);
	return 0;
}

int veloY(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	/*return 5710 * sin(3.74767*ky);*/
	ippsMulC_64f(ky,3.74767,temp,length);
	vdSin(length,temp,out);
	ippsMulC_64f_I(5710,out,length);
	return 0;
}

int veloZ(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	/*return  .10 * sin(3.3*kz);*/
	ippsMulC_64f(kz,3.3,temp,length);
	vdSin(length,temp,out);
	ippsMulC_64f_I(0.1,out,length);
	return 0;
}

int fx(Ipp64f * field, Ipp64f *vy, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out) {
	/*return -1 / (11538.5) * (vy * field[2] - vz*field[1]);*/
	ippsMulC_64f(vy,field[2],temp,length);
	ippsMulC_64f(vz,field[1],&temp[length],length);
	ippsMulC_64f_I(-1,&temp[length],length);
	ippsAdd_64f_I(&temp[length],temp,length);
	ippsMulC_64f(temp,1/(11538.5),out,length);//+1 to run back in time
	
	return 0;
	
}
int fy(Ipp64f * field, Ipp64f *vx, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out) {
	//return -1 / (11538.5) * (vz * field[0] - vx*field[2]);
	ippsMulC_64f(vz,field[0],temp,length);
	ippsMulC_64f(vx,field[2],&temp[length],length);
	ippsMulC_64f_I(-1,&temp[length],length);
	ippsAdd_64f_I(&temp[length],temp,length);
	ippsMulC_64f(temp,1/(11538.5),out,length);//+1 to run back in time
	
	return 0;
}
int fz(Ipp64f * field, Ipp64f *vx, Ipp64f *vy, int length, Ipp64f *temp, Ipp64f *out) {
	//return -1 / (11538.5) * (vx * field[1] - vy*field[0]);
	ippsMulC_64f(vx,field[1],temp,length);
	ippsMulC_64f(vy,field[0],&temp[length],length);
	ippsMulC_64f_I(-1,&temp[length],length);
	ippsAdd_64f_I(&temp[length],temp,length);
	ippsMulC_64f(temp,1/(11538.5),out,length);//+1 to run back in time
	
	return 0;
}

