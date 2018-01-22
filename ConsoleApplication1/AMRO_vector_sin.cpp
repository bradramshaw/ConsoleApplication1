#include "stdafx.h"
using namespace std;


int veloX(double *kx, double *ky, double *kz, int length, double *out);
int veloY(double *kx, double *ky, double *kz, int length, double *out);
int veloZ(double *kx, double *ky, double *kz, int length, double *out);

double fx(double * field, double vy, double vz);
double fy(double * field, double vx, double vz);
double fz(double * field, double vx, double vy);

double kx(double h, double * field, double vx, double vy, double vz);
double ky(double h, double * field, double vx, double vy, double vz);
double kz(double h, double * field, double vx, double vy, double vz);

int _tmain(int argc, _TCHAR* argv[])
{

	double field[3] = { 5.00853, 5.96893, 1.37392 };
	double final = 10;
    long steps = 3;
	double h = final / steps;

	DataExtractor extractor("starts.dat");
	double * starts = extractor.getDataArray();
	int nPoints = floor((extractor.getNumberOfLines())/3);	
	
	double *output = new double[nPoints*steps*3];
	

	std::clock_t startT;
	double duration;

	while (true) {
			
		startT = std::clock();

		for (int j = 0; j < nPoints; j++) {
			output[0 * nPoints + j] = starts[j * 3 + 0];
			output[1 * nPoints + j] = starts[j * 3 + 1];
			output[2 * nPoints + j] = starts[j * 3 + 2];
		}

			
		double *vx = new double[nPoints];
		double *vy = new double[nPoints];
		double *vz = new double[nPoints];

			for (int i = 1; i < steps; i++) {

				double *argx0 = &output[nPoints * (i - 1 + 0)];
				double *argy0 = &output[nPoints * (i - 1 + 1)];
				double *argz0 = &output[nPoints * (i - 1 + 2)];
							
				veloX(argx0, argy0, argz0, nPoints, vx);
				veloY(argx0, argy0, argz0, nPoints, vy);
				veloZ(argx0, argy0, argz0, nPoints, vz);

				//double k1x0 = kx(h, field, vx0, vy0, vz0);
				//double k1y0 = ky(h, field, vx0, vy0, vz0);
				//double k1z0 = kz(h, field, vx0, vy0, vz0);

				//double argx1 = output[3 * steps * j + 3 * (i - 1) + 0] + h*k1x0 / 2;
				//double argy1 = output[3 * steps * j + 3 * (i - 1) + 1] + h*k1y0 / 2;
				//double argz1 = output[3 * steps * j + 3 * (i - 1) + 2] + h*k1z0 / 2;
				//double vx1 = veloX(argx1, argy1, argz1);
				//double vy1 = veloY(argx1, argy1, argz1);
				//double vz1 = veloZ(argx1, argy1, argz1);
				//double k2x0 = kx(h, field, vx1, vy1, vz1);
				//double k2y0 = ky(h, field, vx1, vy1, vz1);
				//double k2z0 = kz(h, field, vx1, vy1, vz1);

				//double argx2 = output[3 * steps * j + 3 * (i - 1) + 0] + h*k2x0 / 2;
				//double argy2 = output[3 * steps * j + 3 * (i - 1) + 1] + h*k2y0 / 2;
				//double argz2 = output[3 * steps * j + 3 * (i - 1) + 2] + h*k2z0 / 2;
				//double vx2 = veloX(argx2, argy2, argz2);
				//double vy2 = veloY(argx2, argy2, argz2);
				//double vz2 = veloZ(argx2, argy2, argz2); 
				//double k3x0 = kx(h, field, vx2, vy2, vz2);
				//double k3y0 = ky(h, field, vx2, vy2, vz2);
				//double k3z0 = kz(h, field, vx2, vy2, vz2);


				//double argx3 = output[3 * steps * j + 3 * (i - 1) + 0] + h*k3x0;
				//double argy3 = output[3 * steps * j + 3 * (i - 1) + 1] + h*k3y0;
				//double argz3 = output[3 * steps * j + 3 * (i - 1) + 2] + h*k3z0;
				//double vx3 = veloX(argx3, argy3, argz3);
				//double vy3 = veloY(argx3, argy3, argz3);
				//double vz3 = veloZ(argx3, argy3, argz3); 
				//double k4x0 = kx(h, field, vx3, vy3, vz3);
				//double k4y0 = ky(h, field, vx3, vy3, vz3);
				//double k4z0 = kz(h, field, vx3, vy3, vz3);

				//output[3 * steps * j + 3 * i + 0] = output[3 * steps * j + 3 * (i - 1) + 0] + h / 6*(k1x0 + 2 * k2x0 + 2 * k3x0 + k4x0);
				//output[3 * steps * j + 3 * i + 1] = output[3 * steps * j + 3 * (i - 1) + 1] + h / 6*(k1y0 + 2 * k2y0 + 2 * k3y0 + k4y0);
				//output[3 * steps * j + 3 * i + 2] = output[3 * steps * j + 3 * (i - 1) + 2] + h / 6*(k1z0 + 2 * k2z0 + 2 * k3z0 + k4z0);

			}

		


		duration = (std::clock() - startT) / (double)CLOCKS_PER_SEC;
		cout << "time: " << duration << endl;

		ofstream fout;
		fout.open("vx.dat");
		fout.precision(15);	

		for (int i = 0; i < nPoints; i++) {

			fout << vx[i] << endl;
		}

		fout.close();

		
		char a;
		cin >> a;
	
	}
	delete[] output;
	return 0;
}

int veloX(double *kx, double *ky, double *kz, int length, double *out) {
	/*return 5710 * sin(3.74767*kx);*/
	vdSin(length,kx,out);
	return 0;
}

int veloY(double *kx, double *ky, double *kz, int length, double *out) {
	/*return 5710 * sin(3.74767*ky);*/
	return 0;
}

int veloZ(double *kx, double *ky, double *kz, int length, double *out) {
	/*return  .10 * sin(3.3*kz);*/
	return 0;
}

double fx(double * field, double vy, double vz) {
	return -1 / (11538.5) * (vy * field[2] - vz*field[1]);
	
}
double fy(double * field, double vx, double vz) {
	return -1 / (11538.5) * (vz * field[0] - vx*field[2]);

}
double fz(double * field, double vx, double vy) {
	return -1 / (11538.5) * (vx * field[1] - vy*field[0]);
}

double kx(double h, double * field, double vx, double vy, double vz) {
	return	fx(field, vy, vz);
}
double ky(double h, double * field, double vx, double vy, double vz) {
	return	fy(field, vx, vz);
}
double kz(double h, double * field, double vx, double vy, double vz) {
	return	fz(field, vx, vy);
}
