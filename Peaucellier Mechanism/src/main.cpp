//OpenCV
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

using namespace std;
using namespace cv;

const double PI = 3.141592;

int main(){

	float a = 180; //30
	float b = 390; //80
	float d = 110; //40

	float Qc = abs( ( (b*b) - (d*d) - 4*(a*a) )/(4*a) );
	float Rc = sqrt( (d*d) - (Qc*Qc) );

	float O4x = 450;
	float O4y = 350;
	float O2x = O4x + a;
	float O2y = O4y;
	float Qx = O4x + 2*a;
	float Qy = O4y;
	float Rx = O4x + 2*a + Qc;
	float Ry = O4y - Rc;
	float Sx = O4x + 2*a + Qc;
	float Sy = O4y + Rc;
	float Px = O4x + 2*a + 2*Qc;
	float Py = O4y;

	float Betamax = acos( abs(b-d)/(2*a) );

	float signal = 1;
	double teta = 0, ang = 2*PI/180;
	vector<float> pontosy;
	while(1){

		Mat simulation(700,1366,CV_8UC3,Scalar(255,255,255) );
		for(size_t i=0;i<pontosy.size();i++){
			circle(simulation,Point( int(Px),int(pontosy.at(i)) ),1,Scalar(0,0,0),-1,8,0);
		}
		line(simulation,Point(int(O2x),int(O2y)),Point(int(Qx),int(Qy)),Scalar(255,0,0),1,8,0);
		line(simulation,Point(int(Qx),int(Qy)),Point(int(Sx),int(Sy)),Scalar(0,0,255),1,8,0);
		line(simulation,Point(int(Sx),int(Sy)),Point(int(Px),int(Py)),Scalar(0,0,255),1,8,0);
		line(simulation,Point(int(Qx),int(Qy)),Point(int(Rx),int(Ry)),Scalar(0,0,255),1,8,0);
		line(simulation,Point(int(Rx),int(Ry)),Point(int(Px),int(Py)),Scalar(0,0,255),1,8,0);
		line(simulation,Point(int(O4x),int(O4y)),Point(int(Rx),int(Ry)),Scalar(0,255,0),1,8,0);
		line(simulation,Point(int(O4x),int(O4y)),Point(int(Sx),int(Sy)),Scalar(0,255,0),1,8,0);
		circle(simulation,Point( int(O4x),int(O4y) ),4,Scalar(0,0,0),-1,8,0);
		circle(simulation,Point( int(O2x),int(O2y) ),4,Scalar(0,0,0),-1,8,0);
		circle(simulation,Point( int(Qx),int(Qy) ),4,Scalar(0,0,0),1,8,0);
		circle(simulation,Point( int(Rx),int(Ry) ),4,Scalar(0,0,0),1,8,0);
		circle(simulation,Point( int(Sx),int(Sy) ),4,Scalar(0,0,255),1,8,0);
		circle(simulation,Point( int(Px),int(Py) ),4,Scalar(0,0,0),1,8,0);
		imshow("MECANISMO",simulation);
		waitKey(0);

		teta = teta + ang;
		float tempbeta = float((PI/2) - ( ( PI-abs(teta) )/2 ));
		if( tempbeta >= Betamax ){
			ang = (-1)*ang;
			teta = teta + ang;
		}

		float O4q = float(sqrt( 2*a*a - 2*a*a*cos( PI - abs(teta) ) ));
		float O4p = abs( ( (b*b) - (d*d) )/O4q );
		float Beta = abs(acos( O4q/(2*a) ));
		float Alfa = abs(acos( ( (O4q*O4q) + (b*b) - (d*d) )/( 2*O4q*b ) ));

		Qx = float(cos(teta) * a); Qx = Qx + O2x;
		Qy = float(sin(teta) * a); Qy = Qy + O2y;
		Rx = O4x + b*cos( Beta+Alfa );
		Sx = O4x + b*cos( Beta-Alfa );
		Px = O4x + O4p*cos(Beta);
		if( Beta==0 ){
			signal = signal*(-1);
		}
		Py = O4y + signal*O4p*sin(Beta);
		pontosy.push_back( Py );
		Ry = O4y + signal*b*sin( Beta+Alfa );
		Sy = O4y + signal*b*sin( Beta-Alfa );
	}
	return 0;
}
