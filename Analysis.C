//#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "ana112.C"

int Analysis(int id){
	ana112 t(id);
	t.Loop();
	cout << endl;
	cout << "complete Loop()" << endl;
	cout << endl;

	t.Ped_mip();
	cout << endl;
	cout << "complete Ped_mip()" << endl;
	cout << endl;

	for(int i=0; i<6; i++){
		t.Slewing(i,10);
		t.Resolution(i);
		cout << endl;
		cout << "complete ana det"<< i << endl;
		cout << endl;
	}
	t.Save_resolution();
	return 0;
}
