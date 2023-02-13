//#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "ana112.C"

int Check_signal(int id){
	ana112 t(id);
	t.Loop();
	return 0;
}
