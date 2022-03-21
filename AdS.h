#pragma once

double AdSrho(double a){
	double b;
	b = a/(1+a);
	return b;
	}
	
double AdSq(double a){
	double q;
	q = 1-a;
	return q;
	}

double AdSa( double b, double c){
	double a;
	a = b*b + c*c;
	return a;
	}

double DAdSa(double a){
	double b;
	b = -2 + 4*a;
	return b;
	}
