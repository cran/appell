#include <R.h> 
#include <Rinternals.h>

void printstringc(const char *label, long int length)
{
    if(length > 255) {
    	Rf_warning("printstringc: invalid string length");
    } else if(length > 0) {
    	for (long int k = 0; k < length; k++)
    	    Rprintf("%c", label[k]);
    }
}

/* the string length is passed silently by value from Fortran! */
void F77_SUB(printline)(const char *label, long int length)
{
    printstringc(label, length);
    Rprintf("\n");
}

void F77_SUB(printstring)(const char *label, long int length)
{
    printstringc(label, length);
}

void F77_SUB(printdouble)(double *val)
{
    Rprintf("%17.8f", *val);
}   






