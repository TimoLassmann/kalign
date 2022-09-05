
#include "stdio.h"

#include "test.h"
#include "version.h"

int add(int a, int b)
{
        return sub(a,b);//a + b;
}

int sub(int a, int b)
{
        fprintf(stdout,"%s %s\n", KALIGN_PACKAGE_NAME, KALIGN_PACKAGE_VERSION);
        return a -b;
}
