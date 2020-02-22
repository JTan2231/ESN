#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "matrix.h"
#include "map.h"

int main() {
    Parent p;
    initParent(&p, 5);
    parentPrint(&p);
    parentAdd(&p, 3);
    parentAdd(&p, 6);
    parentAdd(&p, 9);
    parentPrint(&p);

    return 0;
}
