#include <stdio.h>
#include <time.h>
#include "generation.h"

int main() {
    srand(time(0));
    for (int i = 0; i < 10; i++)
        printf("%f\n", marsagliaNormal());
}
