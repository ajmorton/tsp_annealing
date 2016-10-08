#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv){


    int n;
    sscanf(argv[1],"%d",&n);

    printf("%d\n", n);

    for(int i = 0; i < n; ++i){
        printf("%d %d\n", rand() % n, rand() % n);
    }


}
