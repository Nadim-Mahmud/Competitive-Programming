/**
*   Name : all permutation generation
*   This code will print all nPr
*   While ara contains n elements
*   a is starring postion of array
*/

int ara[200],mr[200],k;

void permutation(int a,int n){
    int i;
    if(a==n+1){
        for(i=1;i<=n;i++) printf("%c",ara[i]);
        printf("\n");
        return ;
    }
    for(i=1;i<=n;i++) if(!mr[i]){
        mr[i] = 1;
        ara[a] = i;
        permutation(a+1,n);
        mr[i] = 0;
    }
}
